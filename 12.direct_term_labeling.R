library(tidyr)
library(reticulate)
use_python("/home/jathavp/lab/Neuroprotective-Transcription-factors/.venv/bin/python3.11")
#Convert GO list into gene specific with GO lists afterwards


#GO extraction
if(!exists("intersected_combined_df")){
  intersected_combined_df <- read.table("Expression_data.txt", sep = "\t", header = TRUE)
}
if(!exists("all_BP")){  
  source("utils/helpers.R")
  all_BP <- read_GO(go_analysis)
}
id_conversion_intersected <- intersected_combined_df %>%
  dplyr::select(
    ensembl_id, symbol
  )

gprofiler_long <- separate_rows(all_BP,intersections, sep = ",\\s*|\\s+") %>%
  dplyr::left_join(id_conversion_intersected, by = c("intersections" = "ensembl_id")) %>%
  dplyr::select(
    term_name, symbol
  )
terms_in_symbol <- gprofiler_long %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    go_terms = paste(term_name, collapse = ", ")
  ) %>%
  dplyr::mutate(
    num_terms = lengths(strsplit(go_terms, ",\\s*"))
  )


#metabolism and cell d
py_run_string("
import pandas as pd
input_df = r.all_BP
")
py_run_file("go_embeddings/GO_graph_traversal.py")
embeddings <- py$traversal_df %>%
  dplyr::mutate(
    metabolism = ifelse(grepl("metabolism", category) & !grepl("exclude", category), TRUE, FALSE),
    cell_killing = ifelse(grepl("cell_death", category) & !grepl("exclude", category), TRUE, FALSE)
  ) %>%
  dplyr::select(
    term_name, metabolism, cell_killing, category
  )

metabolism_terms <- embeddings %>%
  dplyr::filter(metabolism == TRUE) %>%
  dplyr::pull(term_name)
cell_death_terms <- embeddings %>%
  dplyr::filter(cell_killing == TRUE) %>%
  dplyr::pull(term_name)

metabolism_genes <- terms_in_symbol %>%
  dplyr::mutate(
    split_terms = strsplit(go_terms, ",\\s*")
  ) %>%
  dplyr::mutate(
    metabolism = sapply(split_terms, function(x) any(x %in% metabolism_terms)),
    metabolism_terms = sapply(split_terms, function(x) paste(intersect(x, metabolism_terms), collapse = ", ")),
    num_metabolism_terms = sapply(split_terms, function(x) sum(x %in% metabolism_terms))
  ) %>%
  dplyr::filter(metabolism == TRUE) %>%
  dplyr::select(-split_terms)
cell_death_genes <- terms_in_symbol %>%
  dplyr::mutate(
    split_terms = strsplit(go_terms, ",\\s*")
  ) %>%
  dplyr::mutate(
    cell_killing = sapply(split_terms, function(x) any(x %in% cell_death_terms)),
    cell_killing_terms = sapply(split_terms, function(x) paste(intersect(x, cell_death_terms), collapse = ", ")),
    num_cell_killing_terms = sapply(split_terms, function(x) sum(x %in% cell_death_terms))
  ) %>%
  dplyr::filter(cell_killing == TRUE) %>%
  dplyr::select(-split_terms)

metabolism_genelist <- metabolism_genes%>%
  dplyr::pull(symbol)
cell_death_genelist <- cell_death_genes%>%
  dplyr::pull(symbol)


library(readxl)
#call in the LLM gene annotation file and make clean metabolism and cell death df's with the information available
llm_gene_annotation_filepath <- "LLM_gene_annotations_shakour.xlsx"

sheet_names <- excel_sheets(llm_gene_annotation_filepath)
all_sheets_list <- lapply(sheet_names, function(sheet) {
  read_excel(llm_gene_annotation_filepath, sheet = sheet)
})
names(all_sheets_list) <- sheet_names


llm_metabolism_sheet <- as.data.frame(all_sheets_list[["Metabolism"]]) 
names(llm_metabolism_sheet) <- as.character(llm_metabolism_sheet[2,])
llm_metabolism_sheet <- llm_metabolism_sheet[-c(1, 2), ]%>%
  dplyr::rename(symbol = `Mouse symbol`)
rownames(llm_metabolism_sheet) <- NULL

llm_celldeath_sheet <- as.data.frame(all_sheets_list[["Cell Death"]])
names(llm_celldeath_sheet) <- as.character(llm_celldeath_sheet[2,])
llm_celldeath_sheet <- llm_celldeath_sheet[-c(1, 2), ]%>%
  dplyr::rename(symbol = `Mouse symbol`)
rownames(llm_celldeath_sheet) <- NULL



llm_metabolism_genes <- llm_metabolism_sheet %>%
  dplyr::pull('symbol')
llm_cell_death_genes <- llm_celldeath_sheet %>%
  dplyr::pull('symbol')



#comparsions

intersect(llm_metabolism_genes, metabolism_genelist)
intersect(llm_cell_death_genes, cell_death_genelist)
print("--- METABOLISM DIFFERENCES ---")
print("GO-Only (Over-inclusive broad terms):")
setdiff(llm_metabolism_genes, metabolism_genelist)
print("LLM-Only (Specific/Peripheral terms):")
setdiff(metabolism_genelist, llm_metabolism_genes)
print("--- CELL DEATH DIFFERENCES ---")
print("GO-Only (Over-inclusive broad terms):")
setdiff(llm_cell_death_genes, cell_death_genelist)
print("LLM-Only (Specific/Peripheral terms):")
setdiff(cell_death_genelist, llm_cell_death_genes)

data.frame(
  Total_LLM = length(llm_metabolism_genes),
  Total_Meta = length(metabolism_genelist),
  Shared = length(intersect(llm_metabolism_genes, metabolism_genelist)),
  Only_in_LLM = length(setdiff(llm_metabolism_genes, metabolism_genelist)),
  Only_in_Meta = length(setdiff(metabolism_genelist, llm_metabolism_genes))
)


data.frame(
  Total_LLM = length(llm_cell_death_genes),
  Total_CD = length(cell_death_genelist),
  Shared = length(intersect(llm_cell_death_genes, cell_death_genelist)),
  Only_in_LLM = length(setdiff(llm_cell_death_genes, cell_death_genelist)),
  Only_in_CD = length(setdiff(cell_death_genelist, llm_cell_death_genes))
)


metabolism_comparasion_table <- metabolism_genes %>%
  full_join(
    llm_metabolism_sheet,
    by = "symbol",
    suffix = c(".GO", ".LLM")
  ) %>%
  dplyr::rename(
    metabolism_GO = metabolism
  ) %>%
  dplyr::mutate(
    metabolism_llm = ifelse(!is.na(`metabolism_GO`), TRUE, FALSE)
  ) %>%
  dplyr::select(-go_terms) %>%
  dplyr::left_join(
    terms_in_symbol %>%
      dplyr::select(symbol, go_terms),
    by = "symbol"
  )

cell_death_comparasion_table <- cell_death_genes %>%
  full_join(
    llm_celldeath_sheet,
    by = "symbol",
    suffix = c(".GO", ".LLM")
  ) %>%
  dplyr::rename(
    cell_killing_GO = cell_killing
  ) %>%
  dplyr::mutate(
    cell_killing_llm = ifelse(!is.na(`cell_killing_GO`), TRUE, FALSE)
  ) %>%
  dplyr::select(-go_terms) %>%
  dplyr::left_join(
    terms_in_symbol %>%
      dplyr::select(symbol, go_terms),
    by = "symbol"
  )
 
write(metabolism_genelist, "filtered_genelists/GO_metabolism_genes.txt")
write(cell_death_genelist, "filtered_genelists/GO_cell_death_genes.txt")

write(llm_metabolism_genes, "filtered_genelists/LLM_metabolism_genes.txt")
write(llm_cell_death_genes, "filtered_genelists/LLM_cell_death_genes.txt")