library(dplyr)
library(tidyr)

# 1. Load Files
df <- read.csv("intersected_results.csv", skip = 11, stringsAsFactors = FALSE, check.names = FALSE)
all_upregulated_genes <- readLines("intersected_upregulated_genes.txt")
expression_data <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

expression_data <- expression_data %>%
  mutate(
    is_sig_RGC = !is.na(pvalue_RGC) & !is.na(log2FoldChange_RGC) & pvalue_RGC < 0.05 & abs(log2FoldChange_RGC) > 0.58,
    is_sig_MN  = !is.na(pvalue_MN)  & !is.na(log2FoldChange_MN)  & pvalue_MN  < 0.05 & abs(log2FoldChange_MN)  > 0.58,

    Specificity = case_when(
      is_sig_RGC & is_sig_MN & (sign(log2FoldChange_RGC) != sign(log2FoldChange_MN)) ~ "Opposite",
      is_sig_RGC & is_sig_MN ~ "Both",
      is_sig_RGC ~ "RGC_Only",
      is_sig_MN  ~ "MN_Only"
    ),

    `EAE_Combined_Score` = case_when(
      Specificity == "Opposite" ~ 0,
      Specificity == "Both"     ~ (log2FoldChange_RGC + log2FoldChange_MN) / 2,
      Specificity == "RGC_Only" ~ log2FoldChange_RGC,
      Specificity == "MN_Only"  ~ log2FoldChange_MN
    )
)
# 2. Build the Complete Unique TF-to-Target Edge List
edges <- df %>%
  filter(!is.na(`Transcription factor`) & `Transcription factor` != "") %>%
  select(`Transcription factor`, `Target genes`, `Motif id`, `NES`) %>%
  separate_rows(`Transcription factor`, sep = ",\\s*") %>%
  separate_rows(`Target genes`, sep = ",\\s*") %>%
  mutate(
    `Regulator Gene` = trimws(`Transcription factor`),
    `Target Gene` = trimws(`Target genes`),
    interaction = paste0("regulates via ", trimws(`Motif id`)),
    Assembly = "mm9"
  ) %>%
  filter(`Regulator Gene` != "" & `Target Gene` != "") %>%
  select(`Regulator Gene`, `Target Gene`, interaction, `Motif id`, NES, Assembly) %>%
  distinct(`Regulator Gene`, `Target Gene`, .keep_all = TRUE)

# 3. Build Node Table (Explicitly labeling the Targets)
all_tfs <- unique(edges$`Regulator Gene`)
background_genes <- unique(trimws(all_upregulated_genes))

# Merge TF nodes with expression data and filter for significance in at least one dataset
tf_nodes <- data.frame(
  id = all_tfs, 
  `Regulatory function` = "Regulator", 
  check.names = FALSE
) %>%
  left_join(
    expression_data,
    by = c("id" = "symbol")
  ) %>%
filter((is_sig_RGC | is_sig_MN) & Specificity != "Opposite") 
  

# NEW: Clean up the edges list so it only contains the filtered TFs
edges <- edges %>%
  filter(`Regulator Gene` %in% tf_nodes$id)

# Merge gene nodes with expression data
gene_nodes <- data.frame(
  id = background_genes, 
  `Regulatory function` = "Regulated",
  check.names = FALSE
) %>%
  left_join(
    expression_data,
    by = c("id" = "symbol")
  ) %>%
  filter(!is.na(ensembl_id))

# Combine. TFs remain "Regulator" if they happen to be in the background list
node_attributes <- bind_rows(tf_nodes, gene_nodes) %>%
  distinct(id, .keep_all = TRUE)

# Export as tab-separated files matching your file setup
write.table(edges, "iregulon_edges.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(node_attributes, "metadata.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# write regulators (one symbol per line)
writeLines(sort(unique(as.character(tf_nodes$id))), "regulators.txt")
