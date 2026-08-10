library(dplyr)
library(tidyr)

reformat_networks <- function(iregulon_results, gene_list_file, expression_file = "Expression_data.txt"){
  
  sep <- ifelse(grepl("\\.tsv$", iregulon_results, ignore.case = TRUE), "\t", ",")
  df <- read.csv(iregulon_results, skip = 11, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  all_upregulated_genes <- readLines(gene_list_file)
  expression_data <<- read.csv(expression_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

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

  edges <- df %>%
    dplyr::filter(!is.na(`Transcription factor`) & `Transcription factor` != "") %>%
    dplyr::select(`Transcription factor`, `Target genes`, `Motif id`, `NES`) %>%
    separate_rows(`Transcription factor`, sep = ",\\s*") %>%
    separate_rows(`Target genes`, sep = ",\\s*") %>%
    mutate(
      `Regulator Gene` = trimws(`Transcription factor`),
      `Target Gene` = trimws(`Target genes`),
      interaction = paste0("regulates via ", trimws(`Motif id`)),
      Assembly = "mm9"
    ) %>%
    dplyr::filter(`Regulator Gene` != "" & `Target Gene` != "") %>%
    dplyr::select(`Regulator Gene`, `Target Gene`, interaction, `Motif id`, NES, Assembly) %>%
    distinct(`Regulator Gene`, `Target Gene`, .keep_all = TRUE)

  all_tfs <- unique(edges$`Regulator Gene`)
  background_genes <- unique(trimws(all_upregulated_genes))

  tf_nodes <- data.frame(
    id = all_tfs, 
    `Regulatory function` = "Regulator", 
    check.names = FALSE
  ) %>%
    left_join(expression_data, by = c("id" = "symbol")) %>%
    dplyr::filter((is_sig_RGC | is_sig_MN) & Specificity != "Opposite") 
    
  edges <- edges %>%
    dplyr::filter(`Regulator Gene` %in% tf_nodes$id)

# Combine background genes (Ensembl IDs) with expression data
  gene_nodes <- data.frame(
    id_ens = background_genes,
    `Regulatory function` = "Regulated",
    check.names = FALSE
  )

  if (any(grepl("^ENSMUSG", background_genes))) {
    # Input list contains Ensembl IDs: join on ensembl_id, keep the symbol
    gene_nodes <- left_join(gene_nodes, expression_data, by = c("id_ens" = "ensembl_id")) %>%
      mutate(id = symbol)
  } else {
    # Input list contains gene symbols: join on symbol, keep id_ens
    gene_nodes <- left_join(gene_nodes, expression_data, by = c("id_ens" = "symbol")) %>%
      mutate(id = id_ens)
  }

  gene_nodes <- gene_nodes %>%
    dplyr::filter(!is.na(id)) %>%
    dplyr::select(-id_ens)
  node_attributes <- bind_rows(tf_nodes, gene_nodes) %>%
    distinct(id, .keep_all = TRUE)


  # Return lists of data frames to use dynamically
  return(list(edges = edges, nodes = node_attributes, tfs = tf_nodes, genes = gene_nodes))
}

results <- "iregulon_results.csv"
list <- "intersected_upregulated_genes.txt"
reformat_result <- reformat_networks(results, list)

write.csv(reformat_result$edges, "reformatted_edges.csv", row.names = FALSE)
write.csv(reformat_result$nodes, "reformatted_nodes.csv", row.names = FALSE)
writeLines(sort(unique(as.character(reformat_result$tfs$id))), "regulators.txt")
