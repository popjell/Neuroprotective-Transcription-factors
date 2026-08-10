# Shared utility functions for the Neuroprotective TFs project
# Source this file from any script that needs these helpers:
#   source("utils/helpers.R")

# Build a deduplicated preranked gene vector for GSEA.
# Ranking metric = signed significance: -log10(pvalue) * sign(log2FoldChange)
# When a symbol appears multiple times, keeps the entry with the largest
# absolute value (most significant).
#   pval     - numeric vector of p-values
#   lfc      - numeric vector of log2 fold changes
#   symbols  - character vector of gene symbols
# Returns a named numeric vector sorted in decreasing order.
build_ranked <- function(pval, lfc, symbols) {
  metric <- -log10(pval) * sign(lfc)
  names(metric) <- symbols
  metric <- metric[!is.na(metric) & is.finite(metric)]
  metric <- tapply(metric, names(metric), function(x) x[which.max(abs(x))])
  sort(metric, decreasing = TRUE)
}

# Load expression data and build MN + RGC ranked gene lists.
# Returns a list with components: exp_map, mn_ranked, rgc_ranked
load_expression <- function(expression_file = "Expression_data.txt") {
  exp_map <- read.csv(expression_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  mn_ranked  <- build_ranked(exp_map$pvalue_MN, exp_map$log2FoldChange_MN, exp_map$symbol)
  rgc_ranked <- build_ranked(exp_map$pvalue_RGC, exp_map$log2FoldChange_RGC, exp_map$symbol)
  list(exp_map = exp_map, mn_ranked = mn_ranked, rgc_ranked = rgc_ranked)
}

# Load reformatted iRegulon edges and build TF -> target list.
# Returns a list with components: edges, all_tfs, tf_targets
load_tf_network <- function(edges_file = "reformatted_edges.csv") {
  edges <- read.csv(edges_file, stringsAsFactors = FALSE, check.names = FALSE)
  all_tfs <- unique(edges$`Regulator Gene`)
  tf_targets <- list()
  for (tf in all_tfs) {
    tf_targets[[tf]] <- unique(edges$`Target Gene`[edges$`Regulator Gene` == tf])
  }
  list(edges = edges, all_tfs = all_tfs, tf_targets = tf_targets)
}



read_GO <- function(go_analysis){

  all_BP <- go_analysis %>%
    dplyr::filter(source == "GO:BP") %>%
    dplyr::arrange(desc(negative_log10_of_adjusted_p_value))%>%
    dplyr::mutate(
      term_name = reorder(term_name, -negative_log10_of_adjusted_p_value)
    ) %>%
    dplyr::select(
      term_name, negative_log10_of_adjusted_p_value,intersections,term_id
    )
  return (all_BP)
}