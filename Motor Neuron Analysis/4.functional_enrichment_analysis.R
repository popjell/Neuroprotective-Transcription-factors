# ssGSEA analysis of MN's (Ensembl IDs & EAE/Healthy Sample Names)
library(GSVA)
library(pheatmap)
library(msigdbr)
library(tidyverse)

# 1. VSTransform the data
vsd_MN <- vst(dds_MN, blind = FALSE)
vst_mat <- assay(vsd_MN)

# 2. Get gene sets (Using ensembl_gene since your matrix uses Ensembl IDs)
genesets <- msigdbr(species = "Mus musculus", collection = "M5", db_species = "MM", subcollection = "GO:BP")
gene_sets_list <- split(genesets$ensembl_gene, genesets$gs_name)

# 3. Run ssGSEA
ssgsea_params <- ssgseaParam(expr = vst_mat, geneSets = gene_sets_list, alpha = 0.25, normalize = FALSE)
ssgsea_results <- gsva(ssgsea_params)

# 4. Map columns (Using grep since "EAE" and "Healthy" are in the column names)
eae_cols    <- grep("EAE", colnames(ssgsea_results))
naive_cols  <- grep("Healthy", colnames(ssgsea_results)) 

# 5. Statistical significance (FIX: Dropped the > 5000 filter that breaks VST data)
p_values <- apply(ssgsea_results, 1, function(row) {
  t_test <- t.test(row[eae_cols], row[naive_cols], var.equal = TRUE)
  return(t_test$p.value)
})
significant_idx <- p_values < 0.05
ssgsea_sig      <- ssgsea_results[significant_idx, , drop = FALSE]

# 6. Calculate enrichment delta
eae_means        <- rowMeans(ssgsea_sig[, eae_cols, drop = FALSE])
healthy_means    <- rowMeans(ssgsea_sig[, naive_cols, drop = FALSE])
enrichment_delta <- eae_means - healthy_means

# 7. Isolate true top 50 upregulated terms by delta size
upregulated_idx  <- enrichment_delta > 0
ssgsea_up        <- ssgsea_sig[upregulated_idx, , drop = FALSE]
final_deltas     <- enrichment_delta[upregulated_idx]

top_50_enriched  <- names(sort(final_deltas, decreasing = TRUE))[1:50]
heatmap_mat      <- ssgsea_up[top_50_enriched, , drop = FALSE]

# 8. Clean up row names for readability
rownames(heatmap_mat) <- rownames(heatmap_mat) %>%
  gsub("^GOBP_|^GOMF_|^GOCC_", "", .) %>%
  gsub("_", " ", .) %>%
  tolower()

# 9. Build annotation metadata
annotation_col <- data.frame(Group = colData(vsd_MN)$condition)
rownames(annotation_col) <- colnames(heatmap_mat)

# 10. Generate the naturally clustered heatmap
pheatmap(
  heatmap_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  scale = "row", 
  main = "Top Upregulated GO:BP Terms in Motor Neurons (ssGSEA)",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  cluster_cols = TRUE, 
  cluster_rows = TRUE,   # FIX: Rows cluster by biological pattern, not alphabet
  fontsize_row = 7,
  border_color = NA
)