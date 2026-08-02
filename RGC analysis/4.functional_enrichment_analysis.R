

#ssGSEA analysis of RGC's
library(GSVA)
library(pheatmap)
library(GSEABase)
library(msigdbr)
library(ExperimentHub)

#VSTransform the data for ssGSEA input
vsd_RGC <- vst(dds_RGC, blind = FALSE)
vst_mat <- assay(vsd_RGC)

#Get gene sets from MSigDB for mouse GO:BP
genesets <- msigdbr(species = "Mus musculus", collection = "M5", db_species = "MM", subcollection = "GO:BP")
gene_sets_list <- split(genesets$gene_symbol, genesets$gs_name)
# Run ssGSEA directly on your VST matrix
ssgsea_params <- ssgseaParam(expr = vst_mat, geneSets = gene_sets_list, alpha = 0.25, normalize = FALSE)
ssgsea_results <- gsva(ssgsea_params)

# Map sample groups using metadata indices
eae_cols    <- which(coldata_RGC$condition == "EAE")
naive_cols  <- which(coldata_RGC$condition == "Healthy")

# Step 1: Calculate statistical significance across ALL pathways (skip the >5000 filter)
p_values <- apply(ssgsea_results, 1, function(row) {
  t_test <- t.test(row[eae_cols], row[naive_cols], var.equal = TRUE)
  return(t_test$p.value)
})
significant_idx <- p_values < 0.05
ssgsea_sig      <- ssgsea_results[significant_idx, , drop = FALSE]

# Step 2: Calculate the actual enrichment delta (EAE change relative to Healthy)
eae_means        <- rowMeans(ssgsea_sig[, eae_cols, drop = FALSE])
healthy_means    <- rowMeans(ssgsea_sig[, naive_cols, drop = FALSE])
enrichment_delta <- eae_means - healthy_means

# Step 3: Filter for upregulated terms and isolate a TRUE top 50 by delta size
upregulated_idx  <- enrichment_delta > 0
ssgsea_up        <- ssgsea_sig[upregulated_idx, , drop = FALSE]
final_deltas     <- enrichment_delta[upregulated_idx]

top_50_enriched  <- names(sort(final_deltas, decreasing = TRUE))[1:50]
heatmap_mat      <- ssgsea_up[top_50_enriched, , drop = FALSE]

# Clean up row names for readability
rownames(heatmap_mat) <- rownames(heatmap_mat) %>%
  gsub("^GOBP_|^GOMF_|^GOCC_", "", .) %>%
  gsub("_", " ", .) %>%
  tolower()

# Construct annotation columns
annotation_col <- data.frame(Group = colData(vsd_RGC)$condition)
rownames(annotation_col) <- colnames(heatmap_mat)

# Plot the real, naturally clustered top 50 pathways
pheatmap(
  heatmap_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  scale = "row", 
  main = "Top Upregulated GO:BP Terms in RGCs (ssGSEA)",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  cluster_cols = TRUE, 
  cluster_rows = TRUE, # Let them cluster by signature, not alphabet
  fontsize_row = 7,
  border_color = NA
)