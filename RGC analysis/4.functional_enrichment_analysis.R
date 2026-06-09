

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
gene_sets_list <- split(genesets$ensembl_gene, genesets$gs_name)

#Run ssGSEA
ssgsea_params <- ssgseaParam(expr = vst_mat, geneSets = gene_sets_list, alpha = 0.25, normalize = FALSE)
ssgsea_results <- gsva(ssgsea_params)

#get columns for EAE and Healthy samples
eae_cols <- grep("EAE", colnames(ssgsea_results))
naive_cols <- grep("Healthy", colnames(ssgsea_results)) 

# Step 1: Filter for baseline Enrichment Score > 5000 in EAE
high_scores <- rowMeans(ssgsea_results[, eae_cols]) > 5000
ssgsea_filtered <- ssgsea_results[high_scores, ]

# Step 2: Filter for statistical significance (p < 0.05)
p_values <- apply(ssgsea_filtered, 1, function(row) {
t_test <- t.test(row[eae_cols], row[naive_cols], var.equal = TRUE)
return(t_test$p.value)
})
significant_idx <- p_values < 0.05
ssgsea_sig <- ssgsea_filtered[significant_idx, ]

# Step 3: Calculate the actual enrichment delta (EAE change relative to Healthy)
eae_means <- rowMeans(ssgsea_sig[, eae_cols])
healthy_means <- rowMeans(ssgsea_sig[, naive_cols])
enrichment_delta <- eae_means - healthy_means

# Step 4: Keep only upregulated terms and sort by the biggest delta
upregulated_idx <- enrichment_delta > 0
ssgsea_up <- ssgsea_sig[upregulated_idx, ]
final_deltas <- enrichment_delta[upregulated_idx]

top_50_enriched <- names(sort(final_deltas, decreasing = TRUE))[1:50]
heatmap_mat <- ssgsea_up[top_50_enriched, ]

# Clean up prefixes across all namespaces
rownames(heatmap_mat) <- rownames(heatmap_mat) %>%
gsub("^GOBP_|^GOMF_|^GOCC_", "", .) %>%
gsub("_", " ", .) %>%
tolower()

heatmap_mat <- heatmap_mat[order(rownames(heatmap_mat)), ]

annotation_col <- data.frame(Group = colData(vsd_RGC)$condition)
rownames(annotation_col) <- colnames(heatmap_mat)

pheatmap(
heatmap_mat,
color = colorRampPalette(c("blue", "white", "red"))(50),
scale = "row", 
main = "Top 50 upregulated GO:BP Terms in RGC Neurons (ssGSEA)",
annotation_col = annotation_col,
show_colnames = TRUE,
cluster_cols = TRUE, 
cluster_rows = FALSE,
fontsize_row = 7,
  border_color = NA
  )
