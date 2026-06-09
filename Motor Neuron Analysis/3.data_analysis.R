library(pheatmap)


# 1. Estimate size factors (required for depth normalization)
dds_MN <- estimateSizeFactors(dds_MN)

# 2. Use normTransform instead of vst() for depth-normalized log2 counts
log_norm_MN <- normTransform(dds_MN)

#PCA plot
plotPCA(log_norm_MN, intgroup = "condition") +
  labs(title="PCA plot (Depth-Normalized Log2)",
       subtitle = "EAE differential expression in retinal ganglion cells",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

pca_data <- plotPCA(log_norm_MN, intgroup = "condition", returnData = TRUE)

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.7, hjust = 0.5, size = 3) +
  labs(
    title = "PCA plot with sample labels (Depth-Normalized Log2)",
    subtitle = "EAE differential expression in retinal ganglion cells",
    caption = paste0("produced on ", Sys.time())
  ) +
  theme_bw()



#Filtering based on p-value and log2 fold change 
MN_upregulated_candidates <- MN_results.df |> 
dplyr::filter(pvalue < 0.05 & log2FoldChange > 0)

#Volcano plot
ggplot(MN_results.df) +
  aes(y=-log10(pvalue), x=log2FoldChange, text = paste("Symbol:", symbol), color = pvalue < 0.05 & log2FoldChange > 0) +
  geom_point(size=2) +
  labs(title="Volcano plot",
  subtitle = "EAE diferential expression in motor neurons",
  caption=paste0("produced on ", Sys.time()),
  color = "Upregulated") +
  theme_bw()


#get top genes based on log2 fold change for heatmap
top_genes <- MN_upregulated_candidates |> 
dplyr::filter(log2FoldChange > 10) |>
arrange(desc(log2FoldChange))

#vst analysis for heatmap
vsd_MN <- vst(dds_MN, blind = FALSE)
vst_mat <- assay(vsd_MN)

vst_mat_top50 <- assay(vsd_MN)[top_genes$ensembl_id, ]

#replace ensemble ids with symbols where available, otherwise keep ensemble ids
rownames(vst_mat_top50) <- ifelse(!is.na(top_genes$symbol) & top_genes$symbol != "", 
top_genes$symbol, 
top_genes$ensembl_id)

annotation_col <- data.frame(
  Condition = colData(vsd_MN)$condition,
  row.names = colnames(vst_mat_top50)
)

ann_colors <- list(
  Condition = c(EAE = "#ff756b", Healthy = "#00BFC4")
)
pheatmap(
  vst_mat_top50,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  scale = "row",                   
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "euclidean",   
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,  
  show_colnames = FALSE,           
  main = "Top Upregulated Genes in Motor Neurons"
)
