library(pheatmap)


# 1. Estimate size factors (required for depth normalization)
dds_RGC <- estimateSizeFactors(dds_RGC)

# 2. Use normTransform instead of vst() for depth-normalized log2 counts
log_norm_RGC <- normTransform(dds_RGC)

plotPCA(log_norm_RGC, intgroup = "condition") +
  labs(title="PCA plot (Depth-Normalized Log2)",
       subtitle = "EAE differential expression in retinal ganglion cells",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

pca_data <- plotPCA(log_norm_RGC, intgroup = "condition", returnData = TRUE)

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), vjust = -0.7, hjust = 0.5, size = 3) +
  labs(
    title = "PCA plot with sample labels (Depth-Normalized Log2)",
    subtitle = "EAE differential expression in retinal ganglion cells",
    caption = paste0("produced on ", Sys.time())
  ) +
  theme_bw()


#Filtering based on p-value and log2 fold change in paper
RGC_upregulated_candidates <- RGC_results.df |> 
dplyr::filter(pvalue < 0.05 & log2FoldChange > 0)

#Volcano plot
ggplot(RGC_results.df) +
  aes(y=-log10(pvalue), x=log2FoldChange, text = paste("Symbol:", symbol), color = pvalue < 0.05 & log2FoldChange > 0) +
  geom_point(size=2) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
  subtitle = "EAE diferential expression in RGC's",
  caption=paste0("produced on ", Sys.time()),
  color = "Upregulated") +
  theme_bw()



top_genes <- RGC_upregulated_candidates |> 
dplyr::filter(log2FoldChange > 2) |>
arrange(desc(log2FoldChange))

vsd_RGC <- vst(dds_RGC, blind = FALSE)
vst_mat_top50 <- assay(vsd_RGC)[top_genes$ensembl_id, ]

#replace ensemble ids with symbols where available, otherwise keep ensemble ids
rownames(vst_mat_top50) <- ifelse(!is.na(top_genes$symbol) & top_genes$symbol != "", 
top_genes$symbol, 
top_genes$ensembl_id)

annotation_col <- data.frame(
  Condition = colData(vsd_RGC)$condition,
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
  main = "Top 50 Upregulated Genes in RGC's"
)
