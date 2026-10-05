# Compare MS and control cells within each neuronal subtype.
library(Seurat)
library(ggplot2)

x <- readRDS("data/neuron_counts.rds")
meta <- x$meta
gene_map <- x$genes
meta$subtype <- as.character(meta$cell_type)
meta$subtype[grepl("^EN-L2-3", meta$subtype)] <- "EN-L2-3"
meta$subtype[meta$subtype == "IN-PV"] <- "IN-PVALB"

# Use the same basic cell QC as the pseudobulk analysis.
keep <- meta$genes >= 500 & meta$UMIs >= 1000
meta <- meta[keep, ]
counts <- x$counts[, keep]
rm(x)
gc()

cells <- CreateSeuratObject(counts, meta.data = meta)
cells <- NormalizeData(cells, verbose = FALSE)
Idents(cells) <- "subtype"

markers <- list()
for (subtype in unique(cells$subtype)) {
  result <- FindMarkers(cells, ident.1 = "MS", ident.2 = "Control",
                        group.by = "diagnosis", subset.ident = subtype)
  result$gene_id <- rownames(result)
  result$gene <- gene_map$gene_symbol[match(result$gene_id,
                                             gene_map$ensembl_gene_id)]
  result$subtype <- subtype
  markers[[subtype]] <- result
  message(subtype, ": ", nrow(result), " genes tested")
}

de <- do.call(rbind, markers)
write.csv(de, "results/cell_level_markers.csv", row.names = FALSE)

summary <- aggregate(cbind(tested = !is.na(de$p_val),
                           nominal = !is.na(de$p_val) & de$p_val < 0.05,
                           adjusted = !is.na(de$p_val_adj) & de$p_val_adj < 0.05) ~
                       subtype, data = de, sum)
write.csv(summary, "results/cell_level_marker_counts.csv", row.names = FALSE)

p <- ggplot(summary, aes(subtype, adjusted)) +
  geom_col(fill = "#3B82F6") + coord_flip() + theme_minimal() +
  labs(title = "Schirmer cell-level markers by neuronal subtype",
       subtitle = "Seurat default Bonferroni-adjusted p < 0.05",
       x = NULL, y = "Significant genes")
ggsave("plots/cell_level_marker_counts.png", p, width = 8, height = 5, dpi = 200)
