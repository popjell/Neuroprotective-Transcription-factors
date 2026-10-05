# Compare MS and control nuclei within each author-labelled neuronal subtype.
library(Seurat)
library(ggplot2)

x <- readRDS("data/new_regions.rds")
meta <- x$meta

is_neuron <- grepl("(^L[2-6]|-IN$|-EN$|^PC$|-NEU$)", meta$regional_subtype)
meta <- meta[is_neuron, ]
counts <- x$counts[, rownames(meta)]
rm(x)
gc()

cells <- CreateSeuratObject(counts, meta.data = meta)
cells <- NormalizeData(cells, verbose = FALSE)
Idents(cells) <- "regional_subtype"

markers <- list()
for (subtype in unique(cells$regional_subtype)) {
  result <- FindMarkers(cells, ident.1 = "MS", ident.2 = "Control",
                        group.by = "condition", subset.ident = subtype)
  result$gene <- rownames(result)
  result$subtype <- subtype
  result$region <- unique(meta$cns_region[meta$regional_subtype == subtype])
  markers[[subtype]] <- result
  message(subtype, ": ", nrow(result), " genes tested")
}

de <- do.call(rbind, markers)
write.csv(de, "results/cell_level_markers.csv", row.names = FALSE)

summary <- aggregate(cbind(tested = !is.na(de$p_val),
                           nominal = !is.na(de$p_val) & de$p_val < 0.05,
                           adjusted = !is.na(de$p_val_adj) & de$p_val_adj < 0.05) ~
                       subtype, data = de, sum)
summary$region <- de$region[match(summary$subtype, de$subtype)]
write.csv(summary, "results/cell_level_marker_counts.csv", row.names = FALSE)

p <- ggplot(summary, aes(subtype, adjusted)) +
  geom_col(fill = "#3B82F6") + coord_flip() + theme_minimal() +
  labs(title = "Trobisch cell-level markers by neuronal subtype",
       subtitle = "Seurat default Bonferroni-adjusted p < 0.05",
       x = NULL, y = "Significant genes")
ggsave("plots/cell_level_marker_counts.png", p, width = 8, height = 5, dpi = 200)
