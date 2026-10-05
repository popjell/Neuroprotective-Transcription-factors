# Plot the authors' full integrated atlas annotations on their published UMAP.
library(ggplot2)
library(ggrepel)

cells <- read.delim("data/meta.tsv")
coords <- read.delim(gzfile("data/UMAP.coords.tsv.gz"), header = FALSE,
                     col.names = c("cellId", "UMAP_1", "UMAP_2"))
cells$UMAP_1 <- coords$UMAP_1[match(cells$cellId, coords$cellId)]
cells$UMAP_2 <- coords$UMAP_2[match(cells$cellId, coords$cellId)]

labels <- aggregate(cbind(UMAP_1, UMAP_2) ~ cell_subtypes_intg,
                    data = cells, median)
p <- ggplot(cells, aes(UMAP_1, UMAP_2, color = cell_subtypes_intg)) +
  geom_point(size = 0.15, alpha = 0.75) +
  geom_label_repel(data = labels, aes(label = cell_subtypes_intg),
                   max.overlaps = Inf, size = 2.7, show.legend = FALSE) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Trobisch full atlas with author annotations", color = "Cell type")
ggsave("plots/full_atlas_umap.png", p, width = 11, height = 8, dpi = 200)
