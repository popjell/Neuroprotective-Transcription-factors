# Replot the saved UMAP with the neuronal groups used for pseudobulk.
library(ggplot2)

cells <- read.csv("results/cell_metadata.csv")
cells$subtype <- sub("^EN-L2-3.*", "EN-L2-3", cells$subtype)
cells$subtype[cells$subtype == "IN-PV"] <- "IN-PVALB"

p <- ggplot(cells, aes(UMAP_1, UMAP_2, color = subtype)) +
  geom_point(size = 0.25) +
  theme_classic() +
  labs(title = "Schirmer neuronal subtypes used for pseudobulk", color = "Subtype")
ggsave("plots/umap.png", p, width = 9, height = 6, dpi = 200)
