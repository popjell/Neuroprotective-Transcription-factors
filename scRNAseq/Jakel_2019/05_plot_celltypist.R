# Plot CellTypist labels on the existing full-atlas UMAP.
library(ggplot2)

umap <- read.csv("results/cell_metadata.csv")
labels <- read.csv("results/celltypist_annotations.csv")
umap$celltypist_label <- labels$predicted_labels[match(umap$cell, labels$cell)]

# Collapse detailed reference labels into broad cell classes.
x <- umap$celltypist_label
umap$celltypist_broad <- "Other"
umap$celltypist_broad[grepl("^(L[2-6]|InN|Chandelier)", x)] <- "Neuron"
umap$celltypist_broad[grepl("^Astro", x)] <- "Astrocyte"
umap$celltypist_broad[grepl("^Oligo", x)] <- "Oligodendrocyte"
umap$celltypist_broad[grepl("^(OPC|COP)", x)] <- "OPC/COP"
umap$celltypist_broad[grepl("^(Micro|Macro|Myeloid|T |B )", x)] <- "Immune"
umap$celltypist_broad[grepl("^Endo", x)] <- "Endothelial"
umap$celltypist_broad[grepl("^(PC |SMC|VLMC)", x)] <- "Vascular"
write.csv(umap, "results/cell_metadata.csv", row.names = FALSE)

p <- ggplot(umap, aes(UMAP_1, UMAP_2, color = celltypist_broad)) +
  geom_point(size = 0.25) + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Jakel: CellTypist annotation", color = NULL) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA))
ggsave("plots/umap_celltypist.png", p, width = 10, height = 7, dpi = 200)

# Show the CellTypist neuronal specializations separately.
umap$celltypist_subtype <- NA_character_
umap$celltypist_subtype[grepl("^L2-3", x)] <- "EX-L2-3"
umap$celltypist_subtype[grepl("^L2-4", x)] <- "EX-L2-4"
umap$celltypist_subtype[grepl("^L3-5", x)] <- "EX-L3-5"
umap$celltypist_subtype[grepl("^L5", x)] <- "EX-L5-6"
umap$celltypist_subtype[grepl("^L6 FEZF2", x)] <- "EX-L6-FEZF2"
umap$celltypist_subtype[grepl("^L6 OPRK1", x)] <- "EX-L6-OPRK1"
umap$celltypist_subtype[grepl("^InN ADARB2", x)] <- "IN-ADARB2"
umap$celltypist_subtype[grepl("^InN LAMP5", x)] <- "IN-LAMP5"
umap$celltypist_subtype[grepl("^InN PVALB", x)] <- "IN-PVALB"
umap$celltypist_subtype[grepl("^InN SST", x)] <- "IN-SST"
umap$celltypist_subtype[grepl("^InN VIP", x)] <- "IN-VIP"
umap$celltypist_subtype[grepl("^InN LHX6", x)] <- "IN-LHX6"
write.csv(umap, "results/cell_metadata.csv", row.names = FALSE)

neurons <- umap[!is.na(umap$celltypist_subtype), ]
p <- ggplot(neurons, aes(UMAP_1, UMAP_2, color = celltypist_subtype)) +
  geom_point(size = 0.4) + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Jakel: CellTypist neuronal subtypes", color = NULL) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA))
ggsave("plots/umap_celltypist_neurons.png", p, width = 10, height = 7, dpi = 200)
