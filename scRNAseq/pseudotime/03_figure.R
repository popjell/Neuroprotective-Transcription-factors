# Figure: module score along the author pseudotime (A line, B gene heatmap) and the nuclei on the UMAP (C pseudotime, D module score), for one gene set
library(gridExtra)
source("common.R")

nuclei <- read.csv(shared("nuclei.csv"))
umap <- read.csv(shared("umap.csv"))
scores <- read.csv(res("module_scores.csv"))
deciles <- read.csv(res("decile_summary.csv"))
gene_means <- read.csv(res("decile_gene_means.csv"))
genes <- read.csv(res("geneset_genes.csv"))
stopifnot(identical(nuclei$nucleus, umap$nucleus), identical(nuclei$nucleus, scores$nucleus))
n_module <- length(unique(gene_means$gene))
n_donors <- sum(table(nuclei$donor) >= min_cells)

# Panel A: module score per pseudotime decile
pA <- ggplot(deciles, aes(mean_pseudotime, mean_module)) +
  geom_ribbon(aes(ymin = mean_module - se_module, ymax = mean_module + se_module), fill = orange, alpha = 0.25) +
  geom_line(color = orange, linewidth = 1) + geom_point(color = orange, size = 2.4) +
  labs(tag = "A", title = "Module score along pseudotime", x = "Pseudotime", y = "Module score (UCell)")

# Panel B: module genes x decile, row z-score, genes ordered by their correlation with pseudotime
wide <- tapply(gene_means$mean_expression, list(gene_means$gene, gene_means$decile), mean)
z <- t(scale(t(wide)))
gene_order <- rownames(wide)[order(cor(t(wide), 1:10, method = "spearman"), decreasing = TRUE)]
heat <- data.frame(gene = factor(rep(rownames(z), 10), levels = rev(gene_order)), decile = rep(1:10, each = nrow(z)), z = as.vector(z))
pB <- ggplot(heat, aes(decile, gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = blue, mid = "white", high = orange, midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "Row\nz-score") +
  scale_x_continuous(breaks = 1:10, expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(tag = "B", title = "Module genes along pseudotime", x = "Pseudotime decile (early to late)", y = NULL) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic", size = if (n_module > 25) 7 else 9))

# Panels C and D: UMAP in hexagons, mean per hexagon (hexagons with fewer than 3 nuclei are left out)
umap$pseudotime <- nuclei$author_pseudotime
umap$module <- pmin(scores$module_ucell, quantile(scores$module_ucell, 0.99))
mean_min3 <- function(z) if (length(z) >= 3) mean(z) else NA_real_
umap_panel <- function(value, label, tag, title) {
  ggplot(umap, aes(UMAP1, UMAP2, z = .data[[value]])) +
    stat_summary_hex(fun = mean_min3, bins = 32, na.rm = TRUE) +
    scale_fill_viridis_c(option = "D", name = label, guide = guide_colourbar(barheight = unit(2.6, "cm"), barwidth = unit(0.3, "cm"))) +
    labs(tag = tag, title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}
pC <- umap_panel("pseudotime", "Mean\npseudotime", "C", "UMAP: mean pseudotime per hexagon")
pD <- umap_panel("module", "Mean\nmodule\nscore", "D", "UMAP: mean module score per hexagon")

caption <- wrap(paste0("EN-L2-3 nuclei from Schirmer 2019 (", nrow(nuclei), " nuclei, ", n_donors, " donors), ordered by the pseudotime of the authors. ",
                        "Module: ", n_module, " of ", nrow(genes), " genes of ", geneset_arg, " (genes absent from the data or detected in fewer than ", 100 * min_detect,
                        "% of nuclei are left out; see results/", geneset_name, "/geneset_genes.csv). ",
                        "A, module score (UCell), mean +/- s.e.m. per pseudotime decile (about ", round(nrow(nuclei) / 10), " nuclei each). ",
                        "B, mean log-normalised expression per decile as a row z-score, genes ordered by their correlation with pseudotime. ",
                        "C, D, UMAP of the nuclei (Monocle 3, aligned across sequencing batches) in hexagons, hexagons with fewer than 3 nuclei omitted; ",
                        "colour is the mean pseudotime (C) or the mean module score (D)."))
top_height <- max(4.4, 0.25 * n_module + 1.5)
p <- arrangeGrob(pA, pB, pC, pD, layout_matrix = rbind(c(1, 2), c(3, 4)), widths = c(1, 1.1), heights = c(top_height, 4.8),
                 bottom = grid::textGrob(caption, x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 8)))
save_fig(p, "module_pseudotime", 12, top_height + 4.8 + 0.6)
