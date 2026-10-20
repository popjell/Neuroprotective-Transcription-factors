# Figure: z of GENE's association with each program against matched genes, with donor-bootstrap interval and BH q
source("common.R")

tab <- read.csv(res("program_specificity.csv"))
tab <- tab[order(-tab$z, tab$program), ]
cells <- readRDS("data/cortex_neurons.rds")$cells
donor_condition <- table(unique(cells[, c("donor", "condition")])$condition)
n_null_genes <- tab$n_null[1]

tab$program <- factor(tab$program, levels = rev(tab$program))
tab$label <- paste("q =", ifelse(tab$q_BH < 0.001, "<0.001", formatC(tab$q_BH, format = "f", digits = 3)))
x_text <- ceiling(max(tab$z_boot_hi)) + 0.4
x_min <- min(-2, floor(min(tab$z_boot_lo)))

p <- ggplot(tab, aes(z, program)) +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_vline(xintercept = 1.96, linetype = "dashed") +
  geom_errorbar(aes(xmin = z_boot_lo, xmax = z_boot_hi), orientation = "y", width = 0.25) +
  geom_point(shape = 21, size = 3.5, fill = orange) +
  geom_text(aes(x = x_text, label = label), hjust = 0, size = 3.5) +
  scale_y_discrete(labels = setNames(sprintf("%s (n = %d)", tab$program, tab$n_subtypes), tab$program)) +
  scale_x_continuous(limits = c(x_min, x_text + 1.4), breaks = seq(x_min, ceiling(max(tab$z_boot_hi)), 1)) +
  labs(title = paste(gene, "vs expression-matched genes: z by", set_label),
       subtitle = paste0("Dot = z; bar = donor-bootstrap 95% interval; q = BH over ", nrow(tab), " programs; dashed line z = 1.96."),
       x = paste("z against", n_null_genes, "expression-matched genes"), y = NULL,
       caption = wrap("Trobisch cortex neurons (Schirmer 2019 tissue),", length(unique(cells$subtype)), "subtypes,", donor_condition[["Control"]], "control and", donor_condition[["MS"]],
                      "MS donors; donors are the unit and cells are never pooled across donors. For each program (n = subtypes tested), the statistic is the median across subtypes of the",
                      gene, "partial Spearman r with the UCell program score (adjusting for nUMI and genes detected; donors combined by random effects). z compares it with", n_null_genes,
                      "genes outside every gene set that are closest to", gene, "in mean expression and fraction of cells detected, the same genes in every subtype.",
                      "Bars: 95% interval over", n_boot, "resamples of the donors, matched genes recomputed in each. q: Benjamini-Hochberg over the two-sided empirical p (share of matched genes with |r| at least that of",
                      paste0(gene, "; smallest possible p ", signif(1 / (n_null_genes + 1), 2), "). Subtypes share donors and programs share genes, so q is approximate.",
                             " Correlational; no direction or effect is claimed. Condition is confounded with batch and brain bank."))) +
  theme(panel.grid.minor = element_blank()) +
  caption_theme
save_fig(p, "program_specificity", 9, 3 + 0.4 * nrow(tab))
