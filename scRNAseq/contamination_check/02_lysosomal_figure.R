# Figure G: violins of the lysosomal lipid-processing score and PSAP in cortical neurons, control vs MS.
# Needs data/cortex.rds, ../../Gene_Lists_Rank1_and_Rank2_Pathways.xlsx (sheet Rank2, 36 mouse symbols) and ../../rank_genes.csv (human mapping).
# Score: log-normalised expression log1p(count / nUMI * 1e4); each gene detected (>= 1 count) in these neurons is min-max scaled to 0-1 across
# neurons, then averaged per nucleus.
# Violins: points are nuclei, open circles are tissue-sample means, bars are group means of sample means. Test: Welch t on donor means
# (nuclei pooled per donor, MS371 = CR_ms3 + CR_ms4; 4 control vs 3 MS donors); stars at p < 0.05, 0.01, 0.001. Exploratory.
source("common.R")
library(patchwork)

full <- readRDS("data/cortex.rds")
cells <- full$cells[full$cells$broad == "Neuron", ]
counts <- full$counts[, cells$cell]
cat("neurons:", nrow(cells), " donors:", length(unique(cells$donor)), "\n")

# The 36 mouse symbols -> human symbols (all Rank 2 mappings, ambiguous ones included)
mouse <- as.data.frame(readxl::read_excel(rank_xlsx, sheet = "Rank2", col_names = FALSE))[, 1]
mouse <- mouse[which(mouse == "Gene") + seq_len(36)]
rank <- read.csv(rank_genes_csv)
stopifnot(length(mouse) == 36, all(mouse %in% rank$input_symbol))
rank <- rank[rank$rank == "Rank 2" & rank$input_symbol %in% mouse & !is.na(rank$human_symbol), ]
set_genes <- unique(rank$human_symbol)
x <- counts[rownames(counts) %in% set_genes, ]
lognorm <- function(m) log1p(t(t(as.matrix(m)) / cells$nUMI * norm_scale))
ln <- lognorm(x)
detected <- rowSums(x) > 0
score_genes <- data.frame(gene = rownames(x), mouse_symbols = sapply(rownames(x), function(g) paste(unique(rank$input_symbol[rank$human_symbol == g]), collapse = "/")),
                          detected_in_neurons = detected, nuclei_detecting = rowSums(x > 0))
write.csv(score_genes, res("lysosomal_score_genes.csv"), row.names = FALSE)
cat("36 mouse symbols map to", length(set_genes), "human genes;", sum(set_genes %in% rownames(counts)), "in the matrix;", sum(detected), "detected in neurons\n")

# Per-gene 0-1 scaling over neurons, then the mean per nucleus
scaled <- t(apply(ln[detected, , drop = FALSE], 1, function(v) (v - min(v)) / (max(v) - min(v))))
cells$score <- colMeans(scaled)
cells$psap <- as.numeric(lognorm(counts["PSAP", , drop = FALSE]))

# Figure G: sample means, donor means and the Welch test
sample_means <- do.call(rbind, lapply(split(cells, cells$sample_id), function(d)
  data.frame(sample_id = d$sample_id[1], donor = d$donor[1], condition = d$condition[1], n_nuclei = nrow(d), score = mean(d$score), psap = mean(d$psap))))
donor_means <- do.call(rbind, lapply(split(cells, cells$donor), function(d)
  data.frame(donor = d$donor[1], condition = d$condition[1], n_nuclei = nrow(d), score = mean(d$score), psap = mean(d$psap))))
write.csv(sample_means, res("lysosomal_sample_means.csv"), row.names = FALSE)
tests <- do.call(rbind, lapply(c(score = "score", psap = "psap"), function(v) {
  tt <- t.test(donor_means[[v]] ~ donor_means$condition)
  data.frame(variable = v, control_donors = sum(donor_means$condition == "Control"), ms_donors = sum(donor_means$condition == "MS"),
             mean_control = tt$estimate[[1]], mean_ms = tt$estimate[[2]], diff_ms_minus_control = tt$estimate[[2]] - tt$estimate[[1]],
             welch_t = tt$statistic[[1]], p = tt$p.value)
}))
write.csv(tests, res("lysosomal_figure_tests.csv"), row.names = FALSE)
print(tests, row.names = FALSE)

violin <- function(v, ylab) {
  p <- tests$p[tests$variable == v]
  stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  top <- max(cells[[v]])
  gm <- aggregate(sample_means[[v]], list(condition = sample_means$condition), mean)
  ggplot(cells, aes(condition, .data[[v]], fill = condition)) +
    geom_violin(colour = NA, alpha = 0.35, scale = "width") +
    geom_jitter(aes(colour = condition), width = 0.15, size = 0.2, alpha = 0.3) +
    geom_point(data = sample_means, aes(y = .data[[v]]), shape = 1, size = 3, stroke = 1) +
    geom_errorbar(data = gm, aes(y = x, ymin = x, ymax = x), width = 0.5, linewidth = 0.9) +
    annotate("segment", x = 1, xend = 2, y = top * 1.05, yend = top * 1.05) +
    annotate("text", x = 1.5, y = top * 1.08, label = stars, size = 5) +
    scale_fill_manual(values = c(Control = blue, MS = orange), guide = "none") +
    scale_colour_manual(values = c(Control = blue, MS = orange), guide = "none") +
    labs(x = NULL, y = ylab)
}
p_g1 <- violin("score", "Lysosomal lipid-processing score (0-1 scaled mean)") + labs(title = "Score in cortical neurons")
p_g2 <- violin("psap", "PSAP, log-normalised RNA") + labs(title = "PSAP in cortical neurons")
save_fig(p_g1 + p_g2 + plot_annotation(tag_levels = list(c("G", ""))), "lysosomal_G_violins", 8, 4.8)
