# Score Rank 1 and Rank 2 gene sets with the UCell formula.
library(Matrix)
library(ggplot2)
library(UCell)

x <- readRDS("data/neuron_counts.rds")
counts <- x$counts
meta <- x$meta
gene_map <- x$genes
rm(x)

# Use the same cell QC and subtype labels as the Schirmer pseudobulk analysis.
keep <- meta$genes >= 500 & meta$UMIs >= 1000
meta <- meta[keep, ]
counts <- counts[, keep]
meta$subtype <- as.character(meta$cell_type)
meta$subtype[grepl("^EN-L2-3", meta$subtype)] <- "EN-L2-3"
meta$subtype[meta$subtype == "IN-PV"] <- "IN-PVALB"

# Use only the human genes selected for the primary rank-list analyses.
rank_genes <- read.csv("data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_ensembl", "human_symbol")])
rownames(counts) <- gene_map$ensembl_gene_id
rank1 <- intersect(rank_genes$human_ensembl[rank_genes$rank == "Rank 1"],
                    rownames(counts))
rank2 <- intersect(rank_genes$human_ensembl[rank_genes$rank == "Rank 2"],
                    rownames(counts))

# Rank genes within each cell, capping ranks at UCell's default of 1,500.
scores <- ScoreSignatures_UCell(
  matrix = counts, features = list(Rank1 = rank1, Rank2 = rank2),
  maxRank = 1500, chunk.size = 100, ncores = 1)
cell_scores <- data.frame(Rank1 = scores[, "Rank1_UCell"],
                          Rank2 = scores[, "Rank2_UCell"])

# Average cell scores within each patient and neuronal subtype.
cell_scores$patient <- meta$patient
cell_scores$diagnosis <- meta$diagnosis
cell_scores$subtype <- meta$subtype
donor_scores <- aggregate(cbind(Rank1, Rank2) ~ patient + diagnosis + subtype,
                          data = cell_scores, FUN = mean)
cell_counts <- aggregate(rep(1, nrow(meta)),
                         by = list(patient = meta$patient, subtype = meta$subtype),
                         FUN = sum)
donor_scores$cells <- cell_counts$x[match(
  paste(donor_scores$patient, donor_scores$subtype),
  paste(cell_counts$patient, cell_counts$subtype))]
write.csv(donor_scores, "results/ucell_scores_by_donor.csv", row.names = FALSE)

# Compare donor-level scores between MS and control within each subtype.
long <- reshape(donor_scores,
                varying = c("Rank1", "Rank2"), v.names = "score",
                timevar = "rank_list", times = c("Rank 1", "Rank 2"),
                direction = "long")
tests <- list()
for (rank_list in unique(long$rank_list)) {
  for (subtype in unique(long$subtype)) {
    d <- long[long$rank_list == rank_list & long$subtype == subtype, ]
    if (length(unique(d$diagnosis)) < 2) next
    test <- wilcox.test(score ~ diagnosis, data = d, exact = FALSE)
    tests[[paste(rank_list, subtype)]] <- data.frame(
      rank_list = rank_list, subtype = subtype,
      control_donors = sum(d$diagnosis == "Control"),
      ms_donors = sum(d$diagnosis == "MS"),
      control_mean = mean(d$score[d$diagnosis == "Control"]),
      ms_mean = mean(d$score[d$diagnosis == "MS"]),
      ms_minus_control = mean(d$score[d$diagnosis == "MS"]) -
        mean(d$score[d$diagnosis == "Control"]),
      p_value = test$p.value)
  }
}
tests <- do.call(rbind, tests)
tests$q_within_rank <- ave(tests$p_value, tests$rank_list,
                           FUN = function(p) p.adjust(p, method = "BH"))
tests$q_global <- p.adjust(tests$p_value, method = "BH")
write.csv(tests, "results/ucell_rank_list_tests.csv", row.names = FALSE)

plot_data <- long
plot_data$rank_list <- factor(plot_data$rank_list, levels = c("Rank 1", "Rank 2"))
p <- ggplot(plot_data, aes(diagnosis, score, color = diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.5) +
  facet_grid(rank_list ~ subtype) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Schirmer UCell scores by donor and neuronal subtype",
       subtitle = "Cell scores averaged per donor; higher score means higher within-cell ranks",
       x = NULL, y = "Mean donor UCell score")
ggsave("plots/ucell_rank_scores.png", p, width = 15, height = 7, dpi = 180)
