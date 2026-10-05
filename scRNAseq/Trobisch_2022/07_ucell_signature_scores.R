# Score Rank 1 and Rank 2 genes in Trobisch's new neuronal samples.
library(Matrix)
library(UCell)
library(ggplot2)

x <- readRDS("data/new_regions.rds")
counts <- x$counts
meta <- x$meta
genes <- readLines("data/gene_names.txt")
if (length(genes) != nrow(counts)) stop("Gene names do not match the count matrix")
rownames(counts) <- genes

# Use only the six author-labelled neuronal populations in the new regions.
neurons <- c("GRAN-EN", "FRMPD4-IN", "KCNC2-IN", "PC",
             "DORSAL-NEU", "VENTRAL-NEU")
keep <- meta$regional_subtype %in% neurons
meta <- meta[keep, ]
counts <- counts[, keep]
if (!all(colSums(counts) == meta$nCount_RNA)) stop("Counts are not raw RNA counts")

# Use the same primary Rank 1 and Rank 2 genes as the Schirmer analysis.
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
rank1 <- intersect(rank_genes$human_symbol[rank_genes$rank == "Rank 1"], genes)
rank2 <- intersect(rank_genes$human_symbol[rank_genes$rank == "Rank 2"], genes)

# Calculate within-cell expression ranks using UCell's default 1,500-gene cap.
scores <- ScoreSignatures_UCell(
  matrix = counts, features = list(Rank1 = rank1, Rank2 = rank2),
  maxRank = 1500, chunk.size = 100, ncores = 1)
cell_scores <- data.frame(
  Rank1 = scores[, "Rank1_UCell"],
  Rank2 = scores[, "Rank2_UCell"],
  cell = rownames(scores))
meta$cell <- rownames(meta)
cell_scores <- merge(cell_scores, meta[, c("cell", "sample_id", "condition",
                                           "cns_region", "regional_subtype")],
                     by = "cell", sort = FALSE)

# Average cells so each donor contributes one score per neuronal subtype.
donor_scores <- aggregate(cbind(Rank1, Rank2) ~ sample_id + condition + cns_region +
                            regional_subtype,
                          data = cell_scores, FUN = mean)
names(donor_scores)[names(donor_scores) == "regional_subtype"] <- "subtype"
cell_counts <- aggregate(rep(1, nrow(cell_scores)),
                         by = list(sample_id = cell_scores$sample_id,
                                   cns_region = cell_scores$cns_region,
                                   subtype = cell_scores$regional_subtype), FUN = sum)
donor_scores$cells <- cell_counts$x[match(
  paste(donor_scores$sample_id, donor_scores$cns_region, donor_scores$subtype),
  paste(cell_counts$sample_id, cell_counts$cns_region, cell_counts$subtype))]
write.csv(donor_scores, "results/ucell_scores_by_donor.csv", row.names = FALSE)

# Test MS versus control using donors, not cells, as the replicates.
long <- reshape(donor_scores,
                varying = c("Rank1", "Rank2"), v.names = "score",
                timevar = "rank_list", times = c("Rank 1", "Rank 2"),
                direction = "long")
tests <- list()
for (rank_list in unique(long$rank_list)) {
  for (subtype in unique(long$subtype)) {
    d <- long[long$rank_list == rank_list & long$subtype == subtype, ]
    if (length(unique(d$condition)) < 2) next
    test <- wilcox.test(score ~ condition, data = d, exact = FALSE)
    tests[[paste(rank_list, subtype)]] <- data.frame(
      rank_list = rank_list, subtype = subtype,
      region = unique(d$cns_region), control_donors = sum(d$condition == "Control"),
      ms_donors = sum(d$condition == "MS"),
      control_mean = mean(d$score[d$condition == "Control"]),
      ms_mean = mean(d$score[d$condition == "MS"]),
      ms_minus_control = mean(d$score[d$condition == "MS"]) -
        mean(d$score[d$condition == "Control"]),
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
p <- ggplot(plot_data, aes(condition, score, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.8, size = 2) +
  facet_grid(rank_list ~ subtype, scales = "free_x") +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "Trobisch UCell scores by donor and neuronal subtype",
       subtitle = "New cerebellum and spinal-cord samples; each point is one donor",
       x = NULL, y = "Mean donor UCell score")
ggsave("plots/ucell_rank_scores.png", p, width = 13, height = 6, dpi = 180)
