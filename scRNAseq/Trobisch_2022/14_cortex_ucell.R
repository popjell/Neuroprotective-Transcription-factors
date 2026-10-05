# Score both rank lists in the integrated cortical labels.
# These cortical donors overlap Schirmer; this is not independent replication.
library(Matrix)
library(UCell)
library(ggplot2)

x <- readRDS("data/cortex_neurons.rds")
counts <- x$counts
meta <- x$meta
stopifnot(all(colSums(counts) == meta$nCount_RNA))
rm(x)
gc()
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
sets <- lapply(sets, intersect, y = rownames(counts))

cell_scores <- ScoreSignatures_UCell(
  matrix = counts, features = sets, maxRank = 1500,
  chunk.size = 100, ncores = 1)
cell_scores <- as.data.frame(cell_scores)
names(cell_scores) <- names(sets)
cell_scores$cellId <- rownames(cell_scores)
cell_scores <- merge(cell_scores, meta[, c("cellId", "sample_id", "condition",
                                           "cell_subtypes_intg")],
                     by = "cellId", sort = FALSE)
names(cell_scores)[names(cell_scores) == "cell_subtypes_intg"] <- "subtype"

donors <- aggregate(cbind(`Rank 1`, `Rank 2`) ~ sample_id + condition + subtype,
                    cell_scores, mean)
cells <- aggregate(cellId ~ sample_id + subtype, cell_scores, length)
donors$cells <- cells$cellId[match(paste(donors$sample_id, donors$subtype),
                                  paste(cells$sample_id, cells$subtype))]
long <- reshape(donors, varying = c("Rank 1", "Rank 2"), v.names = "score",
                timevar = "rank", times = c("Rank 1", "Rank 2"), direction = "long")
rownames(long) <- NULL
write.csv(long, "results/cortex_integrated_ucell_donors.csv", row.names = FALSE)

tests <- list()
for (subtype in sort(unique(long$subtype))) {
  for (rank_name in names(sets)) {
    d <- long[long$subtype == subtype & long$rank == rank_name, ]
    if (sum(d$condition == "Control") < 2 || sum(d$condition == "MS") < 2) next
    model <- lm(score ~ factor(condition, levels = c("Control", "MS")), d)
    tests[[paste(subtype, rank_name)]] <- data.frame(
      subtype = subtype, rank = rank_name,
      control_donors = sum(d$condition == "Control"),
      ms_donors = sum(d$condition == "MS"),
      ms_minus_control = unname(coef(model)[2]),
      p_value = summary(model)$coefficients[2, 4])
  }
}
tests <- do.call(rbind, tests)
tests$q_global <- p.adjust(tests$p_value, method = "BH")
write.csv(tests, "results/cortex_integrated_ucell_tests.csv", row.names = FALSE)

# Give equal weight to labels represented in every cortical donor.
coverage <- aggregate(sample_id ~ subtype, long, function(z) length(unique(z)))
complete <- coverage$subtype[coverage$sample_id == length(unique(long$sample_id))]
balanced <- long[long$subtype %in% complete, ]
broad <- aggregate(score ~ sample_id + condition + rank, balanced, mean)
broad_tests <- list()
for (rank_name in names(sets)) {
  d <- broad[broad$rank == rank_name, ]
  model <- lm(score ~ factor(condition, levels = c("Control", "MS")), d)
  broad_tests[[rank_name]] <- data.frame(
    rank = rank_name, donors = nrow(d), labels = length(complete),
    ms_minus_control = unname(coef(model)[2]),
    p_value = summary(model)$coefficients[2, 4])
}
broad_tests <- do.call(rbind, broad_tests)
broad_tests$q_two_lists <- p.adjust(broad_tests$p_value, method = "BH")
write.csv(broad_tests, "results/cortex_integrated_broad_ucell.csv", row.names = FALSE)
write.csv(broad, "results/cortex_integrated_broad_ucell_donors.csv", row.names = FALSE)

p <- ggplot(broad, aes(condition, score, color = condition)) +
  geom_point(position = position_jitter(width = 0.08), size = 3) +
  facet_wrap(~rank, scales = "free_y") + theme_minimal() +
  labs(title = "Cortical neuronal rank scores by donor",
       subtitle = "Eight donors shared with Schirmer; 18 labels weighted equally",
       x = NULL, y = "Mean UCell score")
ggsave("plots/cortex_broad_ucell.png", p, width = 7, height = 4, dpi = 180)

vip <- long[long$subtype == "VIP-IN", ]
p <- ggplot(vip, aes(condition, score, color = condition)) +
  geom_point(position = position_jitter(width = 0.08), size = 3) +
  facet_wrap(~rank, scales = "free_y") + theme_minimal() +
  labs(title = "Cortical VIP-IN rank scores by donor",
       subtitle = "Same eight donors represented in the Schirmer cohort",
       x = NULL, y = "Mean UCell score")
ggsave("plots/cortex_vip_ucell.png", p, width = 7, height = 4, dpi = 180)

p <- ggplot(long, aes(condition, score, color = condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 1.3) +
  facet_grid(rank ~ subtype, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Trobisch cortical UCell scores",
       subtitle = "One point per donor; cortex overlaps Schirmer",
       x = NULL, y = "Mean UCell score")
ggsave("plots/cortex_integrated_ucell.png", p, width = 19, height = 7, dpi = 180)
