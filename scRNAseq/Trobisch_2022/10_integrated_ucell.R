# Score Rank 1 and Rank 2 in every integrated neuronal group with raw counts.
library(Matrix)
library(UCell)
library(ggplot2)

x <- readRDS("data/new_regions.rds")
atlas <- read.delim("data/meta.tsv")
neurons <- sort(unique(atlas$cell_subtypes_intg[
  grepl("-EN$|-IN$|^PC$|^motoEN$", atlas$cell_subtypes_intg)]))
keep <- x$meta$cell_subtypes_intg %in% neurons
meta <- x$meta[keep, ]
counts <- x$counts[, rownames(meta), drop = FALSE]
stopifnot(all(colSums(counts) == meta$nCount_RNA))
rm(x, atlas)
gc()

rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
sets <- lapply(sets, intersect, y = rownames(counts))

# UCell scores cells; tests below average scores by donor first.
cell_scores <- ScoreSignatures_UCell(
  matrix = counts, features = sets, maxRank = 1500,
  chunk.size = 100, ncores = 1)
cell_scores <- as.data.frame(cell_scores)
names(cell_scores) <- names(sets)
cell_scores$cellId <- rownames(cell_scores)
cell_scores <- merge(cell_scores, meta[, c("cellId", "sample_id", "condition",
                                           "cns_region", "cell_subtypes_intg")],
                     by = "cellId", sort = FALSE)
names(cell_scores)[names(cell_scores) == "cell_subtypes_intg"] <- "subtype"

donors <- aggregate(cbind(`Rank 1`, `Rank 2`) ~ sample_id + condition +
  cns_region + subtype, cell_scores, mean)
cells <- aggregate(cellId ~ sample_id + subtype, cell_scores, length)
donors$cells <- cells$cellId[match(paste(donors$sample_id, donors$subtype),
                                  paste(cells$sample_id, cells$subtype))]
long <- reshape(donors, varying = c("Rank 1", "Rank 2"), v.names = "score",
                timevar = "rank", times = c("Rank 1", "Rank 2"), direction = "long")
rownames(long) <- NULL
write.csv(long, "results/integrated_neuron_ucell_donors.csv", row.names = FALSE)

tests <- list()
for (region in unique(long$cns_region)) {
  for (subtype in unique(long$subtype[long$cns_region == region])) {
    for (rank_name in names(sets)) {
      d <- long[long$cns_region == region & long$subtype == subtype &
                  long$rank == rank_name, ]
      if (sum(d$condition == "Control") < 2 || sum(d$condition == "MS") < 2) next
      model <- lm(score ~ factor(condition, levels = c("Control", "MS")), d)
      tests[[paste(region, subtype, rank_name)]] <- data.frame(
        region = region, subtype = subtype, rank = rank_name,
        control_donors = sum(d$condition == "Control"),
        ms_donors = sum(d$condition == "MS"),
        ms_minus_control = unname(coef(model)[2]),
        p_value = summary(model)$coefficients[2, 4])
    }
  }
}
tests <- do.call(rbind, tests)
tests$q_global <- p.adjust(tests$p_value, method = "BH")
write.csv(tests, "results/integrated_neuron_ucell_tests.csv", row.names = FALSE)

# Use only labels seen in every donor of their region for the broad score.
# All other labels remain in the subtype-level results above.
coverage <- aggregate(sample_id ~ cns_region + subtype, long,
                      function(z) length(unique(z)))
region_donors <- aggregate(sample_id ~ cns_region, long,
                           function(z) length(unique(z)))
coverage$complete <- coverage$sample_id == region_donors$sample_id[
  match(coverage$cns_region, region_donors$cns_region)]
complete <- paste(coverage$cns_region[coverage$complete],
                  coverage$subtype[coverage$complete])
balanced <- long[paste(long$cns_region, long$subtype) %in% complete, ]
broad <- aggregate(score ~ sample_id + condition + cns_region + rank,
                   balanced, mean)
broad_tests <- list()
for (rank_name in names(sets)) {
  d <- broad[broad$rank == rank_name, ]
  model <- lm(score ~ factor(condition, levels = c("Control", "MS")) +
                factor(cns_region), d)
  broad_tests[[rank_name]] <- data.frame(
    rank = rank_name, donors = nrow(d),
    ms_minus_control = unname(coef(model)[2]),
    p_value = summary(model)$coefficients[2, 4])
}
broad_tests <- do.call(rbind, broad_tests)
broad_tests$q_two_lists <- p.adjust(broad_tests$p_value, method = "BH")
write.csv(broad_tests, "results/integrated_neuron_broad_ucell.csv", row.names = FALSE)

p <- ggplot(long, aes(condition, score, color = condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 1.5) +
  facet_grid(rank ~ subtype, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Trobisch UCell scores across integrated neuronal groups",
       subtitle = "One point per donor; new-region nuclei only",
       x = NULL, y = "Mean UCell score")
ggsave("plots/integrated_neuron_ucell.png", p, width = 16, height = 6, dpi = 180)
