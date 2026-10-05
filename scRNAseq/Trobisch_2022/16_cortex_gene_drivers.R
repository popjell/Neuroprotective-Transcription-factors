# Check whether a few genes explain the cortical UCell difference.
library(Matrix)
library(UCell)

x <- readRDS("data/cortex_neurons.rds")
counts <- x$counts
meta <- x$meta
rm(x)
gc()
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
r1 <- intersect(rank_genes$human_symbol[rank_genes$rank == "Rank 1"], rownames(counts))
r2 <- intersect(rank_genes$human_symbol[rank_genes$rank == "Rank 2"], rownames(counts))
sets <- list(
  R1_full = r1,
  R1_no_PSAP = setdiff(r1, "PSAP"),
  R1_no_PSAP_PREX1 = setdiff(r1, c("PSAP", "PREX1")),
  R2_full = r2,
  R2_no_PSAP = setdiff(r2, "PSAP"),
  R2_no_ITM2C = setdiff(r2, "ITM2C"),
  R2_no_PSAP_ITM2C = setdiff(r2, c("PSAP", "ITM2C")),
  R2_no_PSAP_ITM2C_CNDP2 = setdiff(r2, c("PSAP", "ITM2C", "CNDP2")))

scores <- ScoreSignatures_UCell(matrix = counts, features = sets,
                                maxRank = 1500, chunk.size = 100, ncores = 1)
scores <- as.data.frame(scores)
names(scores) <- names(sets)
scores$cellId <- rownames(scores)
scores <- merge(scores, meta[, c("cellId", "sample_id", "condition",
                                "cell_subtypes_intg")], by = "cellId", sort = FALSE)
names(scores)[names(scores) == "cell_subtypes_intg"] <- "subtype"

coverage <- aggregate(cellId ~ sample_id + subtype, scores, length)
complete <- aggregate(sample_id ~ subtype, coverage,
                      function(z) length(unique(z)))
complete <- complete$subtype[complete$sample_id == length(unique(scores$sample_id))]
scores <- scores[scores$subtype %in% complete, ]

tests <- list()
donors <- list()
for (name in names(sets)) {
  d <- aggregate(scores[[name]],
                 by = list(sample_id = scores$sample_id,
                           condition = scores$condition,
                           subtype = scores$subtype), mean)
  names(d)[4] <- "score"
  broad <- aggregate(score ~ sample_id + condition, d, mean)
  broad$set <- name
  donors[[name]] <- broad
  group <- factor(broad$condition, levels = c("Control", "MS"))
  model <- lm(broad$score ~ group)
  tests[[name]] <- data.frame(
    set = name, genes_in_matrix = length(sets[[name]]),
    ms_minus_control = unname(coef(model)[2]),
    p_value = summary(model)$coefficients[2, 4])
}
tests <- do.call(rbind, tests)
tests$q_all_variants <- p.adjust(tests$p_value, method = "BH")
write.csv(tests, "results/cortex_gene_driver_tests.csv", row.names = FALSE)
write.csv(do.call(rbind, donors), "results/cortex_gene_driver_donors.csv",
          row.names = FALSE)
