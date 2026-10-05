# Compare Rank 1 and Rank 2 expression using patients as replicates.
library(Matrix)
library(edgeR)
library(ggplot2)

x <- readRDS("results/pseudobulk_counts.rds")
counts <- x$counts
meta <- x$metadata
stopifnot(all(counts@x == round(counts@x)))
rank_genes <- read.csv("data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_ensembl", "human_symbol")])
sets <- split(rank_genes$human_ensembl, rank_genes$rank)

# Add a correction across every tested subtype-gene comparison.
de <- read.csv("results/deseq2_results.csv")
de$q_global <- p.adjust(de$pvalue, method = "BH")
rank_de <- merge(rank_genes, de, by.x = "human_ensembl", by.y = "gene_id")
write.csv(rank_de, "results/rank_gene_global.csv", row.names = FALSE)
rm(de, rank_de)

scores <- list()
tests <- list()
camera_results <- list()
for (subtype in sort(unique(meta$subtype))) {
  use <- meta$subtype == subtype & meta$cells >= 20
  d <- meta[use, ]
  if (sum(d$diagnosis == "Control") < 2 || sum(d$diagnosis == "MS") < 2) next
  y <- counts[, use, drop = FALSE]
  group <- factor(d$diagnosis, levels = c("Control", "MS"))
  keep <- filterByExpr(y, group = group)
  z <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
  logcpm <- cpm(z, log = TRUE, prior.count = 1)
  design <- model.matrix(~group)
  v <- limma::voom(z, design, plot = FALSE)
  index <- lapply(sets, function(s) which(rownames(v$E) %in% s))
  index <- index[lengths(index) >= 3]
  if (length(index)) {
    cam <- limma::camera(v, index, design, contrast = 2)
    camera_results[[subtype]] <- data.frame(
      subtype = subtype, rank = rownames(cam), genes_tested = cam$NGenes,
      direction = cam$Direction, p_value = cam$PValue)
  }
  for (rank_name in names(sets)) {
    genes <- intersect(sets[[rank_name]], rownames(logcpm))
    if (length(genes) < 3) next
    score <- colMeans(logcpm[genes, , drop = FALSE])
    scores[[paste(subtype, rank_name)]] <- data.frame(
      patient = d$patient, subtype = subtype, rank = rank_name,
      condition = d$diagnosis, cells = d$cells, score = score,
      genes_tested = length(genes))
    model <- lm(score ~ group)
    tests[[paste(subtype, rank_name)]] <- data.frame(
      subtype = subtype, rank = rank_name, genes_tested = length(genes),
      control_patients = sum(group == "Control"), ms_patients = sum(group == "MS"),
      ms_minus_control = unname(coef(model)[2]),
      p_value = summary(model)$coefficients[2, 4])
  }
}

scores <- do.call(rbind, scores)
tests <- do.call(rbind, tests)
tests$q_global <- p.adjust(tests$p_value, method = "BH")
camera_results <- do.call(rbind, camera_results)
camera_results$q_global <- p.adjust(camera_results$p_value, method = "BH")
write.csv(scores, "results/donor_rank_scores.csv", row.names = FALSE)
write.csv(tests, "results/donor_rank_tests.csv", row.names = FALSE)
write.csv(camera_results, "results/donor_rank_camera.csv", row.names = FALSE)

# Check whether any nominal subtype result changes direction when one patient is removed.
sensitivity <- list()
for (i in seq_len(nrow(tests))) {
  if (tests$p_value[i] >= 0.05) next
  d <- scores[scores$subtype == tests$subtype[i] & scores$rank == tests$rank[i], ]
  effects <- sapply(seq_len(nrow(d)), function(j) {
    z <- d[-j, ]
    unname(coef(lm(score ~ factor(condition, levels = c("Control", "MS")), z))[2])
  })
  sensitivity[[i]] <- data.frame(
    subtype = tests$subtype[i], rank = tests$rank[i],
    patients = nrow(d), same_direction = sum(sign(effects) ==
      sign(tests$ms_minus_control[i])),
    smallest_effect = min(effects), largest_effect = max(effects))
}
write.csv(do.call(rbind, sensitivity), "results/donor_rank_leave_one_out.csv",
          row.names = FALSE)

# Center each subtype on its control patients, then weight subtypes equally.
broad_results <- list()
for (rank_name in names(sets)) {
  a <- scores[scores$rank == rank_name, ]
  a$centered <- a$score - ave(seq_len(nrow(a)), a$subtype,
    FUN = function(i) mean(a$score[i][a$condition[i] == "Control"]))
  donor <- aggregate(centered ~ patient + condition, a, mean)
  model <- lm(centered ~ factor(condition, levels = c("Control", "MS")), donor)
  broad <- data.frame(rank = rank_name, patients = nrow(donor),
    ms_minus_control = unname(coef(model)[2]),
    p_value = summary(model)$coefficients[2, 4])
  broad_results[[rank_name]] <- broad
}
broad_results <- do.call(rbind, broad_results)
broad_results$q_two_lists <- p.adjust(broad_results$p_value, method = "BH")
write.csv(broad_results, "results/broad_neuron_rank_tests.csv", row.names = FALSE)

p <- ggplot(scores, aes(condition, score, color = condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 1.5) +
  facet_grid(rank ~ subtype, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Schirmer rank-gene expression by patient",
       subtitle = "One point per patient; each subtype is tested separately",
       x = NULL, y = "Mean normalized expression")
ggsave("plots/donor_rank_scores.png", p, width = 14, height = 6, dpi = 180)
