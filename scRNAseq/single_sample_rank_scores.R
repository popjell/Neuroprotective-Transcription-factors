# Score Rank 1 and Rank 2 for every neuronal subtype using donor pseudobulk.
# Run from this folder: Rscript single_sample_rank_scores.R Schirmer
# Or: Rscript single_sample_rank_scores.R Trobisch
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, args %in% c("Schirmer", "Trobisch"))
study <- args[1]
.libPaths(c(".r_lib", .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(GSVA)
  library(singscore)
  library(AUCell)
  library(ggplot2)
})

if (study == "Schirmer") {
  x <- readRDS("Schirmer_2019/results/pseudobulk_counts.rds")
  counts <- x$counts
  info <- x$metadata
  info$donor <- info$patient
  info$condition <- info$diagnosis
  info$region <- "Cortex"
  rank_genes <- read.csv("Schirmer_2019/data/rank_genes.csv")
  gene_id <- "human_ensembl"
  min_cells <- 20
  out_dir <- "Schirmer_2019"
} else {
  cortex <- readRDS("Trobisch_2022/results/cortex_integrated_pseudobulk.rds")
  other <- readRDS("Trobisch_2022/results/integrated_neuron_pseudobulk.rds")
  stopifnot(identical(rownames(cortex$counts), rownames(other$counts)))
  cortex$metadata$cns_region <- "Cortex"
  counts <- cbind(cortex$counts, other$counts)
  info <- rbind(cortex$metadata, other$metadata)
  info$donor <- info$sample_id
  info$region <- info$cns_region
  rank_genes <- read.csv("Schirmer_2019/data/rank_genes.csv")
  gene_id <- "human_symbol"
  min_cells <- 10
  out_dir <- "Trobisch_2022"
  rm(cortex, other)
}
stopifnot(identical(colnames(counts), rownames(info)))
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", gene_id)])
sets <- split(rank_genes[[gene_id]], rank_genes$rank)
methods <- c("ssGSEA", "GSVA", "singscore", "AUCell", "PLAGE", "zscore", "minmax")
labels <- sort(unique(info$subtype))
scores <- list()
tests <- list()
set.seed(1)

for (subtype in labels) {
  # A donor-subtype sample needs enough nuclei to give a stable profile.
  use <- info$subtype == subtype & info$cells >= min_cells
  d <- info[use, ]
  group <- factor(d$condition, levels = c("Control", "MS"))
  design <- if (!nrow(d)) NULL else if (length(unique(d$region)) > 1)
    model.matrix(~factor(d$region) + group) else model.matrix(~group)
  status <- if (sum(group == "Control") < 2 || sum(group == "MS") < 2)
    "Too few donors" else if (qr(design)$rank < ncol(design) ||
                              nrow(d) - ncol(design) < 2)
    "Disease effect not estimable" else "Tested"

  for (rank_name in names(sets)) for (method in methods) {
    key <- paste(subtype, rank_name, method)
    tests[[key]] <- data.frame(
      study = study, subtype = subtype, rank = rank_name, method = method,
      control_donors = sum(group == "Control"), ms_donors = sum(group == "MS"),
      regions = paste(sort(unique(d$region)), collapse = ", "),
      genes_tested = NA_integer_, ms_minus_control = NA_real_,
      p_value = NA_real_, status = status)
  }
  if (status != "Tested") next

  y <- counts[, use, drop = FALSE]
  # Keep genes expressed in any region with at least four usable donors.
  by_region <- split(seq_len(nrow(d)), d$region)
  keep <- Reduce(`|`, lapply(by_region, function(i) {
    if (length(i) < 4) return(rep(FALSE, nrow(y)))
    filterByExpr(y[, i, drop = FALSE], group = group[i])
  }))
  if (sum(keep) < 2) {
    for (rank_name in names(sets)) for (method in methods)
      tests[[paste(subtype, rank_name, method)]]$status <- "Too few expressed genes"
    next
  }
  y <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
  expr <- as.matrix(cpm(y, log = TRUE, prior.count = 1))
  valid <- lapply(sets, intersect, y = rownames(expr))
  for (rank_name in names(sets)) for (method in methods) {
    key <- paste(subtype, rank_name, method)
    tests[[key]]$genes_tested <- length(valid[[rank_name]])
    if (length(valid[[rank_name]]) < 3)
      tests[[key]]$status <- "Fewer than three expressed list genes"
  }
  valid <- valid[lengths(valid) >= 3]
  if (!length(valid)) next

  # All seven methods score the same normalized donor expression matrix.
  score_mats <- list(
    ssGSEA = gsva(ssgseaParam(expr, valid, minSize = 3), verbose = FALSE),
    GSVA = gsva(gsvaParam(expr, valid, minSize = 3, kcdf = "Gaussian"),
                verbose = FALSE),
    PLAGE = gsva(plageParam(expr, valid, minSize = 3), verbose = FALSE),
    zscore = gsva(zscoreParam(expr, valid, minSize = 3), verbose = FALSE)
  )
  ranked <- rankGenes(expr)
  score_mats$singscore <- do.call(rbind, lapply(valid, function(g)
    simpleScore(ranked, upSet = g)$TotalScore))
  rownames(score_mats$singscore) <- names(valid)
  rankings <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
  auc <- AUCell_calcAUC(valid, rankings,
                        aucMaxRank = ceiling(nrow(expr) * 0.05), verbose = FALSE)
  score_mats$AUCell <- as.matrix(SummarizedExperiment::assay(auc))
  scaled <- t(apply(expr, 1, function(z) {
    span <- diff(range(z))
    if (span == 0) rep(0, length(z)) else (z - min(z)) / span
  }))
  rownames(scaled) <- rownames(expr)
  colnames(scaled) <- colnames(expr)
  score_mats$minmax <- do.call(rbind, lapply(valid, function(g)
    colMeans(scaled[g, , drop = FALSE])))
  rownames(score_mats$minmax) <- names(valid)

  for (method in methods) for (rank_name in names(valid)) {
    score <- as.numeric(score_mats[[method]][rank_name, ])
    if (method == "PLAGE") {
      mean_expr <- colMeans(expr[valid[[rank_name]], , drop = FALSE])
      if (is.finite(cor(score, mean_expr)) && cor(score, mean_expr) < 0)
        score <- -score
    }
    key <- paste(subtype, rank_name, method)
    if (!all(is.finite(score)) || sd(score) == 0) {
      tests[[key]]$status <- "Score unavailable or constant"
      next
    }
    scores[[key]] <- data.frame(
      study = study, subtype = subtype, rank = rank_name, method = method,
      donor = d$donor, region = d$region, condition = d$condition,
      cells = d$cells, genes_tested = length(valid[[rank_name]]), score = score)
    model <- if (length(unique(d$region)) > 1)
      lm(score ~ factor(d$region) + group) else lm(score ~ group)
    tests[[key]]$ms_minus_control <- unname(coef(model)["groupMS"])
    tests[[key]]$p_value <- summary(model)$coefficients["groupMS", 4]
  }
  message("Finished ", study, " ", subtype)
  rm(y, expr, scaled, ranked, rankings, auc, score_mats)
  gc(verbose = FALSE)
}

tests <- do.call(rbind, tests)
tests$q_method <- ave(tests$p_value, tests$method,
                      FUN = function(p) p.adjust(p, method = "BH"))
tests$q_all_methods <- p.adjust(tests$p_value, method = "BH")
rownames(tests) <- NULL
write.csv(tests, file.path(out_dir, "results/single_sample_rank_tests.csv"),
          row.names = FALSE)
if (length(scores)) {
  scores <- do.call(rbind, scores)
  rownames(scores) <- NULL
  write.csv(scores, file.path(out_dir, "results/single_sample_rank_scores.csv"),
            row.names = FALSE)
}

# Recheck Trobisch scores within each region, without mixing donor populations.
if (study == "Trobisch" && is.data.frame(scores)) {
  regional <- list()
  for (region_name in sort(unique(info$region)))
    for (subtype in labels) for (rank_name in names(sets)) for (method in methods) {
      d <- scores[scores$region == region_name & scores$subtype == subtype &
                    scores$rank == rank_name & scores$method == method, ]
      group <- factor(d$condition, levels = c("Control", "MS"))
      status <- if (sum(group == "Control") < 2 || sum(group == "MS") < 2)
        "Too few donors" else if (sd(d$score) == 0)
        "Constant score" else "Tested"
      effect <- p <- NA_real_
      if (status == "Tested") {
        fit <- lm(score ~ group, data = d)
        effect <- unname(coef(fit)["groupMS"])
        p <- summary(fit)$coefficients["groupMS", 4]
      }
      regional[[paste(region_name, subtype, rank_name, method)]] <- data.frame(
        region = region_name, subtype = subtype, rank = rank_name,
        method = method, control_donors = sum(group == "Control"),
        ms_donors = sum(group == "MS"), genes_tested = if (nrow(d))
          d$genes_tested[1] else NA_integer_, ms_minus_control = effect,
        p_value = p, status = status)
    }
  regional <- do.call(rbind, regional)
  regional$q_region_method <- ave(regional$p_value,
                                  interaction(regional$region, regional$method),
                                  FUN = function(p) p.adjust(p, method = "BH"))
  regional$q_all_regions_methods <- p.adjust(regional$p_value, method = "BH")
  rownames(regional) <- NULL
  write.csv(regional,
            file.path(out_dir, "results/single_sample_rank_regional_tests.csv"),
            row.names = FALSE)
}

# Plot the score-test status for every subtype and method.
tests$result <- ifelse(is.na(tests$q_method), "Untestable",
                       ifelse(tests$q_method < 0.05 & tests$ms_minus_control > 0,
                              "MS-up q < 0.05",
                              ifelse(tests$q_method < 0.05,
                                     "MS-down q < 0.05", "q >= 0.05")))
tests$subtype <- factor(tests$subtype, levels = rev(labels))
p <- ggplot(tests, aes(rank, subtype, fill = result)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = c("MS-up q < 0.05" = "#28689C",
                               "MS-down q < 0.05" = "#B34A3C",
                               "q >= 0.05" = "#C7D3DD",
                               "Untestable" = "#F0F0F0")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste(study, "donor-level rank-list scores"),
       subtitle = "FDR within each scoring method across testable subtype-list comparisons",
       fill = NULL)
ggsave(file.path(out_dir, "plots/single_sample_rank_tests.png"), p,
       width = 15, height = 7, dpi = 180)
