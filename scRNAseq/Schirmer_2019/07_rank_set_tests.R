# Test Rank 1 and Rank 2 in the existing Schirmer pseudobulk data.
library(Matrix)
library(edgeR)
library(limma)
library(fgsea)

x <- readRDS("results/pseudobulk_counts.rds")
counts <- x$counts
meta <- x$metadata
de <- read.csv("results/deseq2_results.csv")
rank_genes <- read.csv("data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_ensembl")])
sets <- list(
  `Rank 1` = unique(rank_genes$human_ensembl[rank_genes$rank == "Rank 1"]),
  `Rank 2` = unique(rank_genes$human_ensembl[rank_genes$rank == "Rank 2"])
)

out <- list()
for (subtype in unique(meta$subtype)) {
  use <- meta$subtype == subtype & meta$cells >= 20
  info <- meta[use, ]
  info$condition <- factor(info$diagnosis, levels = c("Control", "MS"))
  y <- counts[, use, drop = FALSE]
  y <- y[rowSums(y >= 10) >= 2, , drop = FALSE]
  design <- model.matrix(~condition, info)
  gene_sets <- lapply(sets, function(s) intersect(s, rownames(y)))

  # Rank all expressed genes by the signed DESeq2 statistic; positive means MS-up.
  d <- de[de$subtype == subtype & is.finite(de$stat), ]
  ranks <- d$stat
  names(ranks) <- d$gene_id
  fg <- fgsea(gene_sets, ranks, minSize = 1, maxSize = 500,
              BPPARAM = BiocParallel::SerialParam())
  fg$subtype <- subtype
  fg$method <- "fgsea"
  fg$direction <- ifelse(fg$NES > 0, "MS-up", "MS-down")
  out[[paste(subtype, "fgsea")]] <- data.frame(
    method = fg$method, rank = fg$pathway, subtype = subtype,
    genes_tested = fg$size, direction = fg$direction,
    effect = fg$NES, p_value = fg$pval)

  # Test each set against other genes while accounting for gene correlation.
  v <- voom(DGEList(y), design, plot = FALSE)
  index <- lapply(gene_sets, function(s) which(rownames(v$E) %in% s))
  cam <- camera(v, index, design, contrast = 2)
  out[[paste(subtype, "camera")]] <- data.frame(
    method = "camera", rank = rownames(cam), subtype = subtype,
    genes_tested = cam$NGenes, direction = paste0("MS-", cam$Direction),
    effect = NA_real_, p_value = cam$PValue)

  # Test whether genes in each set show an overall effect.
  set.seed(1)
  roast_result <- mroast(v, index, design, contrast = 2, nrot = 9999)
  out[[paste(subtype, "roast")]] <- data.frame(
    method = "roast", rank = rownames(roast_result), subtype = subtype,
    genes_tested = roast_result$NGenes, direction = roast_result$Direction,
    effect = NA_real_, p_value = roast_result$PValue)

  # Compare each patient's average normalized expression of the set.
  logcpm <- cpm(DGEList(y), log = TRUE, prior.count = 1)
  for (rank_name in names(gene_sets)) {
    s <- gene_sets[[rank_name]]
    score <- colMeans(logcpm[s, , drop = FALSE])
    test <- wilcox.test(score ~ info$diagnosis, exact = FALSE)
    out[[paste(subtype, "score", rank_name)]] <- data.frame(
      method = "pseudobulk_mean_logCPM", rank = rank_name,
      subtype = subtype, genes_tested = length(s),
      direction = ifelse(mean(score[info$diagnosis == "MS"]) >
                           mean(score[info$diagnosis == "Control"]), "MS-up", "MS-down"),
      effect = mean(score[info$diagnosis == "MS"]) -
        mean(score[info$diagnosis == "Control"]), p_value = test$p.value)
  }
}

results <- do.call(rbind, out)
results$q_within_method <- ave(results$p_value, results$method,
                               FUN = function(p) p.adjust(p, method = "BH"))
results$q_all_analyses <- p.adjust(results$p_value, method = "BH")
write.csv(results, "results/rank_set_method_comparison.csv", row.names = FALSE)
