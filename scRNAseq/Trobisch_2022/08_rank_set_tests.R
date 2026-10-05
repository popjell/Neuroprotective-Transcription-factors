# Test Rank 1 and Rank 2 in the existing Trobisch pseudobulk data.
library(Matrix)
library(edgeR)
library(limma)
library(fgsea)

x <- readRDS("results/pseudobulk_counts.rds")
counts <- x$counts
meta <- x$metadata
de <- read.csv("results/deseq2_results.csv")
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- list(
  `Rank 1` = unique(rank_genes$human_symbol[rank_genes$rank == "Rank 1"]),
  `Rank 2` = unique(rank_genes$human_symbol[rank_genes$rank == "Rank 2"])
)

out <- list()
for (region in unique(meta$cns_region)) {
  for (subtype in unique(meta$subtype[meta$cns_region == region])) {
    use <- meta$cns_region == region & meta$subtype == subtype
    info <- meta[use, ]
    info$condition <- factor(info$condition, levels = c("Control", "MS"))
    info$sex <- factor(info$sex)
    y <- counts[, use, drop = FALSE]
    y <- y[rowSums(y >= 10) >= 2, , drop = FALSE]
    design <- model.matrix(~sex + condition, info)
    gene_sets <- lapply(sets, function(s) intersect(s, rownames(y)))
    testable <- names(gene_sets)[lengths(gene_sets) >= 3]

    # Rank all genes by the sex-adjusted DESeq2 statistic; positive means MS-up.
    d <- de[de$region == region & de$subtype == subtype & is.finite(de$stat), ]
    ranks <- d$stat
    names(ranks) <- d$gene
    if (length(testable)) {
      fg <- fgsea(gene_sets[testable], ranks, minSize = 3, maxSize = 500,
                  BPPARAM = BiocParallel::SerialParam())
    } else fg <- data.frame()
    for (rank_name in names(gene_sets)) {
      if (rank_name %in% testable && rank_name %in% fg$pathway) {
        f <- fg[fg$pathway == rank_name, ]
        out[[paste(region, subtype, "fgsea", rank_name)]] <- data.frame(
          method = "fgsea", rank = rank_name, region = region, subtype = subtype,
          genes_tested = length(gene_sets[[rank_name]]),
          direction = ifelse(f$NES > 0, "MS-up", "MS-down"),
          effect = f$NES, p_value = f$pval)
      } else {
        out[[paste(region, subtype, "fgsea", rank_name)]] <- data.frame(
          method = "fgsea", rank = rank_name, region = region, subtype = subtype,
          genes_tested = length(gene_sets[[rank_name]]),
          direction = "Not testable (<3 genes)", effect = NA_real_, p_value = NA_real_)
      }
    }

    # Compare the set with other genes while accounting for gene correlation.
    v <- voom(DGEList(y), design, plot = FALSE)
    index <- lapply(gene_sets[testable], function(s) which(rownames(v$E) %in% s))
    if (length(testable)) cam <- camera(v, index, design, contrast = "conditionMS")
    for (rank_name in names(gene_sets)) {
      if (rank_name %in% testable) {
        out[[paste(region, subtype, "camera", rank_name)]] <- data.frame(
          method = "camera", rank = rank_name, region = region, subtype = subtype,
          genes_tested = cam[rank_name, "NGenes"],
          direction = paste0("MS-", cam[rank_name, "Direction"]),
          effect = NA_real_, p_value = cam[rank_name, "PValue"])
      } else {
        out[[paste(region, subtype, "camera", rank_name)]] <- data.frame(
          method = "camera", rank = rank_name, region = region, subtype = subtype,
          genes_tested = length(gene_sets[[rank_name]]),
          direction = "Not testable (<3 genes)", effect = NA_real_, p_value = NA_real_)
      }
    }

    # Test whether genes in each set show an overall effect.
    set.seed(1)
    if (length(testable)) {
      roast_result <- mroast(v, index, design, contrast = "conditionMS", nrot = 9999)
    }
    for (rank_name in names(gene_sets)) {
      if (rank_name %in% testable) {
        out[[paste(region, subtype, "roast", rank_name)]] <- data.frame(
          method = "roast", rank = rank_name, region = region, subtype = subtype,
          genes_tested = roast_result[rank_name, "NGenes"],
          direction = roast_result[rank_name, "Direction"], effect = NA_real_,
          p_value = roast_result[rank_name, "PValue"])
      } else {
        out[[paste(region, subtype, "roast", rank_name)]] <- data.frame(
          method = "roast", rank = rank_name, region = region, subtype = subtype,
          genes_tested = length(gene_sets[[rank_name]]),
          direction = "Not testable (<3 genes)", effect = NA_real_, p_value = NA_real_)
      }
    }

    # Compare patient-level average normalized expression of each gene set.
    logcpm <- cpm(DGEList(y), log = TRUE, prior.count = 1)
    for (rank_name in names(gene_sets)) {
      s <- gene_sets[[rank_name]]
      if (length(s) >= 3) {
        score <- colMeans(logcpm[s, , drop = FALSE])
        test <- wilcox.test(score ~ info$condition, exact = FALSE)
        direction <- ifelse(mean(score[info$condition == "MS"]) >
                              mean(score[info$condition == "Control"]), "MS-up", "MS-down")
        effect <- mean(score[info$condition == "MS"]) -
          mean(score[info$condition == "Control"])
        p_value <- test$p.value
      } else {
        direction <- "Not testable (<3 genes)"
        effect <- p_value <- NA_real_
      }
      out[[paste(region, subtype, "score", rank_name)]] <- data.frame(
        method = "pseudobulk_mean_logCPM", rank = rank_name,
        region = region, subtype = subtype, genes_tested = length(s),
        direction = direction, effect = effect, p_value = p_value)
    }
  }
}

# Add the existing patient-averaged UCell tests for the same six populations.
ucell <- read.csv("results/ucell_rank_list_tests.csv")
names(ucell)[names(ucell) == "rank_list"] <- "rank"
ucell$method <- "UCell_existing"
ucell$genes_tested <- vapply(seq_len(nrow(ucell)), function(i) {
  length(intersect(sets[[ucell$rank[i]]], rownames(counts)))
}, integer(1))
ucell$direction <- ifelse(ucell$ms_minus_control > 0, "MS-up", "MS-down")
ucell$effect <- ucell$ms_minus_control
ucell$region <- ucell$region
ucell <- ucell[, c("method", "rank", "region", "subtype", "genes_tested",
                   "direction", "effect", "p_value")]

results <- rbind(do.call(rbind, out), ucell)
results$q_within_method <- ave(results$p_value, results$method,
                               FUN = function(p) p.adjust(p, method = "BH"))
results$q_all_analyses <- p.adjust(results$p_value, method = "BH")
write.csv(results, "results/rank_set_method_comparison.csv", row.names = FALSE)
