# limma-voom with retina as a blocking factor (two-pass duplicateCorrelation), used by Analysis B and D.
suppressPackageStartupMessages({ library(limma); library(edgeR) })

# counts: genes x units; info: one row per unit (retina, label, condition); design: model matrix for the units
fit_voom_block <- function(counts, info, design, quality_weights = TRUE) {
  y <- DGEList(counts[, info$unit])
  keep <- filterByExpr(y, group = info$condition)  # default min.count 10, min.total.count 15
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  block <- info$retina
  vfun <- if (quality_weights) function(...) voomWithQualityWeights(y, design, plot = FALSE, ...) else function(...) voom(y, design, plot = FALSE, ...)
  v <- vfun()
  cor1 <- duplicateCorrelation(v, design, block = block)$consensus.correlation
  v <- vfun(block = block, correlation = cor1)
  cor2 <- duplicateCorrelation(v, design, block = block)$consensus.correlation
  fit <- eBayes(lmFit(v, design, block = block, correlation = cor2))
  list(fit = fit, v = v, y = y, cor = cor2, n_genes = sum(keep), design = design, block = block)
}

# condition coefficient table: log2FC, SE, t, p, BH padj
coef_table <- function(m, coef) {
  tt <- topTable(m$fit, coef = coef, number = Inf, sort.by = "none")
  se <- sqrt(m$fit$s2.post) * m$fit$stdev.unscaled[, coef]
  data.frame(gene = rownames(tt), log2FC = tt$logFC, se = se, t = tt$t, p = tt$P.Value, padj = tt$adj.P.Val, ave_expr = tt$AveExpr, row.names = NULL)
}

# keep labels that have eligible units (>= min_cells) in at least min_per_cond retinas of each condition
eligible_units <- function(info, min_cells, min_per_cond) {
  info <- info[info$n_cells >= min_cells, ]
  ok <- sapply(split(info, info$label), function(d) sum(d$condition == "Control") >= min_per_cond && sum(d$condition == "MS") >= min_per_cond)
  info[info$label %in% names(ok)[ok], ]
}
