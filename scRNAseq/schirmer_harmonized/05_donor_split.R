# Donor split of the frozen Schirmer models: independent donors (not in Trobisch) and the 7 shared donors only (descriptive technical concordance).
# Needs results/labels.csv, results/schirmer_neurons.rds, results/primary_upperEN/psap.csv and results/pooled_limma/psap.csv from 00-03.
# The two model functions below are copied unchanged from 02_primary_upperEN.R (fit_upper) and 03_pooled_limma.R (pooled), plus a donor-count guard.
source("common.R")
library(edgeR)
library(limma)

out_dir <- res("donor_split")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "psap.csv"))

# upper_EN pseudobulk per donor (donors with at least min_cells cells), then the frozen DESeq2 model
fit_upper <- function(counts, cells) {
  d <- aggregate(cell ~ donor + condition, cells, length)
  names(d)[3] <- "n_cells"
  e <- d[d$n_cells >= min_cells, ]
  e$condition <- factor(e$condition, levels = c("Control", "MS"))
  n <- table(e$condition)
  out <- list(n = n, testable = min(n) >= min_donors)
  if (!out$testable) return(out)
  fit <- fit_deseq(pseudobulk(counts[, cells$cell], cells$donor)[, e$donor], e)
  out$res <- fit$res
  out
}

# Pseudobulk donor x group (>= min_cells), filterByExpr + TMM, ~coarse_group + condition with block = donor
pooled <- function(counts, cells) {
  cells <- cells[!is.na(cells$coarse_group) & cells$coarse_group %in% c("upper_EN", "deep_EN", "IN"), ]
  s <- aggregate(cell ~ donor + coarse_group + condition, cells, length)
  names(s)[4] <- "n_cells"
  s <- s[s$n_cells >= min_cells, ]
  s$id <- paste(s$donor, s$coarse_group, sep = "|")
  s$condition <- factor(s$condition, levels = c("Control", "MS"))
  s$coarse_group <- factor(s$coarse_group, levels = c("upper_EN", "deep_EN", "IN"))
  donors <- unique(s[, c("donor", "condition")])
  n <- table(donors$condition)
  out <- list(n = n, testable = min(n) >= min_donors)
  if (!out$testable) return(out)
  mat <- pseudobulk(counts[, cells$cell], paste(cells$donor, cells$coarse_group, sep = "|"))[, s$id]
  y <- DGEList(mat, group = s$condition)
  y <- y[filterByExpr(y, group = s$condition), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~coarse_group + condition, s)
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = s$donor)$consensus.correlation
  v <- voom(y, design, block = s$donor, correlation = cor1)
  cor2 <- duplicateCorrelation(v, design, block = s$donor)$consensus.correlation
  fit <- eBayes(lmFit(v, design, block = s$donor, correlation = cor2))
  tt <- topTable(fit, coef = "conditionMS", number = Inf, sort.by = "none", confint = TRUE)
  tt$gene <- rownames(tt)
  c(out, list(res = tt))
}

s <- readRDS(res("schirmer_neurons.rds"))
lab <- read.csv(res("labels.csv"))
stopifnot(identical(lab$cell, colnames(s$counts)))
up <- lab[lab$coarse_group == "upper_EN", ]

# One PSAP row per model for a set of donors
split_rows <- function(subset, donors) {
  fu <- fit_upper(s$counts, up[up$donor %in% donors, ])
  fp <- pooled(s$counts, lab[lab$donor %in% donors, ])
  row <- function(model, f, col) {
    r <- data.frame(dataset = "Schirmer", subset = subset, model = model, n_ctrl = f$n[["Control"]], n_ms = f$n[["MS"]],
                    log2FC = NA, ci_low = NA, ci_high = NA, p = NA, status = "not testable")
    if (!f$testable) return(r)
    x <- f$res[f$res$gene == gene, ]
    stopifnot(nrow(x) == 1)
    if (model == "primary_upperEN_deseq2") {
      r[, c("log2FC", "ci_low", "ci_high", "p")] <- c(x$log2FoldChange, x$log2FoldChange - 1.96 * x$lfcSE, x$log2FoldChange + 1.96 * x$lfcSE, x$pvalue)
    } else {
      r[, c("log2FC", "ci_low", "ci_high", "p")] <- c(x$logFC, x$CI.L, x$CI.R, x$P.Value)
    }
    r$status <- "fitted"
    r
  }
  list(rows = rbind(row("primary_upperEN_deseq2", fu), row("pooled_limma_voom", fp)), fu = fu, fp = fp)
}

# Check: all 19 donors reproduce the frozen 02 and 03 PSAP estimates
all_donors <- split_rows("all", unique(lab$donor))$rows
frozen_primary <- read.csv(res("primary_upperEN", "psap.csv"))[1, ]
frozen_pooled <- read.csv(res("pooled_limma", "psap.csv"))[1, ]
stopifnot(isTRUE(all.equal(all_donors$log2FC[1], frozen_primary$log2FC)), isTRUE(all.equal(all_donors$log2FC[2], frozen_pooled$logFC)))

split <- rbind(split_rows("independent", setdiff(unique(lab$donor), shared_donors))$rows,
               split_rows("shared_only", shared_donors)$rows)
print(split)
write.csv(split, file.path(out_dir, "psap.csv"), row.names = FALSE)
cat("donor split written;", sum(split$status == "fitted"), "of", nrow(split), "models fitted\n")
