# Secondary pooled model: all neurons, pseudobulk per donor x coarse group, limma-voom with donor blocking (duplicateCorrelation, two rounds).
# Needs results/labels.csv and results/*_neurons.rds from 00-01. Runs on Schirmer (all 19 donors) and Trobisch (comparator).
source("common.R")
library(edgeR)
library(limma)

out_dir <- res("pooled_limma")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "psap.csv"))

# Pseudobulk donor x group (>= min_cells), filterByExpr + TMM, ~coarse_group + condition with block = donor
pooled <- function(counts, cells) {
  cells <- cells[!is.na(cells$coarse_group) & cells$coarse_group %in% c("upper_EN", "deep_EN", "IN"), ]
  s <- aggregate(cell ~ donor + coarse_group + condition, cells, length)
  names(s)[4] <- "n_cells"
  s <- s[s$n_cells >= min_cells, ]
  s$id <- paste(s$donor, s$coarse_group, sep = "|")
  s$condition <- factor(s$condition, levels = c("Control", "MS"))
  s$coarse_group <- factor(s$coarse_group, levels = c("upper_EN", "deep_EN", "IN"))
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
  list(res = tt, cor = cor2, samples = s)
}

s <- readRDS(res("schirmer_neurons.rds"))
lab <- read.csv(res("labels.csv"))
stopifnot(identical(lab$cell, colnames(s$counts)))
t <- readRDS(res("trobisch_neurons.rds"))
t$cells$coarse_group <- coarse_group(t$cells$label)
fits <- list(Schirmer = pooled(s$counts, lab), Trobisch = pooled(t$counts, t$cells))

results <- do.call(rbind, lapply(names(fits), function(d) cbind(dataset = d, fits[[d]]$res)))
write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)

psap <- do.call(rbind, lapply(names(fits), function(d) {
  f <- fits[[d]]
  x <- f$res[f$res$gene == gene, ]
  stopifnot(nrow(x) == 1)
  donors <- unique(f$samples[, c("donor", "condition")])
  data.frame(dataset = ifelse(d == "Trobisch", "Trobisch discovery cohort, for effect-size comparison only", "Schirmer, all 19 donors (7 shared with Trobisch)"),
             control_donors = sum(donors$condition == "Control"), ms_donors = sum(donors$condition == "MS"), n_pseudobulks = nrow(f$samples),
             logFC = x$logFC, lower = x$CI.L, upper = x$CI.R, p = x$P.Value, consensus_correlation = f$cor,
             direction = ifelse(x$logFC > 0, "higher in MS", "lower in MS"))
}))
print(psap)
write.csv(psap, file.path(out_dir, "psap.csv"), row.names = FALSE)
cat("pooled limma:", paste(names(fits), sapply(fits, function(f) nrow(f$samples)), "pseudobulks", collapse = "; "), "\n")
