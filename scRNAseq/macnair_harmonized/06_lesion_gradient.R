# Secondary: PSAP gradient control GM -> NAGM -> GML in upper_EN, one pseudobulk per donor x tissue, limma-voom ~tissue with donor blocking (two rounds).
# Needs results/pseudobulk_sample_group.rds and results/samples.csv from 00-01, and FROZEN_PLAN.md.
source("common.R")
library(edgeR)
library(limma)

out_dir <- res("lesion_gradient")
dir.create(out_dir, showWarnings = FALSE)
out <- file.path(out_dir, "psap.csv")
no_overwrite(out)

# upper_EN pseudobulk per donor x tissue (samples summed), at least min_cells cells
pb <- load_pb()
up <- pb$cols[pb$cols$group == "upper_EN", ]
up$key <- paste(up$donor, up$tissue, sep = "|")
up$n_samples <- 1
s <- aggregate(cbind(n_cells, n_samples) ~ donor + tissue + condition + key, up, sum)
s <- s[s$n_cells >= min_cells, ]
s$tissue <- factor(s$tissue, levels = c("control GM", "NAGM", "GML"))
mat <- pseudobulk(pb$counts[, up$id], up$key)[, s$key]
y <- DGEList(mat, group = s$tissue)
y <- calcNormFactors(y[filterByExpr(y, group = s$tissue), , keep.lib.sizes = FALSE], method = "TMM")

# voom with donor blocking; contrasts against control GM and GML vs NAGM
design <- model.matrix(~tissue, s)
vb <- voom_block(y, design, s$donor)
cm <- cbind("NAGM vs control GM" = c(0, 1, 0), "GML vs control GM" = c(0, 0, 1), "GML vs NAGM" = c(0, -1, 1))
rownames(cm) <- colnames(design)
fit <- eBayes(contrasts.fit(vb$fit, cm))
tt <- do.call(rbind, lapply(colnames(cm), function(k) {
  x <- topTable(fit, coef = k, number = Inf, sort.by = "none", confint = TRUE)
  cbind(contrast = k, gene = rownames(x), x)
}))
x <- tt[tt$gene == gene, ]
stopifnot(nrow(x) == 3)
n <- table(s$tissue)
psap <- data.frame(contrast = x$contrast, control_gm_pseudobulks = n[["control GM"]], nagm_pseudobulks = n[["NAGM"]], gml_pseudobulks = n[["GML"]],
                   logFC = x$logFC, lower = x$CI.L, upper = x$CI.R, p = x$P.Value, consensus_correlation = vb$cor,
                   direction = ifelse(x$logFC > 0, "higher in first tissue", "lower in first tissue"),
                   note = "contrasts with control GM are confounded with brain bank; GML vs NAGM is within MS donors")
print(psap)
write.csv(psap, out, row.names = FALSE)
write.csv(tt, file.path(out_dir, "results.csv"), row.names = FALSE)
write.csv(s[, c("donor", "tissue", "condition", "n_cells", "n_samples")], file.path(out_dir, "pseudobulks.csv"), row.names = FALSE)
cat("gradient:", nrow(s), "pseudobulks (", paste(names(n), n, collapse = ", "), "),", nrow(y), "genes\n")
