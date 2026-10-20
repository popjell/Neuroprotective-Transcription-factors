# Secondary pooled model: all neurons, pseudobulk per donor x coarse group, limma-voom ~coarse_group + condition with donor blocking (two rounds).
# Needs results/pseudobulk_sample_group.rds and results/samples.csv from 00-01, and FROZEN_PLAN.md.
source("common.R")
library(edgeR)
library(limma)

out_dir <- res("pooled_limma")
dir.create(out_dir, showWarnings = FALSE)
out <- file.path(out_dir, "psap.csv")
no_overwrite(out)

# Pseudobulk donor x group (at least min_cells), filterByExpr + TMM
pb <- load_pb()
key <- paste(pb$cols$donor, pb$cols$group, sep = "|")
s <- aggregate(n_cells ~ donor + group + condition, pb$cols, sum)
s$id <- paste(s$donor, s$group, sep = "|")
s <- s[s$n_cells >= min_cells, ]
s$condition <- factor(s$condition, levels = c("Control", "MS"))
s$coarse_group <- factor(s$group, levels = c("upper_EN", "deep_EN", "IN"))
mat <- pseudobulk(pb$counts, key)[, s$id]
y <- DGEList(mat, group = s$condition)
y <- calcNormFactors(y[filterByExpr(y, group = s$condition), , keep.lib.sizes = FALSE], method = "TMM")

# voom with donor blocking
design <- model.matrix(~coarse_group + condition, s)
vb <- voom_block(y, design, s$donor)
tt <- topTable(eBayes(vb$fit), coef = "conditionMS", number = Inf, sort.by = "none", confint = TRUE)
tt$gene <- rownames(tt)
x <- tt[tt$gene == gene, ]
stopifnot(nrow(x) == 1)
donors <- unique(s[, c("donor", "condition")])
psap <- data.frame(dataset = "Macnair GM, all neurons", control_donors = sum(donors$condition == "Control"), ms_donors = sum(donors$condition == "MS"),
                   n_pseudobulks = nrow(s), logFC = x$logFC, lower = x$CI.L, upper = x$CI.R, p = x$P.Value, consensus_correlation = vb$cor,
                   direction = ifelse(x$logFC > 0, "higher in MS", "lower in MS"),
                   note = "control and MS donors come from different brain banks; comparison is confounded with bank")
print(psap)
write.csv(psap, out, row.names = FALSE)
write.csv(tt, file.path(out_dir, "results.csv"), row.names = FALSE)
cat("pooled limma:", nrow(s), "pseudobulks,", nrow(y), "genes\n")
