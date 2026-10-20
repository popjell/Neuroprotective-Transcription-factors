# Analysis C: composition-standardized pooled test. Per retina, subtype log2 CPM combined with FIXED weights (mean subtype proportion over all 22 retinas), then limma and a permutation test.
# Needs 01_prepare.R and 03_analysis_B.R output (gene filter and eligible units from B).
source("common_v2.R")
suppressPackageStartupMessages({ library(limma); library(edgeR) })

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
B <- readRDS(file.path(out_root, "B_objects.rds"))
info <- B$info
ret <- pb$retina
g31 <- read_gene_list(pb$genes)

# log2 CPM per eligible unit on the genes kept in B
y <- calcNormFactors(DGEList(pb$pb22$pb[B$genes, info$unit]))
lcpm <- cpm(y, log = TRUE, prior.count = 1)

# fixed weights: mean proportion of each subtype over the 22 retinas (22-label), renormalized per retina over its eligible subtypes
w_all <- colMeans(pb$prop22)[unique(info$label)]
write.csv(data.frame(subtype = names(w_all), mean_proportion = as.numeric(w_all)), file.path(out_root, "C_weights.csv"), row.names = FALSE)
retinas <- ret$retina
M <- sapply(retinas, function(r) {
  u <- which(info$retina == r)
  w <- w_all[info$label[u]]
  as.numeric(lcpm[, u, drop = FALSE] %*% (w / sum(w)))
})
rownames(M) <- rownames(lcpm)
stopifnot(!anyNA(M))
cond <- ret$condition[match(retinas, ret$retina)]
is_ms <- cond == "MS"
n_ctrl <- sum(!is_ms); n_ms <- sum(is_ms)
log_line("C: ", nrow(M), " genes, ", ncol(M), " retinas (", n_ctrl, " Control, ", n_ms, " MS); subtypes used ", length(w_all), "; units per retina ", min(table(info$retina)), " to ", max(table(info$retina)))

# limma on the 22 standardized retina values
design <- model.matrix(~ cond)
fit <- eBayes(lmFit(M, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tab <- data.frame(gene = rownames(tt), log2FC = tt$logFC, se = sqrt(fit$s2.post) * fit$stdev.unscaled[, 2], t = tt$t, p = tt$P.Value, padj = tt$adj.P.Val, ave_expr = rowMeans(M), row.names = NULL)
write.csv(tab, file.path(out_root, "C_all_genes.csv"), row.names = FALSE)
r <- tab[tab$gene == "PSAP", ]
register_test("C_PSAP", 2, "C", "PSAP", r$log2FC, r$se, r$p, n_ctrl, n_ms, "composition-standardized pooled, limma ~condition, two-sided")
log_line("C: PSAP log2FC ", round(r$log2FC, 4), " SE ", round(r$se, 4), " p ", signif(r$p, 4), " padj ", signif(r$padj, 4))
r2 <- tab[tab$gene == "ITM2C", ]
register_test("C_ITM2C", 3, "C", "ITM2C", r2$log2FC, r2$se, r2$p, n_ctrl, n_ms, "composition-standardized pooled, limma ~condition")

# permutation of condition labels over retinas, PSAP only
pp <- perm_p(M["PSAP", ], is_ms)
write.csv(data.frame(gene = "PSAP", diff_ms_minus_ctrl = as.numeric(pp["obs"]), perm_p = as.numeric(pp["p"]), mode = pp["mode"]), file.path(out_root, "C_psap_permutation.csv"), row.names = FALSE)
register_test("C_PSAP_perm", 3, "C", "PSAP", as.numeric(pp["obs"]), NA, as.numeric(pp["p"]), n_ctrl, n_ms, paste0("permutation of condition labels across retinas (", pp["mode"], "), two-sided on difference of means"))
log_line("C: PSAP permutation p ", pp["p"], " (", pp["mode"], ")")

# 31-gene set: per-retina mean standardized log2 CPM over the set genes present; estimate = mean gene log2FC; test = t-test of retina set means
gs <- intersect(g31, rownames(M))
sm <- colMeans(M[gs, , drop = FALSE])
tt2 <- t.test(sm[is_ms], sm[!is_ms], var.equal = TRUE)
register_test("C_set31", 3, "C", "set31", mean(tab$log2FC[tab$gene %in% gs]), NA, tt2$p.value, n_ctrl, n_ms, paste0(length(gs), " genes; estimate = mean log2FC; p = t-test on per-retina set mean; set-mean difference ", round(diff(rev(tt2$estimate)) * -1, 4)))
write.csv(data.frame(retina = retinas, condition = cond, psap_standardized_log2cpm = M["PSAP", ], set31_mean = sm), file.path(out_root, "C_retina_values.csv"), row.names = FALSE)
log_line("C: 31-gene set (", length(gs), " genes) mean log2FC ", round(mean(tab$log2FC[tab$gene %in% gs]), 4), "; t-test p ", signif(tt2$p.value, 3))
