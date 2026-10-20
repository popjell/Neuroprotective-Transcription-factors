# Analysis A (Tier 3): composition diagnostics, resilient/susceptible fractions, and PSAP vs resilient fraction across retinas. Needs 01_prepare.R output.
source("common_v2.R")
suppressPackageStartupMessages(library(edgeR))

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
ret <- pb$retina
meta <- pb$meta_small
cond <- ret$condition[match(rownames(pb$prop27), ret$retina)]
is_ms <- cond == "MS"
n_ctrl <- sum(!is_ms); n_ms <- sum(is_ms)

# per-retina proportion of each of the 27 subtypes, Mann-Whitney MS vs Control
mw <- do.call(rbind, lapply(colnames(pb$prop27), function(s) {
  v <- pb$prop27[, s]
  p <- suppressWarnings(wilcox.test(v[is_ms], v[!is_ms])$p.value)
  data.frame(subtype = s, mean_prop_control = mean(v[!is_ms]), mean_prop_ms = mean(v[is_ms]), median_diff_ms_minus_ctrl = median(v[is_ms]) - median(v[!is_ms]), p = p)
}))
mw$padj_BH_over_27 <- p.adjust(mw$p, "BH")
write.csv(mw, file.path(out_root, "A_subtype_proportions.csv"), row.names = FALSE)
for (i in seq_len(nrow(mw))) register_test(paste0("A_mwu_", mw$subtype[i]), 3, "A", paste("proportion", mw$subtype[i]), mw$median_diff_ms_minus_ctrl[i], NA, mw$p[i], n_ctrl, n_ms, "Mann-Whitney on per-retina proportion; estimate = difference of medians (MS - Control)")
log_line("A: 27 subtype Mann-Whitney tests; BH<0.05 in ", sum(mw$padj_BH_over_27 < 0.05), " (", paste(mw$subtype[mw$padj_BH_over_27 < 0.05], collapse = ","), ")")

# resilient and susceptible fractions per retina and shares by condition
frac <- data.frame(retina = rownames(pb$prop27), condition = cond,
                   resilient_fraction = rowSums(pb$prop27[, resilient_groups]), susceptible_fraction = rowSums(pb$prop27[, susceptible_groups]))
write.csv(frac, file.path(out_root, "A_group_fractions_per_retina.csv"), row.names = FALSE)
share <- do.call(rbind, lapply(c("Control", "MS"), function(cn) {
  m <- meta[meta$condition == cn, ]
  data.frame(condition = cn, n_cells = nrow(m), resilient_share = mean(m$group %in% "resilient"), susceptible_share = mean(m$group %in% "susceptible"),
             resilient_fraction_mean_of_retinas = mean(frac$resilient_fraction[frac$condition == cn]), susceptible_fraction_mean_of_retinas = mean(frac$susceptible_fraction[frac$condition == cn]))
}))
write.csv(share, file.path(out_root, "A_group_share_by_condition.csv"), row.names = FALSE)
print(share)
for (g in c("resilient", "susceptible")) {
  v <- frac[[paste0(g, "_fraction")]]
  register_test(paste0("A_mwu_", g, "_fraction"), 3, "A", paste(g, "fraction"), median(v[is_ms]) - median(v[!is_ms]), NA, suppressWarnings(wilcox.test(v[is_ms], v[!is_ms])$p.value), n_ctrl, n_ms, "Mann-Whitney on per-retina fraction; estimate = difference of medians")
}

# variance decomposition of pooled log-CPM PSAP across the 22 retinas (descriptive)
y <- calcNormFactors(DGEList(pb$pb_retina))
lcpm <- cpm(y, log = TRUE, prior.count = 1)
d <- data.frame(psap = lcpm["PSAP", frac$retina], condition = cond, res = frac$resilient_fraction)
r2 <- function(f) c(r2 = summary(lm(f, d))$r.squared, adj_r2 = summary(lm(f, d))$adj.r.squared)
vd <- data.frame(model = c("psap ~ condition", "psap ~ resilient_fraction", "psap ~ condition + resilient_fraction"),
                 rbind(r2(psap ~ condition), r2(psap ~ res), r2(psap ~ condition + res)))
write.csv(vd, file.path(out_root, "A_psap_variance_decomposition.csv"), row.names = FALSE)
print(vd)
ct <- cor.test(d$psap, d$res); cs <- suppressWarnings(cor.test(d$psap, d$res, method = "spearman"))
write.csv(data.frame(pearson_r = ct$estimate, pearson_p = ct$p.value, spearman_rho = cs$estimate, spearman_p = cs$p.value, n_retinas = nrow(d)), file.path(out_root, "A_psap_vs_resilient_fraction.csv"), row.names = FALSE)
register_test("A_cor_psap_resilient_pearson", 3, "A", "PSAP", unname(ct$estimate), NA, ct$p.value, n_ctrl, n_ms, "Pearson r across 22 retinas, pooled log-CPM PSAP vs resilient fraction")
register_test("A_cor_psap_resilient_spearman", 3, "A", "PSAP", unname(cs$estimate), NA, cs$p.value, n_ctrl, n_ms, "Spearman rho across 22 retinas")
log_line("A: PSAP vs resilient fraction Pearson r ", round(ct$estimate, 3), " p ", signif(ct$p.value, 3))
write.csv(data.frame(retina = d$retina <- frac$retina, d), file.path(out_root, "A_psap_pooled_lcpm_per_retina.csv"), row.names = FALSE)
