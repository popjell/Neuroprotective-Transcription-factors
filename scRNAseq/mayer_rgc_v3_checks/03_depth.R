# Depth-adjusted cell-level PSAP models (Tier 3, supportive): condition estimate with and without log total counts, and retina-level depth by condition.
source("../mayer_rgc_v2/common_v2.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

cnt <- load_counts_ro()
x <- cnt$x; genes <- cnt$genes
meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
stopifnot(identical(cnt$cells, meta$cell_barcode), all(meta$nCount_RNA == colSums(x)))
psap <- as.numeric(x[match("PSAP", genes), ])
cells <- data.frame(retina = paste0("S", meta$Sample), condition = factor(ifelse(meta$Condition == "Ctrl", "Control", "MS"), levels = c("Control", "MS")), subtype = factor(meta$celltypes),
                    psap_lognorm = log1p(psap / meta$nCount_RNA * 1e4), psap_detect = as.integer(psap > 0),
                    log_total_counts = log(meta$nCount_RNA), log_n_genes = log(meta$nFeature_RNA))
cells$log_total_counts <- cells$log_total_counts - mean(cells$log_total_counts)  # centered
cells$log_n_genes <- cells$log_n_genes - mean(cells$log_n_genes)
log_line("depth: ", nrow(cells), " nuclei, ", length(unique(cells$retina)), " retinas; PSAP detected in ", round(100 * mean(cells$psap_detect), 1), "% of nuclei")

# retina-level depth by condition, exact permutation
dep <- data.frame(retina = unique(cells$retina))
dep$condition <- cells$condition[match(dep$retina, cells$retina)]
dep$mean_total_counts <- as.numeric(tapply(meta$nCount_RNA, cells$retina, mean)[dep$retina])
dep$mean_n_genes <- as.numeric(tapply(meta$nFeature_RNA, cells$retina, mean)[dep$retina])
dep$n_nuclei <- as.integer(table(cells$retina)[dep$retina])
write.csv(dep, file.path(out_root, "depth_per_retina.csv"), row.names = FALSE)
dt <- do.call(rbind, lapply(c("mean_total_counts", "mean_n_genes"), function(v) {
  y <- dep[[v]]; is_ms <- dep$condition == "MS"
  pp <- perm_p(y, is_ms)
  register_test(paste0("DEPTH_", v), 3, "depth", v, as.numeric(pp["obs"]), NA, as.numeric(pp["p"]), 10, 12, "retina-level mean per nucleus, MS - Control; exact permutation")
  data.frame(metric = v, mean_control = mean(y[!is_ms]), mean_ms = mean(y[is_ms]), diff_ms_minus_ctrl = as.numeric(pp["obs"]), perm_p = as.numeric(pp["p"]), mode = pp["mode"])
}))
write.csv(dt, file.path(out_root, "depth_by_condition.csv"), row.names = FALSE)
print(dt, digits = 4)

# models: M0 no depth, M1 + log total counts, M2 + log total counts + log n genes
forms <- list(M0 = "~ condition + subtype + (1 | retina)", M1 = "~ condition + subtype + log_total_counts + (1 | retina)", M2 = "~ condition + subtype + log_total_counts + log_n_genes + (1 | retina)")
out <- list()
for (mn in names(forms)) {
  fit <- lmer(as.formula(paste("psap_lognorm", forms[[mn]])), data = cells)
  co <- summary(fit)$coefficients["conditionMS", ]
  out[[paste(mn, "lognorm")]] <- data.frame(model = mn, outcome = "PSAP log1p CP10K (lmer)", estimate = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|t|)"], depth_coef = if (mn != "M0") summary(fit)$coefficients["log_total_counts", "Estimate"] else NA, singular = isSingular(fit))
  register_test(paste0("DEPTH_lognorm_", mn), 3, "depth", "PSAP", co["Estimate"], co["Std. Error"], co["Pr(>|t|)"], 10, 12, paste("lmer PSAP log1p CP10K", forms[[mn]]))
  fit <- glmer(as.formula(paste("psap_detect", forms[[mn]])), data = cells, family = binomial)
  co <- summary(fit)$coefficients["conditionMS", ]
  out[[paste(mn, "detect")]] <- data.frame(model = mn, outcome = "PSAP detection (glmer binomial, log odds)", estimate = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|z|)"], depth_coef = if (mn != "M0") summary(fit)$coefficients["log_total_counts", "Estimate"] else NA, singular = isSingular(fit))
  register_test(paste0("DEPTH_detect_", mn), 3, "depth", "PSAP", co["Estimate"], co["Std. Error"], co["Pr(>|z|)"], 10, 12, paste("glmer binomial PSAP detection", forms[[mn]], "; estimate in log odds"))
  log_line("depth ", mn, ": lognorm estimate ", signif(out[[paste(mn, "lognorm")]]$estimate, 4), " p ", signif(out[[paste(mn, "lognorm")]]$p, 3), "; detection log-odds ", signif(out[[paste(mn, "detect")]]$estimate, 4), " p ", signif(out[[paste(mn, "detect")]]$p, 3))
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, file.path(out_root, "depth_adjusted_models.csv"), row.names = FALSE)
print(res, digits = 4)
