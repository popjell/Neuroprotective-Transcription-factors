# Cell-level analysis (supportive, not independent evidence): PSAP, 31-gene UCell module, leave-PSAP-out module; MS vs Control.
# Needs mayer_ucell.R output. Env: TRIM=1 drops the top 5% microglia-score cells and writes to results/celllevel_trimmed.
source("common.R")
suppressPackageStartupMessages(library(lmerTest))

trim <- Sys.getenv("TRIM", "0") == "1"
out_dir <- file.path(res_dir, if (trim) "celllevel_trimmed" else "celllevel")
dir.create(out_dir, showWarnings = FALSE)
log_file <- file.path(out_dir, "run_celllevel.log")
no_overwrite(log_file)
sink(log_file, split = TRUE)

# load
meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
stopifnot(identical(attr(x, "cells")[[1]], meta$cell_barcode))
sc <- read.csv(file.path(res_dir, "contamination", "cell_scores.csv"), stringsAsFactors = FALSE)
stopifnot(identical(sc$cell_barcode, meta$cell_barcode))
psap <- gene_lognorm(x, genes, primary_gene)
cells <- data.frame(sample = sc$sample, condition = ifelse(sc$condition == "Ctrl", "Control", "MS"), celltypes = sc$celltypes,
                    psap_expr = psap$lognorm[primary_gene, ], psap_detect = as.numeric(psap$counts[primary_gene, ] > 0),
                    module31 = sc$module31, module_no_psap = sc$module_no_psap)
if (trim) {
  keep <- sc$microglia_ucell <= quantile(sc$microglia_ucell, trim_quantile)
  cat("contamination trim: dropping", sum(!keep), "of", length(keep), "cells\n")
  cells <- cells[keep, ]
}
cells$condition <- factor(cells$condition, levels = c("Control", "MS"))
scores <- c("psap_expr", "psap_detect", "module31", "module_no_psap")

lmm <- list(); perm <- list(); direction <- list(); wil <- list(); means_out <- list()
groups <- c("all_RGC", sort(unique(cells$celltypes)))
for (grp in groups) {
  d <- if (grp == "all_RGC") cells else cells[cells$celltypes == grp, ]
  n_cells <- table(d$sample)
  ok <- names(n_cells)[n_cells >= min_cells]
  d <- d[d$sample %in% ok, ]
  cond <- tapply(as.character(d$condition), d$sample, unique)
  n_ctrl <- sum(cond == "Control"); n_ms <- sum(cond == "MS")
  cat("\n==", grp, ":", n_ctrl, "Control vs", n_ms, "MS samples with >=", min_cells, "cells;", nrow(d), "cells\n")
  if (n_ctrl < min_donors || n_ms < min_donors) { cat("skipped\n"); next }
  d$sample <- factor(d$sample)
  for (s in scores) {
    # 1. mixed model, sample as random intercept (donor-aware; for psap_detect this is a linear probability model)
    fit <- lmer(as.formula(paste(s, "~ condition + (1 | sample)")), data = d)
    co <- summary(fit)$coefficients["conditionMS", ]
    lmm[[paste(grp, s)]] <- data.frame(group = grp, score = s, n_control_samples = n_ctrl, n_ms_samples = n_ms, estimate_ms_minus_ctrl = co["Estimate"],
                                       se = co["Std. Error"], df = co["df"], p = co["Pr(>|t|)"], singular = isSingular(fit))
    # 2. sample-level means and exact permutation test on the difference of means
    m <- tapply(d[[s]], d$sample, mean)
    is_ms <- cond[names(m)] == "MS"
    perm[[paste(grp, s)]] <- cbind(group = grp, score = s, n_control_samples = n_ctrl, n_ms_samples = n_ms, perm_test(as.numeric(m), is_ms))
    # 3. direction consistency
    mx_ctrl <- max(m[!is_ms]); mn_ms <- min(m[is_ms])
    direction[[paste(grp, s)]] <- data.frame(group = grp, score = s, n_ms = n_ms, n_ms_above_max_control = sum(m[is_ms] > mx_ctrl),
                                             min_ms_over_max_control = mn_ms / mx_ctrl, mean_control = mean(m[!is_ms]), mean_ms = mean(m[is_ms]))
    # cells as replicates (per instruction; NOT donor-aware): Wilcoxon and OLS on cells
    ols <- summary(lm(as.formula(paste(s, "~ condition")), data = d))$coefficients["conditionMS", ]
    wil[[paste(grp, s)]] <- data.frame(group = grp, score = s, n_cells_control = sum(d$condition == "Control"), n_cells_ms = sum(d$condition == "MS"),
                                       ols_estimate = ols["Estimate"], ols_se = ols["Std. Error"], ols_p = ols["Pr(>|t|)"],
                                       wilcoxon_p = wilcox.test(d[[s]][d$condition == "MS"], d[[s]][d$condition == "Control"])$p.value)
    means_out[[paste(grp, s)]] <- data.frame(group = grp, score = s, sample = names(m), condition = cond[names(m)], n_cells = as.integer(n_cells[names(m)]), sample_mean = as.numeric(m))
  }
}

# write
bind <- function(l) { r <- do.call(rbind, l); rownames(r) <- NULL; r }
write.csv(bind(lmm), file.path(out_dir, "mixed_model.csv"), row.names = FALSE)
write.csv(bind(perm), file.path(out_dir, "permutation.csv"), row.names = FALSE)
write.csv(bind(direction), file.path(out_dir, "direction_consistency.csv"), row.names = FALSE)
write.csv(bind(wil), file.path(out_dir, "cells_as_replicates.csv"), row.names = FALSE)
write.csv(bind(means_out), file.path(out_dir, "celllevel_donor_means.csv"), row.names = FALSE)  # unit = sample (retina)
cat("\nLabels: mixed_model, permutation, direction_consistency = supportive, donor-aware (sample as unit). cells_as_replicates = cells treated as independent, not donor-aware, supportive only.\n")
sink()
