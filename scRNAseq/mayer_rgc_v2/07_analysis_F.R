# Analysis F (Tier 3, supportive): cell-level models with subtype and a retina random intercept. Reads the baseline UCell scores read-only.
source("common_v2.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

cnt <- load_counts_ro()
x <- cnt$x; genes <- cnt$genes
sc <- read.csv(file.path(base_res, "contamination", "cell_scores.csv"), stringsAsFactors = FALSE)  # baseline UCell scores (module31 = 15 genes >= 5% detection)
stopifnot(identical(sc$cell_barcode, cnt$cells))
mg <- read.csv(file.path(base_res, "celllevel", "module_genes.csv"), stringsAsFactors = FALSE)
log_line("F: baseline module genes (>= 5% detection): ", sum(mg$in_module), " (", paste(mg$gene[mg$in_module], collapse = ","), ")")

psap_cnt <- as.numeric(x[match("PSAP", genes), ])
cells <- data.frame(retina = sc$sample, condition = factor(ifelse(sc$condition == "Ctrl", "Control", "MS"), levels = c("Control", "MS")), subtype = factor(sc$celltypes),
                    psap_lognorm = log1p(psap_cnt / colSums(x) * 1e4), psap_detect = as.integer(psap_cnt > 0), module31 = sc$module31, module_no_psap = sc$module_no_psap)
n_ret <- tapply(as.character(cells$condition), cells$retina, unique)
nc <- sum(n_ret == "Control"); nm <- sum(n_ret == "MS")
log_line("F: cells ", nrow(cells), "; retinas ", nc, " Control, ", nm, " MS; subtype levels ", nlevels(cells$subtype))

out <- list()
for (s in c("psap_lognorm", "module31", "module_no_psap")) {
  fit <- lmer(as.formula(paste(s, "~ condition + subtype + (1 | retina)")), data = cells)
  co <- summary(fit)$coefficients["conditionMS", ]
  out[[s]] <- data.frame(score = s, model = "lmer", estimate = co["Estimate"], se = co["Std. Error"], df = co["df"], p = co["Pr(>|t|)"], n_control = nc, n_ms = nm, singular = isSingular(fit))
  register_test(paste0("F_", s), 3, "F", s, co["Estimate"], co["Std. Error"], co["Pr(>|t|)"], nc, nm, "cell-level lmer score ~ condition + subtype + (1|retina); supportive")
  log_line("F: ", s, " estimate ", signif(co["Estimate"], 4), " SE ", signif(co["Std. Error"], 4), " p ", signif(co["Pr(>|t|)"], 4))
}
fit <- glmer(psap_detect ~ condition + subtype + (1 | retina), data = cells, family = binomial)
co <- summary(fit)$coefficients["conditionMS", ]
out[["psap_detect"]] <- data.frame(score = "psap_detect", model = "glmer binomial (log odds)", estimate = co["Estimate"], se = co["Std. Error"], df = NA, p = co["Pr(>|z|)"], n_control = nc, n_ms = nm, singular = isSingular(fit))
register_test("F_psap_detect", 3, "F", "psap_detect", co["Estimate"], co["Std. Error"], co["Pr(>|z|)"], nc, nm, "cell-level glmer binomial detection ~ condition + subtype + (1|retina); estimate in log odds; supportive")
log_line("F: psap_detect log-odds ", signif(co["Estimate"], 4), " SE ", signif(co["Std. Error"], 4), " p ", signif(co["Pr(>|z|)"], 4))
write.csv(do.call(rbind, out), file.path(out_root, "F_celllevel_subtype_adjusted.csv"), row.names = FALSE)
