# Sensitivity bound for the Tier 1 PSAP estimate (Analysis B) and retinas needed for 80% power. Normal approximation; SE scales as 1/sqrt(n) from the current 22 retinas.
source("../mayer_rgc_v2/common_v2.R")

bt <- read.csv("../mayer_rgc_v2/results/B_all_genes.csv", stringsAsFactors = FALSE)
r <- bt[bt$gene == "PSAP", ]
n_now <- 22; deltas <- c(0.14, 0.23, 0.40)
bound <- data.frame(estimate = r$log2FC, se = r$se, ci95_low = r$log2FC - 1.96 * r$se, ci95_high = r$log2FC + 1.96 * r$se,
                    ci90_low = r$log2FC - 1.645 * r$se, ci90_high = r$log2FC + 1.645 * r$se, min_detectable_abs_log2FC_80pct_power = 2.8 * r$se)
write.csv(bound, file.path(out_root, "bound_psap_B.csv"), row.names = FALSE)
need <- data.frame(true_log2FC = deltas, retinas_needed = ceiling(n_now * (2.8 * r$se / deltas)^2), current_retinas = n_now, se_at_current = r$se)
need$retinas_needed_if_under_current <- need$retinas_needed <= n_now
write.csv(need, file.path(out_root, "sample_size_psap.csv"), row.names = FALSE)
print(bound, digits = 4); print(need)
log_line("bound: PSAP B log2FC ", round(r$log2FC, 4), " SE ", round(r$se, 4), "; 95% interval ", round(bound$ci95_low, 3), " to ", round(bound$ci95_high, 3), "; 90% ", round(bound$ci90_low, 3), " to ", round(bound$ci90_high, 3), "; MDE ", round(bound$min_detectable_abs_log2FC_80pct_power, 3))
