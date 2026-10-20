# Primary test: PSAP in upper_EN, donor-level pseudobulk DESeq2 (~condition, MS vs Control), Macnair grey matter.
# Needs results/pseudobulk_sample_group.rds and results/samples.csv from 00-01, and FROZEN_PLAN.md.
source("common.R")

out_dir <- res("primary_upperEN")
dir.create(out_dir, showWarnings = FALSE)
out <- file.path(out_dir, "psap.csv")
no_overwrite(out)

# One upper_EN pseudobulk per donor (all GM samples summed); donors with at least min_cells cells
pb <- load_pb()
up <- pb$cols[pb$cols$group == "upper_EN", ]
d <- aggregate(n_cells ~ donor + condition, up, sum)
d$eligible <- d$n_cells >= min_cells
e <- d[d$eligible, ]
e$condition <- factor(e$condition, levels = c("Control", "MS"))
rownames(e) <- e$donor
n <- table(e$condition)
stopifnot(nrow(e) == 34, min(n) >= min_donors)
mat <- pseudobulk(pb$counts[, up$id], up$donor)[, e$donor]

# DESeq2 ~condition, Control reference
fit <- fit_deseq(mat, e, ~condition, c("condition", "MS", "Control"))
e$size_factor <- fit$sizes[e$donor]
e$psap_count <- mat[gene, e$donor]
psap <- cbind(dataset = "Macnair GM, upper_EN", control_donors = n[["Control"]], ms_donors = n[["MS"]], psap_row(fit$res))
psap$direction <- ifelse(psap$log2FC > 0, "higher in MS", "lower in MS")
psap$note <- "control and MS donors come from different brain banks (Edinburgh BB vs UK MS TB); comparison is confounded with bank"
print(psap)
write.csv(psap, out, row.names = FALSE)
write.csv(fit$res, file.path(out_dir, "deseq2.csv"), row.names = FALSE)
write.csv(merge(d, e[, c("donor", "size_factor", "psap_count")], all.x = TRUE), file.path(out_dir, "donors.csv"), row.names = FALSE)
cat("upper_EN Macnair:", n[["Control"]], "control vs", n[["MS"]], "MS donors,", nrow(fit$res), "genes tested\n")
