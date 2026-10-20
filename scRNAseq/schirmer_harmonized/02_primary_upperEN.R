# Primary test: PSAP in upper_EN, donor-level pseudobulk DESeq2 (MS vs Control), Schirmer (all 19 donors) and Trobisch (comparator).
# Needs results/labels.csv and results/*_neurons.rds from 00-01, and FROZEN_PLAN.md.
source("common.R")

out_dir <- res("primary_upperEN")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "psap.csv"))

# upper_EN pseudobulk per donor (donors with at least min_cells cells), then the frozen DESeq2 model
fit_upper <- function(counts, cells) {
  d <- aggregate(cell ~ donor + condition, cells, length)
  names(d)[3] <- "n_cells"
  d$eligible <- d$n_cells >= min_cells
  e <- d[d$eligible, ]
  e$condition <- factor(e$condition, levels = c("Control", "MS"))
  n <- table(e$condition)
  testable <- min(n) >= min_donors
  out <- list(donors = d, n = n, testable = testable)
  if (!testable) return(out)
  mat <- pseudobulk(counts[, cells$cell], cells$donor)[, e$donor]
  fit <- fit_deseq(mat, e)
  e$size_factor <- fit$sizes[e$donor]
  out$res <- fit$res
  out$donors <- merge(d, e[, c("donor", "size_factor")], all.x = TRUE)
  out
}

# Schirmer: cells transferred to upper_EN
s <- readRDS(res("schirmer_neurons.rds"))
lab <- read.csv(res("labels.csv"))
stopifnot(identical(lab$cell, colnames(s$counts)))
up_s <- lab[lab$coarse_group == "upper_EN", ]
fs <- fit_upper(s$counts, up_s)

# Trobisch: author labels in upper_EN (L2-EN, L2-L3-EN, L3-EN)
t <- readRDS(res("trobisch_neurons.rds"))
up_t <- t$cells[coarse_group(t$cells$label) %in% "upper_EN", ]
ft <- fit_upper(t$counts, up_t)

psap <- do.call(rbind, lapply(list(list("Schirmer, all 19 donors (7 shared with Trobisch)", fs), list("Trobisch discovery cohort, for effect-size comparison only", ft)), function(x) {
  f <- x[[2]]
  row <- data.frame(dataset = x[[1]], control_donors = f$n[["Control"]], ms_donors = f$n[["MS"]], testable = f$testable)
  if (!f$testable) return(cbind(row, log2FC = NA, lower = NA, upper = NA, p = NA, direction = "not testable"))
  r <- psap_row(f$res)
  r$direction <- ifelse(r$log2FC > 0, "higher in MS", "lower in MS")
  cbind(row, r)
}))
print(psap)
write.csv(psap, file.path(out_dir, "psap.csv"), row.names = FALSE)

if (fs$testable) {
  write.csv(fs$res, file.path(out_dir, "deseq2.csv"), row.names = FALSE)
  write.csv(transform(fs$donors, dataset = "Schirmer", shared_with_trobisch = donor %in% shared_donors), file.path(out_dir, "donors.csv"), row.names = FALSE)
}
if (ft$testable) {
  write.csv(ft$res, file.path(out_dir, "deseq2_trobisch.csv"), row.names = FALSE)
  write.csv(transform(ft$donors, dataset = "Trobisch"), file.path(out_dir, "donors_trobisch.csv"), row.names = FALSE)
}
cat("upper_EN Schirmer:", fs$n[["Control"]], "control vs", fs$n[["MS"]], "MS donors; Trobisch:", ft$n[["Control"]], "vs", ft$n[["MS"]], "\n")
