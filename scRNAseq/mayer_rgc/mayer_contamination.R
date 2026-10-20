# Contamination check: microglia, astrocyte/Mueller and oligodendrocyte marker expression per sample, MS vs Control. Needs no ambient correction (none applied; the data are plain QC-filtered counts).
source("common.R")

out_dir <- file.path(res_dir, "contamination")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "marker_condition_tests.csv"))

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
stopifnot(identical(attr(x, "cells")[[1]], meta$cell_barcode))
sets <- list(microglia = microglia_markers, astro_muller = astro_markers, oligo = oligo_markers)
ln <- gene_lognorm(x, genes, unlist(sets))$lognorm
sample <- paste0("S", meta$Sample)
condition <- ifelse(meta$Condition == "Ctrl", "Control", "MS")

# per-cell mean expression of each marker set, then per-sample means (sets and single genes)
rows <- c(lapply(names(sets), function(s) list(name = s, v = colMeans(ln[sets[[s]], , drop = FALSE]))),
          lapply(unlist(sets), function(g) list(name = g, v = ln[g, ])))
means <- do.call(rbind, lapply(rows, function(r) {
  m <- tapply(r$v, sample, mean)
  data.frame(marker = r$name, sample = names(m), condition = condition[match(names(m), sample)], mean_lognorm = as.numeric(m))
}))
write.csv(means, file.path(out_dir, "marker_sample_means.csv"), row.names = FALSE)

# MS vs Control at the sample level: Wilcoxon and exact permutation on sample means
tests <- do.call(rbind, lapply(split(means, means$marker), function(d) {
  cbind(marker = d$marker[1], n_control = sum(d$condition == "Control"), n_ms = sum(d$condition == "MS"),
        mean_control = mean(d$mean_lognorm[d$condition == "Control"]), mean_ms = mean(d$mean_lognorm[d$condition == "MS"]),
        wilcoxon_p = wilcox.test(d$mean_lognorm[d$condition == "MS"], d$mean_lognorm[d$condition == "Control"])$p.value,
        perm_test(d$mean_lognorm, d$condition == "MS"))
}))
rownames(tests) <- NULL
write.csv(tests, file.path(out_dir, "marker_condition_tests.csv"), row.names = FALSE)
print(tests, digits = 3)
cat("positive diff = higher in MS; unit = sample (retina); tests are supportive and uncorrected across markers\n")
