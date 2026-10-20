# Part B: decontX on the full cortex object, contamination summaries, DE rerun on decontaminated neuron counts, PSAP survival, marker rule rerun.
# Needs data/cortex.rds, ../cell_death/data/cortex_neurons.rds and ../results/Cortex/<subtype>/deseq2.csv (read only). celda is in ./Rlib.
# decontX: default parameters, counts = raw integers of all 27,577 cortex nuclei, z = broad type (see broad_of in common.R), batch = sample_id.
source("common.R")
library(celda)

full <- readRDS("data/cortex.rds")
cells <- full$cells
cat("cells:", nrow(cells), " broad types:", paste(unique(cells$broad), collapse = ", "), "\n")

out <- "data/decontx.rds"
no_overwrite(out)
dx <- decontX(full$counts, z = cells$broad, batch = cells$sample_id)
stopifnot(identical(colnames(dx$decontXcounts), cells$cell), identical(rownames(dx$decontXcounts), rownames(full$counts)))
dec <- dx$decontXcounts
dec@x <- round(dec@x)
dec <- drop0(dec)
saveRDS(list(contamination = setNames(dx$contamination, cells$cell), counts = dec), out)
cells$contamination <- dx$contamination
rm(dx, full)
gc()

# Contamination per donor and broad type (median and IQR), then neurons by condition with per-donor medians side by side
iqr <- function(v) c(median = median(v), q25 = unname(quantile(v, 0.25)), q75 = unname(quantile(v, 0.75)))
by_donor <- do.call(rbind, lapply(split(cells, list(cells$donor, cells$broad), drop = TRUE), function(d)
  data.frame(donor = d$donor[1], condition = d$condition[1], broad = d$broad[1], n_cells = nrow(d), t(iqr(d$contamination)))))
rownames(by_donor) <- NULL
write.csv(by_donor, res("decontx_contamination_by_donor.csv"), row.names = FALSE)
by_broad <- do.call(rbind, lapply(split(cells, cells$broad), function(d) data.frame(broad = d$broad[1], n_cells = nrow(d), t(iqr(d$contamination)))))
write.csv(by_broad, res("decontx_contamination_by_broad.csv"), row.names = FALSE)
neuron_donor <- by_donor[by_donor$broad == "Neuron", ]
neuron_donor <- neuron_donor[order(neuron_donor$condition, neuron_donor$donor), ]
write.csv(neuron_donor, res("decontx_neuron_by_condition.csv"), row.names = FALSE)
cat("contamination by broad type\n"); print(by_broad, row.names = FALSE)
cat("neuron contamination, per-donor medians\n"); print(neuron_donor, row.names = FALSE)
cat("median of donor medians: Control", median(neuron_donor$median[neuron_donor$condition == "Control"]),
    " MS", median(neuron_donor$median[neuron_donor$condition == "MS"]), "\n")

# Rerun the pipeline on the decontaminated counts of the neurons in cortex_neurons.rds
neurons <- readRDS(neurons_rds)
stopifnot(all(neurons$cells$cell %in% colnames(dec)))
fit <- run_de(dec[, neurons$cells$cell], neurons$cells)
fitted <- fit$status$subtype[fit$status$fitted]
write.csv(fit$status, res("decontx_fit_status.csv"), row.names = FALSE)
for (s in fitted) {
  dir.create(file.path(res_dir, "decontx", s), recursive = TRUE, showWarnings = FALSE)
  write.csv(fit$de[fit$de$subtype == s, ], file.path(res_dir, "decontx", s, "deseq2.csv"), row.names = FALSE)
}

# PSAP: survives if decontaminated padj < 0.05 and log2FC >= 0.5 x the original log2FC
psap <- do.call(rbind, lapply(read.csv("../results/fit_status.csv")$subtype[read.csv("../results/fit_status.csv")$status == "fitted" &
                                                                                 read.csv("../results/fit_status.csv")$region == "Cortex"], function(s) {
  o <- read.csv(file.path(baseline_dir, s, "deseq2.csv"))
  o <- o[o$gene == "PSAP", ]
  d <- fit$de[fit$de$subtype == s & fit$de$gene == "PSAP", ]
  data.frame(subtype = s, original_log2FC = o$log2FoldChange, original_padj = o$padj, original_significant = o$padj < 0.05,
             decontx_log2FC = if (nrow(d)) d$log2FoldChange else NA, decontx_padj = if (nrow(d)) d$padj else NA)
}))
psap$survives <- !is.na(psap$decontx_padj) & psap$decontx_padj < 0.05 & psap$decontx_log2FC >= 0.5 * psap$original_log2FC
write.csv(psap, res("psap_decontx_comparison.csv"), row.names = FALSE)
print(psap, row.names = FALSE)
cat("PSAP survives in", sum(psap$survives & psap$original_significant), "of", sum(psap$original_significant), "originally significant subtypes\n")

# Part A rule on the decontaminated results
mt <- marker_table(fit$de)
write.csv(mt, res("marker_negative_controls_decontx.csv"), row.names = FALSE)
flags <- marker_flags(mt)
write.csv(flags, res("marker_flags_decontx.csv"), row.names = FALSE)
print(flags, row.names = FALSE)
cat("subtype x glial type combinations flagged after decontX:", sum(flags$flag), "of", nrow(flags), "\n")

writeLines(capture.output(sessionInfo()), res("sessionInfo.txt"))
