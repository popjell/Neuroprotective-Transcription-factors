# Rebuild the donor-level pseudobulk DESeq2 (MS vs Control) on uncorrected neuron counts and check it against the existing ../results/Cortex.
# Needs data/cortex.rds, ../cell_death/data/cortex_neurons.rds and ../results/Cortex/<subtype>/deseq2.csv (read only).
source("common.R")

neurons <- readRDS(neurons_rds)
full <- readRDS("data/cortex.rds")

# The neuron cache must be the same counts as the full object
i <- match(neurons$cells$cell, full$cells$cell)
stopifnot(!anyNA(i), identical(rownames(neurons$counts), rownames(full$counts)),
          identical(neurons$counts, full$counts[, i]), identical(neurons$cells$donor, full$cells$donor[i]))

fit <- run_de(neurons$counts, neurons$cells)
fitted <- fit$status$subtype[fit$status$fitted]
write.csv(fit$status, res("rebuild_fit_status.csv"), row.names = FALSE)
for (s in fitted) {
  dir.create(file.path(res_dir, "rebuild", s), recursive = TRUE, showWarnings = FALSE)
  write.csv(fit$de[fit$de$subtype == s, ], file.path(res_dir, "rebuild", s, "deseq2.csv"), row.names = FALSE)
}

# Reproduction check against the baseline: same fitted subtypes, same genes, same log2FC and padj
base_fit <- read.csv("../results/fit_status.csv")
base_fit <- base_fit[base_fit$region == "Cortex" & base_fit$status == "fitted", "subtype"]
stopifnot(setequal(fitted, base_fit))
check <- do.call(rbind, lapply(fitted, function(s) {
  b <- read.csv(file.path(baseline_dir, s, "deseq2.csv"))
  r <- fit$de[fit$de$subtype == s, ]
  j <- match(b$gene, r$gene)
  data.frame(subtype = s, genes_baseline = nrow(b), genes_rebuild = nrow(r), genes_shared = sum(!is.na(j)),
             max_abs_diff_log2FC = max(abs(b$log2FoldChange - r$log2FoldChange[j]), na.rm = TRUE),
             max_abs_diff_padj = max(abs(b$padj - r$padj[j]), na.rm = TRUE),
             psap_log2FC_baseline = b$log2FoldChange[b$gene == "PSAP"], psap_log2FC_rebuild = r$log2FoldChange[r$gene == "PSAP"],
             psap_padj_baseline = b$padj[b$gene == "PSAP"], psap_padj_rebuild = r$padj[r$gene == "PSAP"])
}))
write.csv(check, res("reproduction_check.csv"), row.names = FALSE)
print(check, digits = 6)
stopifnot(all(check$genes_baseline == check$genes_rebuild), all(check$genes_shared == check$genes_baseline),
          abs(check$psap_log2FC_baseline - check$psap_log2FC_rebuild) < 1e-6, abs(check$psap_padj_baseline - check$psap_padj_rebuild) < 1e-6)
cat("reproduction check passed:", length(fitted), "subtypes; PSAP padj < 0.05 in", sum(check$psap_padj_rebuild < 0.05), "\n")
