# Check whether the granule-neuron results change after accounting for age or sex.
library(DESeq2)

x <- readRDS("results/pseudobulk_counts.rds")
info <- x$metadata
keep <- info$subtype == "GRAN-EN" & info$cells >= 20
counts <- as.matrix(x$counts[, keep])
counts <- counts[rowSums(counts >= 10) >= 2, ]
info <- info[keep, ]
cells <- readRDS("data/new_regions.rds")$meta
info$age <- cells$age[match(info$sample_id, cells$sample_id)]
info$age <- (info$age - mean(info$age)) / 10
info$sex <- factor(cells$sex[match(info$sample_id, cells$sample_id)])
info$condition <- factor(info$condition, levels = c("Control", "MS"))

out <- list()
for (model in c("~condition", "~age+condition", "~sex+condition")) {
  dds <- DESeqDataSetFromMatrix(counts, info, design = as.formula(model))
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
  result$gene <- rownames(result)
  result$model <- model
  out[[model]] <- result
}
write.csv(do.call(rbind, out), "results/granule_sensitivity.csv", row.names = FALSE)
