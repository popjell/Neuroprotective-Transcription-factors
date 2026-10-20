# Schirmer 2019 (all cell types): donor-level pseudobulk and DESeq2, MS vs control, one fit per cell type.
# Writes results_schirmer/<cell type>/deseq2.csv in the same layout as the Trobisch results.
library(Matrix)
library(DESeq2)

set.seed(1010)
min_cells <- 5
min_donors <- 3
out <- "results_schirmer"
dir.create(out, showWarnings = FALSE)

meta <- read.delim("../data/Schirmer_2019/meta.tsv")
barcodes <- readLines("../data/Schirmer_2019/raw/barcodes.tsv")
genes <- read.delim("../data/Schirmer_2019/raw/genes.tsv", header = FALSE, col.names = c("ensembl", "symbol"))
stopifnot(all(barcodes %in% meta$cell))
meta <- meta[match(barcodes, meta$cell), ]
counts <- as(readMM("../data/Schirmer_2019/raw/matrix.mtx"), "CsparseMatrix")
stopifnot(ncol(counts) == nrow(meta), nrow(counts) == nrow(genes))

# Samples from the same donor are pooled (e.g. MS6 and MS7 are both donor MS371)
donors <- read.csv("../data/Schirmer_2019/sample_donor.csv")
stopifnot(all(meta$sample %in% donors$sample))
meta$donor <- donors$donor[match(meta$sample, donors$sample)]
meta$condition <- ifelse(grepl("^MS", meta$sample), "MS", "Control")
stopifnot(!any(tapply(meta$condition, meta$donor, function(x) length(unique(x)) > 1)))

# One column per donor x cell type; genes sharing a symbol are summed
group <- paste(meta$cell_type, meta$donor, sep = "|")
indicator <- sparseMatrix(i = seq_along(group), j = as.integer(factor(group)), x = 1, dimnames = list(NULL, levels(factor(group))))
pseudobulk <- rowsum(as.matrix(counts %*% indicator), genes$symbol)
n_cells <- table(group)[colnames(pseudobulk)]
rm(counts)

samples <- data.frame(group = colnames(pseudobulk), n_cells = as.integer(n_cells))
samples$cell_type <- sub("\\|.*", "", samples$group)
samples$donor <- sub(".*\\|", "", samples$group)
samples$condition <- factor(meta$condition[match(samples$donor, meta$donor)], levels = c("Control", "MS"))

fits <- do.call(rbind, lapply(unique(samples$cell_type), function(ct) {
  s <- samples[samples$cell_type == ct & samples$n_cells >= min_cells, ]
  n <- table(s$condition)
  row <- data.frame(cell_type = ct, control_donors = n[["Control"]], ms_donors = n[["MS"]], n_cells = sum(s$n_cells))
  row$status <- if (min(n) >= min_donors) "fitted" else "skipped"
  if (row$status == "fitted") {
    mat <- pseudobulk[, s$group]
    mat <- mat[rowSums(mat >= 10) >= min_donors, ]
    dds <- DESeq(DESeqDataSetFromMatrix(mat, s, ~condition), quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
    res$gene <- rownames(res)
    dir.create(file.path(out, ct), showWarnings = FALSE)
    write.csv(res, file.path(out, ct, "deseq2.csv"), row.names = FALSE)
  }
  cat(ct, row$status, row$control_donors, "vs", row$ms_donors, "\n")
  row
}))
write.csv(fits, file.path(out, "fit_status.csv"), row.names = FALSE)
write.csv(samples, file.path(out, "pseudobulk_samples.csv"), row.names = FALSE)
saveRDS(pseudobulk, file.path(out, "pseudobulk_counts.rds"))
