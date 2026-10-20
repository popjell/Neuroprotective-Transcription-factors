# Validation: normalized_genes.csv and size factors against an independent recomputation of the pooled pseudobulk; writes sessionInfo.txt.
source("common.R")
suppressPackageStartupMessages(library(DESeq2))

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
sample <- paste0("S", meta$Sample)

# independent pseudobulk from the matrix: loop over samples, no sparse indicator trick
ids <- sort(unique(sample))
pb <- sapply(ids, function(s) as.numeric(rowSums(x[, sample == s])))
rownames(pb) <- genes
keep_g <- rowSums(pb >= min_count) >= min_pseudobulks
cat("genes passing filter (independent):", sum(keep_g), "\n")

donors <- read.csv(file.path(res_dir, "all_RGC", "donors.csv"))
norm <- read.csv(file.path(res_dir, "all_RGC", "normalized_genes.csv"))
sf_indep <- sizeFactors(estimateSizeFactors(DESeqDataSetFromMatrix(pb[keep_g, ], data.frame(row.names = ids, one = rep(1, length(ids))), ~1)))
sf_rep <- donors$size_factor[match(ids, donors$sample)]
cat("size factors max abs difference:", max(abs(sf_indep - sf_rep)), "\n")
stopifnot(max(abs(sf_indep - sf_rep)) < 1e-8)

# normalized_count * size_factor must give back the raw summed count
chk <- pb[cbind(match(norm$gene, genes), match(norm$sample, ids))]
back <- norm$normalized_count * donors$size_factor[match(norm$sample, donors$sample)]
cat("normalized_genes rows:", nrow(norm), " max abs difference to raw counts:", max(abs(back - chk)), "\n")
stopifnot(max(abs(back - chk)) < 1e-6)
cat("per-sample count sums of gene-list genes match:", all.equal(tapply(chk, norm$sample, sum), tapply(back, norm$sample, sum)), "\n")
cat("n_cells in donors.csv matches metadata:", all(donors$n_cells == as.integer(table(sample)[donors$sample])), "\n")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("validation passed\n")
