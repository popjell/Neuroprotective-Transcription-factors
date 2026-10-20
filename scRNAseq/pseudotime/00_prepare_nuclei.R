# EN-L2-3 nuclei that have an author pseudotime in Schirmer 2019: metadata and a sparse count matrix (one row per gene symbol)
library(Matrix)
library(data.table)
source("common.R")

meta <- read.delim(file.path(data_dir, "meta.tsv"))
barcodes <- readLines(file.path(data_dir, "raw/barcodes.tsv"))
genes <- read.delim(file.path(data_dir, "raw/genes.tsv"), header = FALSE, col.names = c("ensembl", "symbol"))
stopifnot(identical(meta$cell, barcodes))

# Pooled donors: samples MS6 + MS7 are MS371 and MS2 + MS8 are MS377
donors <- read.csv(file.path(data_dir, "sample_donor.csv"))
meta$donor <- donors$donor[match(meta$sample, donors$sample)]
stopifnot(!anyNA(meta$donor))

keep <- which(meta$cell_type %in% cell_types & is.finite(meta$pseudotime))
nuclei <- data.frame(nucleus = meta$cell[keep], sample = meta$sample[keep], donor = meta$donor[keep], condition = meta$diagnosis[keep],
                     Seqbatch = meta$Seqbatch[keep], UMIs = meta$UMIs[keep], author_pseudotime = meta$pseudotime[keep])
write.csv(nuclei, shared("nuclei.csv"), row.names = FALSE)
cat("nuclei:", nrow(nuclei), " donors with at least", min_cells, "nuclei:", sum(table(nuclei$donor) >= min_cells), "\n")

# Only the chosen nuclei are read from the 1.4 GB matrix; genes sharing a symbol are summed
entries <- fread(file.path(data_dir, "raw/matrix.mtx"), skip = 2, header = FALSE, col.names = c("gene", "cell", "count"))
entries <- entries[cell %in% keep]
counts <- sparseMatrix(i = entries$gene, j = match(entries$cell, keep), x = as.numeric(entries$count),
                       dims = c(nrow(genes), length(keep)), dimnames = list(NULL, nuclei$nucleus))
symbols <- unique(genes$symbol)
counts <- sparseMatrix(i = match(genes$symbol, symbols), j = seq_len(nrow(genes)), x = 1) %*% counts
rownames(counts) <- symbols
saveRDS(counts, shared("counts.rds"))
cat("genes:", nrow(counts), "\n")
