# Shared setup for the Mayer 2026 RGC replication of the Trobisch 31-gene Rank 2 analysis. Run scripts from inside scRNAseq/mayer_rgc/.
# Run with: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript <script>.R
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
set.seed(1010)
in_dir <- "_input"
.libPaths(c(file.path(in_dir, "Rlib"), .libPaths()))  # lmerTest and UCell were installed here (system library not writable)
suppressPackageStartupMessages(library(Matrix))

res_dir <- "results"
min_cells <- 20        # cells per donor per group being fitted
min_donors <- 3        # donors per condition to fit
min_count <- 10        # gene filter: counts >= this ...
min_pseudobulks <- 2   # ... in at least this many pseudobulks
primary_gene <- "PSAP"
secondary_gene <- "ITM2C"

# stop if an output already exists, so a rerun cannot silently replace it
no_overwrite <- function(path) {
  if (file.exists(path)) stop("output exists, remove it yourself or change the name: ", path)
  invisible(path)
}

microglia_markers <- c("P2RY12", "CSF1R", "C1QB", "CD74")
astro_markers <- c("GFAP", "GLUL", "RLBP1")
oligo_markers <- c("MBP", "PLP1")
trim_quantile <- 0.95  # contamination trim: drop cells above this quantile of the microglia UCell score

# the 31-gene Rank 2 list; TBC13 is a non-standard symbol in rank_genes.csv, mapped to TBC1D13 when present in the matrix
read_gene_list <- function(matrix_genes) {
  rk <- read.csv(file.path(in_dir, "rank_genes.csv"), stringsAsFactors = FALSE)
  g <- unique(rk$human_symbol[rk$rank == "Rank 2" & rk$use_primary_analysis == "Y"])
  stopifnot(length(g) == 31)
  if ("TBC13" %in% g && "TBC1D13" %in% matrix_genes) g[g == "TBC13"] <- "TBC1D13"
  g
}

# dense log1p(counts per 10,000) for a few genes; x is the sparse genes x cells matrix
gene_lognorm <- function(x, matrix_genes, g) {
  g <- g[g %in% matrix_genes]
  cnt <- as.matrix(x[match(g, matrix_genes), , drop = FALSE])
  rownames(cnt) <- g
  list(counts = cnt, lognorm = log1p(t(t(cnt) / colSums(x) * 1e4)))
}

# exact permutation test on the difference of group means (all labelings if feasible, else 10,000 random)
perm_test <- function(means, is_ms) {
  n <- length(means); k <- sum(is_ms)
  obs <- mean(means[is_ms]) - mean(means[!is_ms])
  if (choose(n, k) <= 1e6) {
    idx <- combn(n, k)
    s <- colSums(matrix(means[idx], nrow = k))
    d <- s / k - (sum(means) - s) / (n - k)
    p <- mean(abs(d) >= abs(obs) - 1e-12)
    mode <- "exact"
  } else {
    d <- replicate(10000, { i <- sample(n, k); mean(means[i]) - mean(means[-i]) })
    p <- (1 + sum(abs(d) >= abs(obs) - 1e-12)) / 10001
    mode <- "10000 random"
  }
  data.frame(diff_ms_minus_ctrl = obs, perm_p = p, perm_mode = mode)
}

# read the 10x-style mtx.gz once, cache a sparse genes x cells matrix as rds
load_counts <- function() {
  cache <- file.path(in_dir, "counts_dgc.rds")
  if (file.exists(cache)) return(readRDS(cache))
  library(data.table)
  genes <- fread(file.path(in_dir, "RGC_RNA_features.tsv.gz"), header = FALSE)
  cells <- fread(file.path(in_dir, "RGC_RNA_barcodes.tsv.gz"), header = FALSE)
  m <- fread(file.path(in_dir, "RGC_RNA_matrix.mtx.gz"), skip = 2, header = FALSE)  # skip %% comment + dims line
  dims <- scan(file.path(in_dir, "RGC_RNA_matrix.mtx.gz"), what = "", n = 1, skip = 1, sep = "\n", quiet = TRUE)
  dims <- as.integer(strsplit(dims, " ")[[1]])
  x <- sparseMatrix(i = m[[1]], j = m[[2]], x = m[[3]], dims = dims[1:2])
  rm(m); gc()
  stopifnot(nrow(x) == nrow(genes), ncol(x) == nrow(cells))
  attr(x, "genes") <- genes
  attr(x, "cells") <- cells
  saveRDS(x, cache)
  x
}
