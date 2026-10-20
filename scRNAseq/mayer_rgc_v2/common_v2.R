# Shared setup for the mayer_rgc_v2 follow-up. Run scripts from inside scRNAseq/mayer_rgc_v2/. Reads the baseline data and results read-only.
# Env: OUT_ROOT (default results) is where every output of a run goes.
# Header for every output: unit = retina (Sample); 10 control, 12 MS; two MS donors contribute two retinas each; condition is confounded with tissue source (JHU controls, Netherlands Brain Bank MS).
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
set.seed(1010)
base_dir <- "../mayer_rgc"
in_dir <- file.path(base_dir, "_input")      # read-only
base_res <- file.path(base_dir, "results")   # read-only
.libPaths(c("_input/Rlib", file.path(in_dir, "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(Matrix))

out_root <- Sys.getenv("OUT_ROOT", "results")
dir.create(out_root, showWarnings = FALSE)
min_cells <- 20        # cells per unit (retina x label)
min_per_cond <- 3      # retinas per condition
resilient_groups <- c("PG_ON", "RGC26", "RGC13")     # author-defined, 27-label column
susceptible_groups <- c("MG_ON4", "MG_OFF2", "ipRGC")
lysosomal <- c("LAMP1", "LAMP2", "CTSB", "CTSD", "CTSZ", "GBA1", "HEXA", "HEXB", "GLB1", "NPC1", "NPC2", "GRN", "SCARB2", "CD63", "TMEM106B", "ASAH1", "SMPD1", "GALC", "LIPA")

# stop if an output already exists, so a rerun cannot silently replace it
no_overwrite <- function(path) {
  if (file.exists(path)) stop("output exists, remove it yourself or change OUT_ROOT: ", path)
  invisible(path)
}

# one line to console and to run.log
log_line <- function(...) {
  s <- paste0(...)
  cat(s, "\n")
  cat(s, "\n", file = file.path(out_root, "run.log"), append = TRUE)
}

# every hypothesis test goes through here; one row per test in tests_registry.csv
register_test <- function(test_id, tier, analysis, gene_or_set, estimate, se, p, n_control, n_ms, notes = "") {
  f <- file.path(out_root, "tests_registry.csv")
  row <- data.frame(test_id = test_id, tier = tier, analysis = analysis, gene_or_set = gene_or_set, estimate = estimate, se = se, p = p,
                    n_control = n_control, n_ms = n_ms, notes = notes, stringsAsFactors = FALSE)
  stopifnot(tier %in% c(1, 2, 3))
  if (file.exists(f)) {
    old <- read.csv(f, stringsAsFactors = FALSE)
    old <- old[old$test_id == test_id, ]
    if (nrow(old) == 1) {
      # a rerun may re-register a test only if the result is identical; otherwise stop
      same <- isTRUE(all.equal(c(old$estimate, old$se, old$p), c(estimate, se, p), tolerance = 1e-8))
      stopifnot(same)
      return(invisible(row))
    }
    stopifnot(nrow(old) == 0)
  }
  write.table(row, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f), qmethod = "double")
  invisible(row)
}

norm_label <- function(v) toupper(gsub("-", "_", v))

# raw counts (genes x cells) from the baseline cache; read only
load_counts_ro <- function() {
  x <- readRDS(file.path(in_dir, "counts_dgc.rds"))
  list(x = x, genes = attr(x, "genes")[[2]], cells = attr(x, "cells")[[1]])
}

# sum raw counts per unit; unit is a factor over cells (NA = skip)
sum_counts <- function(x, unit) {
  keep <- !is.na(unit)
  f <- droplevels(unit[keep])
  ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  pb <- as.matrix(x[, keep] %*% ind)
  colnames(pb) <- levels(f)
  pb
}

read_gene_list <- function(matrix_genes) {
  rk <- read.csv(file.path(in_dir, "rank_genes.csv"), stringsAsFactors = FALSE)
  g <- unique(rk$human_symbol[rk$rank == "Rank 2" & rk$use_primary_analysis == "Y"])
  stopifnot(length(g) == 31)
  if ("TBC13" %in% g && "TBC1D13" %in% matrix_genes) g[g == "TBC13"] <- "TBC1D13"
  g
}

# exact permutation of condition labels over units (all labelings if <= 1e6, else 100,000 draws); returns two-sided p for the mean difference
perm_p <- function(y, is_ms) {
  n <- length(y); k <- sum(is_ms)
  obs <- mean(y[is_ms]) - mean(y[!is_ms])
  if (choose(n, k) <= 1e6) {
    s <- colSums(matrix(y[combn(n, k)], nrow = k))
    d <- s / k - (sum(y) - s) / (n - k)
    c(p = mean(abs(d) >= abs(obs) - 1e-12), mode = "exact", obs = obs)
  } else {
    d <- replicate(1e5, { i <- sample(n, k); mean(y[i]) - mean(y[-i]) })
    c(p = (1 + sum(abs(d) >= abs(obs) - 1e-12)) / (1e5 + 1), mode = "100000 draws", obs = obs)
  }
}
