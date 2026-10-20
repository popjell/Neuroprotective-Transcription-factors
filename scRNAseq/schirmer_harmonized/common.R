# Shared setup for the Schirmer harmonized-PSAP scripts. Run them from scRNAseq/schirmer_harmonized/ (see run_all.sh).
# Schirmer neurons are relabeled with Trobisch subtypes (SingleR), pooled into upper_EN / deep_EN / IN, and tested for PSAP, MS vs Control.
.libPaths(c("Rlib", .libPaths()))  # SingleR and scuttle live here
library(Matrix)
library(DESeq2)

set.seed(1010)
gene <- "PSAP"
min_cells <- 20    # cells per donor x group
min_donors <- 3    # control donors and MS donors per test
min_count <- 10    # gene filter: counts >= min_count in at least min_samples pseudobulks
min_samples <- 2
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE)
res <- function(...) file.path(res_dir, ...)

schirmer_dir <- "../data/Schirmer_2019"
trobisch_dir <- "../ms-cross-regional/cor"
rank_genes <- "../../rank_genes.csv"

# Trobisch cortex sample -> Schirmer donor ID (from cell_death/00_prepare_cortex.R); these 7 donors are in both datasets
trobisch_donor <- c(CR_ctrl1 = "CO39", CR_ctrl2 = "CO14", CR_ctrl3 = "CO36", CR_ctrl4 = "CO25",
                    CR_ms1 = "MS60", CR_ms2 = "MS229", CR_ms3 = "MS371", CR_ms4 = "MS371")
shared_donors <- unique(unname(trobisch_donor))

# Coarse group from the Trobisch label name only: upper_EN (first layer token L2 or L3), deep_EN (other -EN), IN (-IN), else NA
coarse_group <- function(label) {
  ifelse(grepl("-EN$", label), ifelse(sub("-.*", "", label) %in% c("L2", "L3"), "upper_EN", "deep_EN"),
         ifelse(grepl("-IN$", label), "IN", NA_character_))
}

no_overwrite <- function(file) if (file.exists(file)) stop(file, " exists; move it aside first")

# Sum counts (genes x cells) over cells sharing a group label; returns genes x groups
pseudobulk <- function(counts, group) {
  g <- factor(group)
  indicator <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1, dimnames = list(NULL, levels(g)))
  as.matrix(counts %*% indicator)
}

# DESeq2 ~condition (Control reference), MS vs Control, on genes with min_count in min_samples pseudobulks. s: one row per column of mat.
fit_deseq <- function(mat, s) {
  mat <- mat[rowSums(mat >= min_count) >= min_samples, ]
  dds <- DESeq(DESeqDataSetFromMatrix(mat, s, ~condition), quiet = TRUE)
  out <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
  out$gene <- rownames(out)
  list(res = out, sizes = sizeFactors(dds))
}

# PSAP row of a result table with the 95% interval (lfcSE x 1.96)
psap_row <- function(r, ...) {
  x <- r[r$gene == gene, ]
  stopifnot(nrow(x) == 1)
  data.frame(..., log2FC = x$log2FoldChange, lower = x$log2FoldChange - 1.96 * x$lfcSE, upper = x$log2FoldChange + 1.96 * x$lfcSE, p = x$pvalue)
}
