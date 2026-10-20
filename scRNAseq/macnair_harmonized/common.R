# Shared setup for the Macnair harmonized-PSAP scripts. Run them from scRNAseq/macnair_harmonized/ (see run_all.sh).
# Macnair grey matter neurons are relabeled with Trobisch subtypes (SingleR), pooled into upper_EN / deep_EN / IN, then tested for PSAP.
# Environment variable: MACNAIR_DIR (default "data")
.libPaths(c("../schirmer_harmonized/Rlib", .libPaths()))  # SingleR, scuttle, scran live here
library(Matrix)

set.seed(1010)
gene <- "PSAP"
min_cells <- 20      # cells per donor (or sample) x group
min_donors <- 3      # control donors and MS donors per test; also GML/NAGM donor pairs
min_count <- 10      # gene filter: counts >= min_count in at least min_samples pseudobulks
min_samples <- 2
chunk_records <- 5e7 # mtx records per read
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE)
res <- function(...) file.path(res_dir, ...)

macnair_dir <- Sys.getenv("MACNAIR_DIR", "data")
trobisch_dir <- "../ms-cross-regional/cor"
schirmer_donor_file <- "../data/Schirmer_2019/sample_donor.csv"
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

# The 31 Rank 2 genes (rank_genes.csv, Rank 2 and use_primary_analysis == Y; TBC13 -> TBC1D13)
rank2_genes <- function() {
  rg <- read.csv(rank_genes)
  rg <- rg[rg$rank == "Rank 2" & rg$use_primary_analysis == "Y", ]
  genes <- sub("^TBC13$", "TBC1D13", rg$human_symbol)
  stopifnot(length(genes) == 31)
  genes
}

# Sum counts (genes x cells) over cells sharing a group label; returns genes x groups
pseudobulk <- function(counts, group) {
  g <- factor(group)
  indicator <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1, dimnames = list(NULL, levels(g)))
  as.matrix(counts %*% indicator)
}

# Cached sample x group pseudobulks (genes x "sample|group"): counts plus one row of annotation per column
load_pb <- function() {
  pb <- readRDS(res("pseudobulk_sample_group.rds"))
  s <- read.csv(res("samples.csv"))
  id <- colnames(pb$counts)
  cols <- data.frame(id = id, sample = sub("[|].*", "", id), group = sub(".*[|]", "", id))
  cols$n_cells <- pb$n_cells[cols$id]
  cols[, c("donor", "condition", "tissue")] <- s[match(cols$sample, s$sample), c("donor", "condition", "tissue")]
  stopifnot(!anyNA(cols$donor))
  list(counts = pb$counts, cols = cols)
}

# DESeq2 on genes with min_count in min_samples pseudobulks. s: one row per column of mat. Returns the results table, size factors and the fit.
fit_deseq <- function(mat, s, design, contrast) {
  library(DESeq2)
  mat <- mat[rowSums(mat >= min_count) >= min_samples, ]
  storage.mode(mat) <- "integer"
  dds <- DESeq(DESeqDataSetFromMatrix(mat, s, design), quiet = TRUE)
  out <- as.data.frame(results(dds, contrast = contrast))
  out$gene <- rownames(out)
  list(res = out, sizes = sizeFactors(dds), dds = dds)
}

# PSAP row of a DESeq2 result table with the 95% interval (lfcSE x 1.96)
psap_row <- function(r, ...) {
  x <- r[r$gene == gene, ]
  stopifnot(nrow(x) == 1)
  data.frame(..., log2FC = x$log2FoldChange, lower = x$log2FoldChange - 1.96 * x$lfcSE, upper = x$log2FoldChange + 1.96 * x$lfcSE, p = x$pvalue)
}

# limma-voom with donor blocking (duplicateCorrelation, two rounds). y: DGEList already filtered and normalized.
voom_block <- function(y, design, block) {
  library(limma)
  v <- voom(y, design)
  cor1 <- duplicateCorrelation(v, design, block = block)$consensus.correlation
  v <- voom(y, design, block = block, correlation = cor1)
  cor2 <- duplicateCorrelation(v, design, block = block)$consensus.correlation
  list(fit = lmFit(v, design, block = block, correlation = cor2), cor = cor2)
}
