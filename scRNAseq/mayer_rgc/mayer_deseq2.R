# Retina-level pseudobulk DESeq2, MS vs Control, for all RGCs pooled and per RGC subtype; reports the 31-gene Rank 2 list.
# Unit = Sample (one retina; 2 MS donors have 2 retinas, mapping not in the data). Env: OUT_ROOT (default results), TRIM=1 drops top 5% microglia-score cells.
source("common.R")
suppressPackageStartupMessages(library(DESeq2))

out_root <- Sys.getenv("OUT_ROOT", res_dir)
trim <- Sys.getenv("TRIM", "0") == "1"
trim_quantile <- 0.95  # drop cells above this quantile of the microglia UCell score
dir.create(out_root, showWarnings = FALSE)
log_file <- file.path(out_root, if (trim) "run_deseq2_trimmed.log" else "run_deseq2.log")
no_overwrite(log_file)
sink(log_file, split = TRUE)

# load
meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
stopifnot(identical(attr(x, "cells")[[1]], meta$cell_barcode), anyDuplicated(genes) == 0)
meta$sample <- paste0("S", meta$Sample)
meta$condition <- ifelse(meta$Condition == "Ctrl", "Control", "MS")
stopifnot(all(tapply(meta$condition, meta$sample, function(v) length(unique(v))) == 1))  # one condition per sample
cat("cells:", nrow(meta), " samples:", length(unique(meta$sample)), "\n")

# drop microglia-contaminated cells for the trimmed sensitivity
if (trim) {
  ms <- read.csv(file.path(res_dir, "contamination", "cell_scores.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(ms$cell_barcode, meta$cell_barcode))
  keep <- ms$microglia_ucell <= quantile(ms$microglia_ucell, trim_quantile)
  cat("contamination trim: dropping", sum(!keep), "of", length(keep), "cells\n")
  meta <- meta[keep, ]
  x <- x[, keep]
}

# gene list (TBC13 is a non-standard symbol in rank_genes.csv)
rk <- read.csv(file.path(in_dir, "rank_genes.csv"), stringsAsFactors = FALSE)
g_input <- unique(rk$human_symbol[rk$rank == "Rank 2" & rk$use_primary_analysis == "Y"])
stopifnot(length(g_input) == 31)
g31 <- g_input
if ("TBC13" %in% g31 && "TBC1D13" %in% genes) {
  g31[g31 == "TBC13"] <- "TBC1D13"  # TBC13 is not an HGNC symbol; TBC1D13 is in the matrix
  cat("gene mapping: TBC13 -> TBC1D13\n")
}
cat("31-gene list:", paste(g31, collapse = " "), "\n")
cat("genes in matrix:", sum(g31 %in% genes), "of 31\n")

# groups: pooled, then each celltypes label
groups <- c("all_RGC", sort(unique(meta$celltypes)))
status <- list()
all_genes <- list()

for (grp in groups) {
  cells <- if (grp == "all_RGC") rep(TRUE, nrow(meta)) else meta$celltypes == grp
  n_cells <- table(factor(meta$sample[cells], levels = sort(unique(meta$sample))))
  ok <- names(n_cells)[n_cells >= min_cells]
  cond <- meta$condition[match(ok, meta$sample)]
  n_ctrl <- sum(cond == "Control"); n_ms <- sum(cond == "MS")
  cat("\n== ", grp, ": cells ", sum(cells), "; samples with >=", min_cells, " cells: ", n_ctrl, " Control, ", n_ms, " MS\n", sep = "")
  if (n_ctrl < min_donors || n_ms < min_donors) {
    cat("skipped: fewer than", min_donors, "samples in a condition\n")
    status[[grp]] <- data.frame(group = grp, fitted = "N", n_control = n_ctrl, n_ms = n_ms, reason = "fewer than 3 eligible samples in a condition")
    next
  }

  # pseudobulk: sum raw counts per sample
  use <- cells & meta$sample %in% ok
  f <- factor(meta$sample[use], levels = ok)
  ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), length(ok)))
  pb <- as.matrix(x[, use] %*% ind)
  dimnames(pb) <- list(genes, ok)
  stopifnot(all(pb == round(pb)), sum(pb) == sum(x[, use]))
  keep_g <- rowSums(pb >= min_count) >= min_pseudobulks
  cat("pseudobulks:", ncol(pb), " genes before filter:", nrow(pb), " after (>=", min_count, " counts in >=", min_pseudobulks, " pseudobulks): ", sum(keep_g), "\n", sep = "")

  # DESeq2: ~condition, Control is the reference, positive log2FC = higher in MS
  coldata <- data.frame(sample = ok, condition = factor(cond, levels = c("Control", "MS")), row.names = ok)
  dds <- DESeqDataSetFromMatrix(pb[keep_g, ], coldata, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "MS", "Control"))
  de <- data.frame(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                   ci_low = res$log2FoldChange - 1.96 * res$lfcSE, ci_high = res$log2FoldChange + 1.96 * res$lfcSE,
                   stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  cat("genes tested (non-NA pvalue):", sum(!is.na(de$pvalue)), " padj<0.05 across transcriptome:", sum(de$padj < 0.05, na.rm = TRUE), "\n")

  # the 31 genes, including those removed by the count filter
  gt <- data.frame(gene = g31, input_symbol = g_input, in_matrix = ifelse(g31 %in% genes, "Y", "N"))
  gt$filtered_in <- ifelse(gt$gene %in% de$gene, "Y", "N")
  gt$reason <- ifelse(gt$in_matrix == "N", "not in matrix", ifelse(gt$filtered_in == "N", "removed by count filter", ""))
  gt <- cbind(gt, de[match(gt$gene, de$gene), c("baseMean", "log2FC", "lfcSE", "ci_low", "ci_high", "pvalue", "padj")])
  gt$tested <- ifelse(!is.na(gt$pvalue), "Y", "N")
  gt$padj_within31 <- NA_real_  # SECONDARY: BH across the testable genes of the 31-gene list only
  gt$padj_within31[gt$tested == "Y"] <- p.adjust(gt$pvalue[gt$tested == "Y"], method = "BH")
  gt$group <- grp
  cat("31-gene list: in matrix", sum(gt$in_matrix == "Y"), " passed filter", sum(gt$filtered_in == "Y"), " tested", sum(gt$tested == "Y"), "\n")

  # write
  dir_g <- file.path(out_root, grp)
  dir.create(dir_g, showWarnings = FALSE)
  for (fn in c("deseq2.csv", "genes.csv", "donors.csv", "normalized_genes.csv")) no_overwrite(file.path(dir_g, fn))
  write.csv(de, file.path(dir_g, "deseq2.csv"), row.names = FALSE)
  write.csv(gt, file.path(dir_g, "genes.csv"), row.names = FALSE)
  donors <- data.frame(sample = ok, condition = cond, n_cells = as.integer(n_cells[ok]), size_factor = sizeFactors(dds), row.names = NULL)
  write.csv(donors, file.path(dir_g, "donors.csv"), row.names = FALSE)  # sample = retina; no sex/age/batch covariates exist in the data
  nc <- counts(dds, normalized = TRUE)
  gi <- intersect(g31, rownames(nc))
  norm_long <- do.call(rbind, lapply(gi, function(g) data.frame(sample = ok, condition = cond, gene = g, normalized_count = nc[g, ])))
  write.csv(norm_long, file.path(dir_g, "normalized_genes.csv"), row.names = FALSE)
  all_genes[[grp]] <- gt
  status[[grp]] <- data.frame(group = grp, fitted = "Y", n_control = n_ctrl, n_ms = n_ms, reason = "")
}

# sensitivity covariates: none exist in the metadata
cat("\nsensitivity with covariate: skipped. No sex, age or batch column in the metadata; Sample is the unit itself and condition is nested in it.\n")
cat("padj is BH across the transcriptome within each fit; no correction across subtype groups (", length(groups) - 1, " subtype fits).\n", sep = "")

write.csv(do.call(rbind, status), file.path(out_root, "fit_status.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_genes), file.path(out_root, "all_genes.csv"), row.names = FALSE)
sink()
