# POST-HOC, exploratory: limma roast/fry on the 31-gene set, MS vs Control, retina-level pseudobulk (voom). Not part of the frozen plan.
# Same eligibility and gene filter as mayer_deseq2.R. Writes results/roast/roast.csv.
source("common.R")
suppressPackageStartupMessages({ library(limma); library(edgeR) })

n_rot <- 9999  # rotations for roast
out_dir <- file.path(res_dir, "roast")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "roast.csv"))

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
stopifnot(identical(attr(x, "cells")[[1]], meta$cell_barcode))
meta$sample <- paste0("S", meta$Sample)
meta$condition <- ifelse(meta$Condition == "Ctrl", "Control", "MS")
g31 <- read_gene_list(genes)

out <- list()
for (grp in c("all_RGC", sort(unique(meta$celltypes)))) {
  cells <- if (grp == "all_RGC") rep(TRUE, nrow(meta)) else meta$celltypes == grp
  n_cells <- table(factor(meta$sample[cells], levels = sort(unique(meta$sample))))
  ok <- names(n_cells)[n_cells >= min_cells]
  cond <- meta$condition[match(ok, meta$sample)]
  if (sum(cond == "Control") < min_donors || sum(cond == "MS") < min_donors) next

  # pseudobulk, same filter as DESeq2, TMM + voom
  use <- cells & meta$sample %in% ok
  f <- factor(meta$sample[use], levels = ok)
  pb <- as.matrix(x[, use] %*% sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), length(ok))))
  dimnames(pb) <- list(genes, ok)
  pb <- pb[rowSums(pb >= min_count) >= min_pseudobulks, ]
  design <- model.matrix(~ factor(cond, levels = c("Control", "MS")))
  v <- voom(calcNormFactors(DGEList(pb)), design)

  # set tests; roast is rotation-based (seeded), fry is its fast exact approximation
  sets <- list(set31 = intersect(g31, rownames(pb)), set31_no_psap = setdiff(intersect(g31, rownames(pb)), primary_gene))
  for (nm in names(sets)) {
    idx <- sets[[nm]]
    set.seed(1010)
    r <- roast(v, index = idx, design = design, contrast = 2, nrot = n_rot)
    p <- r$p.value  # rows Down, Up, UpOrDown, Mixed
    fr <- fry(v, index = idx, design = design, contrast = 2)
    out[[paste(grp, nm)]] <- data.frame(group = grp, set = nm, n_control = sum(cond == "Control"), n_ms = sum(cond == "MS"), n_genes = length(idx),
                                        roast_p_up = p["Up", "P.Value"], roast_p_down = p["Down", "P.Value"],
                                        roast_p_two_sided = p["UpOrDown", "P.Value"], roast_p_mixed = p["Mixed", "P.Value"],
                                        fry_direction = fr$Direction, fry_p = fr$PValue, fry_p_mixed = fr$PValue.Mixed)
  }
}
res <- do.call(rbind, out)
rownames(res) <- NULL
write.csv(res, file.path(out_dir, "roast.csv"), row.names = FALSE)
print(res[, c("group", "set", "n_control", "n_ms", "n_genes", "roast_p_up", "roast_p_two_sided", "roast_p_mixed", "fry_direction", "fry_p")], digits = 3)
cat("positive diff = set higher in MS (direction Up); post-hoc, exploratory; no correction across groups\n")
