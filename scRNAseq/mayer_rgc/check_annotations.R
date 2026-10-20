# Check the authors' cell-type labels against marker expression before using them. Writes results/annotation_check/*.csv and a text summary.
source("common.R")

out_dir <- file.path(res_dir, "annotation_check")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "summary.txt"))

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
lib <- colSums(x)
cat("nCount_RNA equals column sums:", all(lib == meta$nCount_RNA), "\n")
cat("nFeature_RNA equals detected genes:", all(diff(x@p) == meta$nFeature_RNA), "\n")

# marker panels (human retina)
markers <- list(
  pan_rgc = c("RBPMS", "RBPMS2", "SNCG", "POU4F1", "POU4F2", "THY1", "NEFL", "SLC17A6", "RBFOX3"),
  photoreceptor = c("RHO", "RCVRN", "PDE6A", "GNAT1", "ARR3"),
  bipolar = c("VSX2", "TRPM1", "GRM6", "PRKCA"),
  amacrine = c("GAD1", "GAD2", "SLC6A9", "TFAP2A"),
  horizontal = c("ONECUT1", "LHX1"),
  microglia = c("P2RY12", "CSF1R", "C1QB", "CD74", "PTPRC"),
  astro_muller = c("GFAP", "GLUL", "RLBP1", "AQP4"),
  oligo = c("MBP", "PLP1", "MOG"),
  endothelial = c("CLDN5", "PECAM1"),
  subtype = c("OPN4", "SPP1", "TBR1", "FOXP2", "CFH", "GPR149", "EOMES", "SATB2")
)
panel <- unique(unlist(markers))
panel <- panel[panel %in% genes]
cnt <- as.matrix(x[match(panel, genes), ])
rownames(cnt) <- panel
lognorm <- log1p(t(t(cnt) / lib * 1e4))  # log1p CP10K

summ <- function(label) {
  g <- split(seq_along(label), label)
  do.call(rbind, lapply(names(g), function(k) {
    i <- g[[k]]
    data.frame(label = k, n = length(i), gene = panel, mean_lognorm = rowMeans(lognorm[, i, drop = FALSE]),
               detected = rowMeans(cnt[, i, drop = FALSE] > 0))
  }))
}
by_type <- summ(meta$celltypes)
by_sub <- summ(meta$cellsubtypes)
write.csv(by_type, file.path(out_dir, "markers_by_celltype.csv"), row.names = FALSE)
write.csv(by_sub, file.path(out_dir, "markers_by_cellsubtype.csv"), row.names = FALSE)

# per-cell RGC score = mean log expression of pan-RGC genes; non-RGC scores likewise
sc <- sapply(markers[1:9], function(m) colMeans(lognorm[intersect(m, panel), , drop = FALSE]))
sc <- as.data.frame(sc)
sc$celltypes <- meta$celltypes
sc$Sample <- meta$Sample
sc$Condition <- meta$Condition
write.csv(aggregate(. ~ celltypes, sc[, c(names(markers)[1:9], "celltypes")], mean), file.path(out_dir, "panel_score_by_celltype.csv"), row.names = FALSE)

sink(file.path(out_dir, "summary.txt"), split = TRUE)
cat("cells:", nrow(meta), "\n")
cat("percent.mt max:", max(meta$percent.mt), " (authors: <5%)\n")
cat("scDblFinder.score >0.5:", sum(meta$scDblFinder.score > 0.5), "\n")
cat("\ncells with RBPMS detected:", round(mean(cnt["RBPMS", ] > 0), 3), "; any of RBPMS/RBPMS2/SNCG/POU4F1/POU4F2:",
    round(mean(colSums(cnt[intersect(c("RBPMS", "RBPMS2", "SNCG", "POU4F1", "POU4F2"), panel), ]) > 0), 3), "\n")
cat("\nper-celltype mean pan-RGC score vs non-RGC scores (log1p CP10K):\n")
print(round(aggregate(. ~ celltypes, sc[, c(names(markers)[1:9], "celltypes")], mean)[, -1], 2) |> cbind(celltype = sort(unique(sc$celltypes))))
cat("\ncells whose non-RGC class score exceeds the pan-RGC score by class (any class):\n")
nonrgc <- setdiff(names(markers)[1:9], "pan_rgc")
worse <- apply(sc[, nonrgc], 1, max) > sc$pan_rgc
print(table(celltype = sc$celltypes, nonrgc_above_rgc = worse))
cat("\nlabel hierarchy consistent (every cellsubtype maps to one celltype):",
    all(rowSums(table(meta$cellsubtypes, meta$celltypes) > 0) == 1), "\n")
cat("\nsubtype specificity checks:\n")
sp <- function(g, lab) { d <- subset(by_sub, gene == g); d[order(-d$mean_lognorm), ][1:5, c("label", "mean_lognorm", "detected")] }
for (g in c("OPN4", "SPP1", "TBR1", "FOXP2")) { cat("top subtypes for", g, ":\n"); print(sp(g)) }
cat("\nlabel x condition (cells):\n")
print(table(meta$celltypes, meta$Condition))
cat("\nlabel x sample (fraction of sample cells):\n")
print(round(prop.table(table(meta$celltypes, meta$Sample), 2), 3))
sink()
cat("wrote", out_dir, "\n")
