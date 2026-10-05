# Cache only cortical neuronal raw counts from the large author Seurat object.
library(Seurat)
library(Matrix)

atlas <- read.delim("data/meta.tsv")
neurons <- unique(atlas$cell_subtypes_intg[
  grepl("-EN$|-IN$|^PC$|^motoEN$", atlas$cell_subtypes_intg)])
atlas <- atlas[atlas$cns_region == "Cortex" &
                 atlas$cell_subtypes_intg %in% neurons, ]

object <- readRDS("data/cor_seurat.rds")
raw <- GetAssayData(object, assay = "RNA", layer = "counts")
cells <- intersect(atlas$cellId, colnames(raw))
if (length(cells) != nrow(atlas))
  stop("Cortical raw counts do not match every integrated neuronal atlas cell")
meta <- atlas[match(cells, atlas$cellId), ]
counts <- raw[, cells, drop = FALSE]
stopifnot(all(colSums(counts) == meta$nCount_RNA))
stopifnot(all(counts@x == round(counts@x)))
rownames(meta) <- cells
rm(object, raw, atlas)
gc()
saveRDS(list(counts = counts, meta = meta), "data/cortex_neurons.rds")
message("Saved raw counts for ", ncol(counts), " cortical neurons")
