# Cache raw cortical counts for the cells that may send signals to neurons.
library(Seurat)
library(Matrix)

atlas <- read.delim("data/meta.tsv")
sender_types <- c("MG", "GM-AS", "WM-AS", "WM-AS-GFPT2", "BG")
atlas <- atlas[atlas$cns_region == "Cortex" &
                 atlas$cell_subtypes_intg %in% sender_types, ]

object <- readRDS("data/cor_seurat.rds")
raw <- GetAssayData(object, assay = "RNA", layer = "counts")
cells <- intersect(atlas$cellId, colnames(raw))
if (length(cells) != nrow(atlas))
  stop("Cortical sender counts do not match the integrated atlas")
meta <- atlas[match(cells, atlas$cellId), ]
counts <- raw[, cells, drop = FALSE]
stopifnot(all(colSums(counts) == meta$nCount_RNA))
rownames(meta) <- cells
rm(object, raw, atlas)
gc()
saveRDS(list(counts = counts, meta = meta), "data/cortex_senders.rds")
message("Saved raw counts for ", ncol(counts), " cortical sender cells")
