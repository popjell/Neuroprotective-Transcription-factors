# Cortex counts (all cell types) + metadata + author UMAP -> data/cortex.rds. Reads ../ms-cross-regional/cor/{seurat.rds,meta.tsv,UMAP.coords.tsv.gz}.
source("common.R")
library(Seurat)

out <- "data/cortex.rds"
no_overwrite(out)
cat("seurat.rds bytes:", file.size(file.path(cortex_dir, "seurat.rds")), "\n")

# Raw integer counts, then drop the Seurat object
counts <- GetAssayData(readRDS(file.path(cortex_dir, "seurat.rds")), assay = "RNA", layer = "counts")
gc()
meta <- read.delim(file.path(cortex_dir, "meta.tsv"), check.names = FALSE)
umap <- read.delim(file.path(cortex_dir, "UMAP.coords.tsv.gz"), header = FALSE, col.names = c("cell", "umap1", "umap2"))
stopifnot(identical(colnames(counts), meta$Cell), all(counts@x == round(counts@x)), all(meta$cell_subtypes %in% names(broad_of)),
          all(meta$sample_id %in% names(donors)), setequal(umap$cell, meta$Cell))

cells <- data.frame(cell = meta$Cell, subtype = meta$cell_subtypes, broad = unname(broad_of[meta$cell_subtypes]), sample_id = meta$sample_id,
                    condition = meta$condition, donor = unname(donors[meta$sample_id]), nUMI = colSums(counts))
cells <- cbind(cells, umap[match(cells$cell, umap$cell), c("umap1", "umap2")])
stopifnot(all(cells$nUMI == meta$nCount_RNA))
saveRDS(list(counts = counts, cells = cells), out)
cat("cells:", nrow(cells), " genes:", nrow(counts), "\n")
print(table(cells$broad, cells$condition))
