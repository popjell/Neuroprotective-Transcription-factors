# Extract raw RNA counts from the authors' cached Seurat objects, one region at a time.
library(Seurat)

metadata <- read.csv("data/new_regions_metadata.csv")
parts <- list()
for (region in c("cer", "sc")) {
  file <- paste0("data/", region, "_seurat.rds")
  if (!file.exists(file)) stop("Download the author Seurat object first: ", file)
  object <- readRDS(file)
  raw <- GetAssayData(object, assay = "RNA", layer = "counts")
  cells <- intersect(metadata$cellId, colnames(raw))
  if (!length(cells)) stop("No matching cells in ", file)
  part_meta <- metadata[match(cells, metadata$cellId), ]
  part_meta$regional_subtype <- as.character(object@meta.data[cells, "cell_subtypes"])
  part_counts <- raw[, cells, drop = FALSE]
  if (!all(colSums(part_counts) == part_meta$nCount_RNA))
    stop("RNA count totals do not match the authors' metadata in ", file)
  parts[[region]] <- list(counts = part_counts, meta = part_meta)
  message(region, ": ", ncol(part_counts), " cells with verified raw RNA counts")
  rm(object, raw, part_counts)
  gc()
}

if (!identical(rownames(parts$cer$counts), rownames(parts$sc$counts)))
  stop("The two RNA assays have different gene lists")
writeLines(rownames(parts$cer$counts), "data/gene_names.txt")
counts <- cbind(parts$cer$counts, parts$sc$counts)
meta <- rbind(parts$cer$meta, parts$sc$meta)
rownames(meta) <- meta$cellId
saveRDS(list(counts = counts, meta = meta), "data/new_regions.rds")
