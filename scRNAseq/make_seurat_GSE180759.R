suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

gse <- "GSE180759"
data_dir <- file.path(getwd(), "GSE180759")
dir.create(data_dir, showWarnings = FALSE)

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180759/suppl"
files <- c(
  expression = "GSE180759_expression_matrix.csv.gz",
  annotation  = "GSE180759_annotation.txt.gz"
)

download_file <- function(url, dest) {
  if (file.exists(dest)) {
    message("Already downloaded: ", dest)
    return(invisible())
  }
  message("Downloading ", dest, " ...")
  download.file(url, dest, mode = "wb", quiet = TRUE)
}

for (nm in names(files)) {
  download_file(file.path(base_url, files[[nm]]), file.path(data_dir, files[[nm]]))
}

expr_file <- file.path(data_dir, files[["expression"]])
anno_file <- file.path(data_dir, files[["annotation"]])

message("Reading expression matrix ...")
con <- gzfile(expr_file, "rt")
barcodes <- strsplit(readLines(con, n = 1), ",")[[1]]
close(con)

expr <- fread(expr_file, header = FALSE, skip = 1, check.names = FALSE)
genes <- expr[[1]]
mat <- expr[, -1, with = FALSE]
rm(expr)

counts <- t(as(as.matrix(mat), "dgCMatrix"))
rm(mat)
rownames(counts) <- barcodes
colnames(counts) <- genes

message("Reading annotation ...")
anno <- fread(anno_file)
anno <- as.data.frame(anno)
rownames(anno) <- anno$nucleus_barcode

keep <- intersect(rownames(counts), rownames(anno))
message("Cells in expression: ", nrow(counts), "; in annotation: ", nrow(anno),
        "; shared: ", length(keep))
counts <- counts[keep, ]
anno <- anno[keep, , drop = FALSE]

seurat <- CreateSeuratObject(counts = counts, meta.data = anno)

out_rds <- file.path(data_dir, "GSE180759_seurat.rds")
saveRDS(seurat, out_rds)
message("Saved Seurat object (", ncol(seurat), " cells, ", nrow(seurat), " genes) to: ", out_rds)
