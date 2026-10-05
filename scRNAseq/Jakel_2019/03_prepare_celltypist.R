# Export the cached sparse matrix for CellTypist.
library(Matrix)

x <- readRDS("data/all_counts.rds")
writeMM(x$counts, "data/celltypist_counts.mtx")
writeLines(rownames(x$counts), "data/celltypist_genes.txt")
writeLines(colnames(x$counts), "data/celltypist_cells.txt")

