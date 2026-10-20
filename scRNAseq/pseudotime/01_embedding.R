# Monocle 3 UMAP of the EN-L2-3 nuclei, aligned across sequencing batches. The author pseudotime is used as is; no trajectory is fitted here.
library(monocle3)
source("common.R")

nuclei <- read.csv(shared("nuclei.csv"))
counts <- readRDS(shared("counts.rds"))

cds <- new_cell_data_set(counts, cell_metadata = data.frame(seq_batch = nuclei$Seqbatch, row.names = nuclei$nucleus),
                         gene_metadata = data.frame(gene_short_name = rownames(counts), row.names = rownames(counts)))
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "seq_batch")
set.seed(1010)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "Aligned")

umap <- reducedDims(cds)[["UMAP"]]
write.csv(data.frame(nucleus = rownames(umap), UMAP1 = umap[, 1], UMAP2 = umap[, 2]), shared("umap.csv"), row.names = FALSE)
cat("UMAP written for", nrow(umap), "nuclei\n")
