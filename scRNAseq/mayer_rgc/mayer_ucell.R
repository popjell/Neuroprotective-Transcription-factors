# UCell scores per RGC for the 31-gene module (genes detected in >=5% of RGCs), the module without PSAP, and the microglia markers. Writes results/contamination/cell_scores.csv.
source("common.R")
suppressPackageStartupMessages(library(UCell))

min_detect <- 0.05  # fraction of RGCs in which a module gene must be detected
out_dir <- file.path(res_dir, "contamination")
dir.create(out_dir, showWarnings = FALSE)
no_overwrite(file.path(out_dir, "cell_scores.csv"))

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
x <- load_counts()
genes <- attr(x, "genes")[[2]]
stopifnot(identical(attr(x, "cells")[[1]], meta$cell_barcode))
g31 <- read_gene_list(genes)

# module genes: detected in >=5% of all RGCs
detect <- rowMeans(as.matrix(x[match(g31, genes), ]) > 0)
module <- data.frame(gene = g31, detected_fraction = detect, in_module = detect >= min_detect)
cat("module genes (>=", min_detect, " detection):", sum(module$in_module), "of 31\n", sep = "")
dir.create(file.path(res_dir, "celllevel"), showWarnings = FALSE)
write.csv(module, file.path(res_dir, "celllevel", "module_genes.csv"), row.names = FALSE)

# UCell on the full matrix, default maxRank = 1500, one core
rownames(x) <- genes
colnames(x) <- meta$cell_barcode
sig <- list(module31 = module$gene[module$in_module],
            module_no_psap = setdiff(module$gene[module$in_module], primary_gene),
            microglia_ucell = microglia_markers)
sc <- ScoreSignatures_UCell(x, features = sig, maxRank = 1500, chunk.size = 2000, ncores = 1)
sc <- as.data.frame(as.matrix(sc))
names(sc) <- sub("_UCell$", "", names(sc))
stopifnot(identical(rownames(sc), meta$cell_barcode))
cat("UCell done for", nrow(sc), "cells\n")
write.csv(data.frame(cell_barcode = meta$cell_barcode, sample = paste0("S", meta$Sample), condition = meta$Condition,
                     celltypes = meta$celltypes, sc, row.names = NULL), file.path(out_dir, "cell_scores.csv"), row.names = FALSE)
