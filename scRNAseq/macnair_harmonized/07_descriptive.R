# Descriptive outputs: the 31 Rank 2 genes from the primary upper_EN fit (BH within the 31) and per-donor normalized upper_EN PSAP with tissue composition.
# Needs results/primary_upperEN/{deseq2,donors}.csv and the pseudobulk cache from 00-03. Exploratory; not part of the primary claim.
source("common.R")

out_dir <- res("descriptive")
dir.create(out_dir, showWarnings = FALSE)
out_rank2 <- file.path(out_dir, "rank2_upperEN.csv")
out_donor <- file.path(out_dir, "donor_psap.csv")
no_overwrite(out_rank2)
no_overwrite(out_donor)

# 31 Rank 2 genes, BH within the genes that were tested
genes <- rank2_genes()
primary <- read.csv(res("primary_upperEN", "deseq2.csv"))
rank2 <- data.frame(gene = genes, primary[match(genes, primary$gene), c("baseMean", "log2FoldChange", "lfcSE", "pvalue")])
rank2$padj_bh31 <- p.adjust(rank2$pvalue, "BH")
rank2$tested <- !is.na(rank2$pvalue)
write.csv(rank2, out_rank2, row.names = FALSE)

# Per-donor normalized PSAP (primary-fit size factors) and upper_EN cells by tissue
donors <- read.csv(res("primary_upperEN", "donors.csv"))
donors <- donors[donors$eligible, ]
pb <- load_pb()
up <- pb$cols[pb$cols$group == "upper_EN", ]
tissue_cells <- as.data.frame.matrix(xtabs(n_cells ~ donor + tissue, up))
names(tissue_cells) <- paste0("cells_", gsub(" ", "_", names(tissue_cells)))
donor_psap <- cbind(donors[, c("donor", "condition", "n_cells", "psap_count", "size_factor")], tissue_cells[donors$donor, , drop = FALSE])
donor_psap$psap_normalized <- donor_psap$psap_count / donor_psap$size_factor
write.csv(donor_psap, out_donor, row.names = FALSE)
writeLines(capture.output(sessionInfo()), res("sessionInfo.txt"))
cat("Rank 2 genes tested:", sum(rank2$tested), "of 31;", nrow(donor_psap), "donors\n")
