# Descriptive outputs: PSAP per Trobisch fine label in Schirmer, and the 31 Rank 2 genes in the upper_EN primary model (BH within the 31).
# Needs results/labels.csv, results/schirmer_neurons.rds and results/primary_upperEN/deseq2.csv from 00-02. Exploratory; not part of the primary claim.
source("common.R")

out_fine <- res("fine_label_psap.csv")
out_rank2 <- res("rank2_upperEN.csv")
no_overwrite(out_fine)
no_overwrite(out_rank2)

s <- readRDS(res("schirmer_neurons.rds"))
lab <- read.csv(res("labels.csv"))
stopifnot(identical(lab$cell, colnames(s$counts)))

# One DESeq2 fit per fine label with at least min_donors control and MS donors holding min_cells cells
fine <- do.call(rbind, lapply(sort(setdiff(unique(lab$label), "unassigned")), function(l) {
  cells <- lab[lab$label == l, ]
  d <- aggregate(cell ~ donor + condition, cells, length)
  e <- d[d$cell >= min_cells, ]
  e$condition <- factor(e$condition, levels = c("Control", "MS"))
  n <- table(e$condition)
  row <- data.frame(label = l, coarse_group = cells$coarse_group[1], n_cells = nrow(cells), control_donors = n[["Control"]], ms_donors = n[["MS"]])
  if (min(n) < min_donors) return(cbind(row, status = "not testable", log2FC = NA, lower = NA, upper = NA, p = NA))
  fit <- fit_deseq(pseudobulk(s$counts[, cells$cell], cells$donor)[, e$donor], e)
  cbind(row, status = "fitted", if (gene %in% fit$res$gene) psap_row(fit$res) else data.frame(log2FC = NA, lower = NA, upper = NA, p = NA))
}))
print(fine)
write.csv(fine, out_fine, row.names = FALSE)

# 31 Rank 2 genes from the primary upper_EN fit; BH within the genes that were tested
rg <- read.csv(rank_genes)
rg <- rg[rg$rank == "Rank 2" & rg$use_primary_analysis == "Y", ]
genes <- sub("^TBC13$", "TBC1D13", rg$human_symbol)
stopifnot(length(genes) == 31)
primary <- read.csv(res("primary_upperEN", "deseq2.csv"))
rank2 <- data.frame(gene = genes, primary[match(genes, primary$gene), c("baseMean", "log2FoldChange", "lfcSE", "pvalue")])
rank2$padj_bh31 <- p.adjust(rank2$pvalue, "BH")
rank2$tested <- !is.na(rank2$pvalue)
write.csv(rank2, out_rank2, row.names = FALSE)
cat("fine labels fitted:", sum(fine$status == "fitted"), "of", nrow(fine), "; Rank 2 genes tested:", sum(rank2$tested), "of 31\n")

writeLines(capture.output(sessionInfo()), res("sessionInfo.txt"))
