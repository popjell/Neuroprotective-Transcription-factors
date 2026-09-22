# Patient-level pseudobulk differential expression with DESeq2.
library(Matrix)
library(DESeq2)
library(ggplot2)

x <- readRDS("data/neuron_counts.rds")
counts <- x$counts
meta <- x$meta
gene_map <- x$genes

# Use the author labels, combining the two layer 2/3 labels.
meta$subtype <- as.character(meta$cell_type)
meta$subtype[grepl("^EN-L2-3", meta$subtype)] <- "EN-L2-3"
meta$subtype[meta$subtype == "IN-PV"] <- "IN-PVALB"

keep <- meta$genes >= 500 & meta$UMIs >= 1000
meta <- meta[keep, ]
counts <- counts[, keep]

# Sum raw counts for every patient and neuronal subtype.
group <- paste(meta$patient, meta$subtype, sep = "__")
samples <- sort(unique(group))
membership <- sparseMatrix(i = seq_along(group), j = match(group, samples), x = 1,
                           dims = c(length(group), length(samples)))
pseudobulk <- counts %*% membership
colnames(pseudobulk) <- samples
sample_info <- meta[match(samples, group), c("patient", "diagnosis", "subtype")]
rownames(sample_info) <- samples
sample_info$cells <- tabulate(match(group, samples), length(samples))
saveRDS(list(counts = pseudobulk, metadata = sample_info),
        "results/pseudobulk_counts.rds")

all_results <- list()
all_pca <- list()

for (subtype in unique(sample_info$subtype)) {
  use <- sample_info$subtype == subtype & sample_info$cells >= 20
  info <- sample_info[use, ]
  info$condition <- factor(info$diagnosis, levels = c("Control", "MS"))
  y <- round(as.matrix(pseudobulk[, use]))
  y <- y[rowSums(y >= 10) >= 2, ]

  dds <- DESeqDataSetFromMatrix(y, info, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
  result$gene_id <- rownames(result)
  result$gene <- gene_map$gene_symbol[match(result$gene_id, gene_map$ensembl_gene_id)]
  result$subtype <- subtype
  all_results[[subtype]] <- result

  vst_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  variable_genes <- order(apply(vst_counts, 1, var), decreasing = TRUE)[1:min(500, nrow(vst_counts))]
  pca <- prcomp(t(vst_counts[variable_genes, ]))
  all_pca[[subtype]] <- data.frame(info, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
}

de <- do.call(rbind, all_results)
pca <- do.call(rbind, all_pca)
write.csv(de, "results/deseq2_results.csv", row.names = FALSE)

p <- ggplot(pca, aes(PC1, PC2, color = condition, label = patient)) +
  geom_point(size = 2) + geom_text(size = 2, vjust = -0.6, show.legend = FALSE) +
  facet_wrap(~subtype, scales = "free") + theme_minimal()
ggsave("plots/pseudobulk_pca.png", p, width = 12, height = 9, dpi = 200)

p <- ggplot(de, aes(log2FoldChange, -log10(pvalue), color = padj < 0.05)) +
  geom_point(size = 0.5, alpha = 0.6) + facet_wrap(~subtype) +
  scale_color_manual(values = c("grey70", "#D73027"), na.value = "grey70") +
  theme_minimal() + labs(color = "FDR < 0.05")
ggsave("plots/volcano.png", p, width = 12, height = 9, dpi = 200)

# Count significant Rank 1 and Rank 2 genes in each subtype.
rank_genes <- read.csv("data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_ensembl", "human_symbol")])
rank_results <- merge(rank_genes, de, by.x = "human_ensembl", by.y = "gene_id")
write.csv(rank_results, "results/rank_gene_results.csv", row.names = FALSE)

rank_results$direction <- ifelse(rank_results$padj < 0.05 & rank_results$log2FoldChange > 0,
                                 "Up", ifelse(rank_results$padj < 0.05, "Down", "Not significant"))
significant_rank <- rank_results[!is.na(rank_results$direction) &
                                   rank_results$direction != "Not significant", ]
bars <- expand.grid(subtype = unique(de$subtype), rank = c("Rank 1", "Rank 2"),
                    direction = c("Up", "Down"))
if (nrow(significant_rank)) {
  observed <- aggregate(human_symbol ~ subtype + rank + direction, significant_rank, length)
  bars <- merge(bars, observed, all.x = TRUE)
  names(bars)[4] <- "genes"
} else {
  bars$genes <- 0
}
bars$genes[is.na(bars$genes)] <- 0
write.csv(bars, "results/rank_gene_counts.csv", row.names = FALSE)

p <- ggplot(bars, aes(subtype, genes, fill = direction)) +
  geom_col(position = "dodge") + facet_wrap(~rank) + coord_flip() +
  scale_y_continuous(limits = c(0, max(1, bars$genes))) +
  theme_minimal() + labs(x = NULL, y = "Significant genes", fill = NULL)
ggsave("plots/rank_gene_barplot.png", p, width = 9, height = 6, dpi = 200)
