# Patient-level pseudobulk using the published neuronal subtypes.
library(Matrix)
library(DESeq2)
library(ggplot2)

counts <- readRDS("data/neuron_counts.rds")
meta <- readRDS("data/neuron_metadata.rds")
rownames(meta) <- meta$barcodes

# Keep QC-passing neurons and use the published subtype labels.
meta$subtype <- meta$expertAnno.l2.old
keep <- grepl("^NEU_", meta$subtype)
meta <- meta[keep, ]
counts <- counts[, rownames(meta)]

# Sum raw counts for every patient and neuronal cluster.
group <- paste(meta$patient, meta$subtype, sep = "__")
samples <- sort(unique(group))
membership <- sparseMatrix(i = seq_along(group), j = match(group, samples), x = 1,
                           dims = c(length(group), length(samples)))
pseudobulk <- counts %*% membership
colnames(pseudobulk) <- samples
sample_info <- meta[match(samples, group), c("patient", "disease", "subtype")]
rownames(sample_info) <- samples
sample_info$cells <- tabulate(match(group, samples), length(samples))
saveRDS(list(counts = pseudobulk, metadata = sample_info),
        "results/pseudobulk_counts.rds")
write.csv(sample_info, "results/pseudobulk_samples.csv", row.names = FALSE)

all_results <- list()
all_pca <- list()

for (subtype in sort(unique(sample_info$subtype))) {
  use <- sample_info$subtype == subtype & sample_info$cells >= 20
  info <- sample_info[use, ]

  # A cluster needs at least two control and two MS patients.
  if (sum(info$disease == "CTRL") < 2 || sum(info$disease == "MS") < 2) next

  info$condition <- factor(info$disease, levels = c("CTRL", "MS"))
  y <- round(as.matrix(pseudobulk[, use]))
  y <- y[rowSums(y >= 10) >= 2, ]

  dds <- DESeqDataSetFromMatrix(y, info, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", "MS", "CTRL")))
  result$gene <- rownames(result)
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

# Count significant Rank 1 and Rank 2 genes in each neuronal cluster.
rank_genes <- read.csv("data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
rank_results <- merge(rank_genes, de, by.x = "human_symbol", by.y = "gene")
write.csv(rank_results, "results/rank_gene_results.csv", row.names = FALSE)

rank_results$direction <- ifelse(rank_results$padj < 0.05 & rank_results$log2FoldChange > 0,
                                 "Up", ifelse(rank_results$padj < 0.05,
                                              "Down", "Not significant"))
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
ggsave("plots/rank_gene_barplot.png", p, width = 8, height = 5, dpi = 200)
