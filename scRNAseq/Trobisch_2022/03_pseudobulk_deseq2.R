# Compare MS and control within each neuronal population and region.
library(Matrix)
library(DESeq2)
library(ggplot2)

x <- readRDS("data/new_regions.rds")
counts <- x$counts
meta <- x$meta
if (!all(colSums(counts) == meta$nCount_RNA)) {
  stop("Pseudobulk requires raw RNA counts; this matrix is not the raw RNA assay.")
}

# Keep nuclei with a neuronal label.
is_neuron <- grepl("(^L[2-6]|-IN$|-EN$|^PC$|-NEU$)", meta$regional_subtype)
meta <- meta[is_neuron, ]
counts <- counts[, rownames(meta)]

# Sum raw counts for each donor, region, and neuronal subtype.
meta$subtype <- meta$regional_subtype
group <- paste(meta$cns_region, meta$subtype, meta$sample_id, sep = "__")
samples <- sort(unique(group))
membership <- sparseMatrix(i = seq_along(group), j = match(group, samples),
                           x = 1, dims = c(length(group), length(samples)))
pseudobulk <- counts %*% membership
colnames(pseudobulk) <- samples
  info <- meta[match(samples, group), c("sample_id", "condition", "cns_region", "subtype", "sex")]
rownames(info) <- samples
info$cells <- tabulate(match(group, samples), length(samples))
write.csv(info, "results/pseudobulk_samples.csv", row.names = FALSE)
saveRDS(list(counts = pseudobulk, metadata = info), "results/pseudobulk_counts.rds")

results_list <- list()
pca_list <- list()
eligible <- list()
for (region in unique(info$cns_region)) {
  for (subtype in unique(info$subtype[info$cns_region == region])) {
    use <- info$cns_region == region & info$subtype == subtype
    d <- info[use, ]
    eligible[[paste(region, subtype)]] <- data.frame(
      region = region, subtype = subtype,
      control_donors = sum(d$condition == "Control"),
      ms_donors = sum(d$condition == "MS"))
    if (sum(d$condition == "Control") < 2 || sum(d$condition == "MS") < 2) next

    # Test every donor with cells in this subtype using raw pseudobulk counts.
    d$condition <- factor(d$condition, levels = c("Control", "MS"))
    d$sex <- factor(d$sex)
    y <- as.matrix(pseudobulk[, use])
    y <- y[rowSums(y >= 10) >= 2, ]
    dds <- DESeqDataSetFromMatrix(y, d, design = ~sex + condition)
    dds <- DESeq(dds, quiet = TRUE)
    result <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
    result$gene <- rownames(result)
    result$region <- region
    result$subtype <- subtype
    results_list[[paste(region, subtype)]] <- result

    transformed <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
    variable <- order(apply(transformed, 1, var), decreasing = TRUE)
    pca <- prcomp(t(transformed[head(variable, 500), ]))
    pca_list[[paste(region, subtype)]] <- data.frame(d, PC1 = pca$x[, 1],
                                                    PC2 = pca$x[, 2])
  }
}
write.csv(do.call(rbind, eligible), "results/test_eligibility.csv", row.names = FALSE)

de <- do.call(rbind, results_list)
de$global_padj <- p.adjust(de$pvalue, method = "BH")
write.csv(de, "results/deseq2_results.csv", row.names = FALSE)
write.csv(de[!is.na(de$padj) & de$padj < 0.05, ],
          "results/significant_genes.csv", row.names = FALSE)

# Count significant genes in each test.
bars <- aggregate(sig ~ subtype,
                  data.frame(subtype = de$subtype,
                             sig = as.integer(!is.na(de$padj) & de$padj < 0.05)), sum)
names(bars)[2] <- "genes"
p <- ggplot(bars, aes(subtype, genes)) +
  geom_col(fill = "#3B82F6", width = 0.6) +
  geom_text(aes(label = genes), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() + labs(title = "Trobisch: MS vs control neuronal DE",
                         subtitle = "Donor-level pseudobulk; FDR < 0.05 within subtype",
                         x = NULL, y = "Significant genes")
ggsave("plots/de_gene_counts.png", p, width = 7, height = 4, dpi = 200)

# Show donor samples and gene-level results.
pca <- do.call(rbind, pca_list)
p <- ggplot(pca, aes(PC1, PC2, color = condition, label = sample_id)) +
  geom_point(size = 3) + geom_text(vjust = -0.7, size = 3, show.legend = FALSE) +
  facet_wrap(~subtype, scales = "free") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) + theme_minimal()
ggsave("plots/pseudobulk_pca.png", p, width = 9, height = 5, dpi = 200)

de$significant <- !is.na(de$padj) & de$padj < 0.05
p <- ggplot(de[!is.na(de$pvalue), ],
            aes(log2FoldChange, -log10(pvalue), color = significant)) +
  geom_point(size = 0.5, alpha = 0.6) + facet_wrap(~subtype) +
  scale_color_manual(values = c("grey70", "#D73027"), na.value = "grey70") +
  theme_minimal() + labs(color = "FDR < 0.05")
ggsave("plots/volcano.png", p, width = 9, height = 5, dpi = 200)

# Check the Rank 1 and Rank 2 gene lists against the DE table.
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
rank_results <- merge(rank_genes, de, by.x = "human_symbol", by.y = "gene")
write.csv(rank_results, "results/rank_gene_results.csv", row.names = FALSE)
p <- ggplot(rank_results[!is.na(rank_results$pvalue), ],
            aes(log2FoldChange, -log10(pvalue),
                              color = significant)) +
  geom_point(size = 1.5) + facet_wrap(~subtype) +
  scale_color_manual(values = c("grey70", "#D73027"), na.value = "grey70") +
  theme_minimal() + labs(title = "Rank 1/2 genes", color = "FDR < 0.05")
ggsave("plots/rank_genes.png", p, width = 9, height = 5, dpi = 200)
