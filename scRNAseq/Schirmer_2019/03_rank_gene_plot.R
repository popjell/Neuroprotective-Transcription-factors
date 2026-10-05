# Plot significant Rank 1 and Rank 2 genes.
library(ggplot2)

r <- read.csv("results/rank_gene_results.csv")
r <- r[!is.na(r$padj) & r$padj < 0.05, ]
r$label <- paste(r$human_symbol, r$subtype, sep = " / ")

p <- ggplot(r, aes(reorder(label, log2FoldChange), log2FoldChange,
                   fill = log2FoldChange > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("#3B82F6", "#D73027"),
                    labels = c("Down", "Up"), name = NULL) +
  theme_minimal() +
  labs(title = "Schirmer: significant Rank 1/2 genes",
       x = NULL, y = "log2 fold change")
ggsave("plots/rank_gene_effects.png", p, width = 9, height = 5, dpi = 200)

# Show how many genes from each rank list could be tested and passed FDR.
genes <- read.csv("data/rank_genes.csv")
genes <- genes[genes$use_primary_analysis == "Y", ]
results <- read.csv("results/rank_gene_results.csv")

summary <- data.frame()
for (rank in c("Rank 1", "Rank 2")) {
  listed <- unique(genes$human_ensembl[genes$rank == rank])
  tested <- unique(results$human_ensembl[results$rank == rank & !is.na(results$padj)])
  significant <- unique(results$human_ensembl[results$rank == rank &
                                                  !is.na(results$padj) & results$padj < 0.05])
  summary <- rbind(summary, data.frame(
    rank = rank,
    status = c("FDR < 0.05", "Tested, not significant", "No adjusted p-value"),
    genes = c(length(significant), length(setdiff(tested, significant)),
              length(setdiff(listed, tested)))
  ))
}
summary$status <- factor(summary$status,
                         levels = c("No adjusted p-value", "Tested, not significant", "FDR < 0.05"))

p <- ggplot(summary, aes(rank, genes, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(genes > 0, genes, "")),
            position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("grey80", "#78A6C8", "#D73027")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 30, 10)) +
  theme_minimal(base_size = 13) +
  labs(title = "Schirmer: overlap with Rank 1 and Rank 2 genes",
       subtitle = "A gene counts as significant if it passes FDR in at least one neuronal subtype",
       x = NULL, y = "Number of genes", fill = NULL)
ggsave("plots/rank_gene_support.png", p, width = 9, height = 5, dpi = 200)
