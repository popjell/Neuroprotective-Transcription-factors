# Does the shared MN + RGC upregulated program (215 genes) move in MS cortex cells?
source("common.R")

program <- human(readLines("../../intersected_upregulated_genes.txt"))
write.csv(data.frame(gene = mouse_only(program)), res("program_not_detected.csv"), row.names = FALSE)
detected <- intersect(program, de$gene)
cat("program genes:", length(program), " detected:", length(detected), "\n")

genes <- de[de$gene %in% detected, ]
write.csv(genes, res("program_gene_effects.csv"), row.names = FALSE)

shift <- set_shift(program, "Shared MN + RGC program")
write.csv(shift, res("program_shift.csv"), row.names = FALSE)

# Figure 1: the 40 program genes with the strongest average signal
strongest <- names(sort(tapply(abs(genes$stat), genes$gene, mean, na.rm = TRUE), decreasing = TRUE))[1:40]
top <- genes[genes$gene %in% strongest, ]
gene_order <- names(sort(tapply(top$log2FoldChange, top$gene, mean, na.rm = TRUE)))
top$gene <- factor(top$gene, levels = gene_order)
top$significant <- !is.na(top$padj) & top$padj < 0.05
limit <- max(abs(top$log2FoldChange), na.rm = TRUE)
p <- ggplot(top, aes(subtype, gene, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(significant, "*", "")), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = blue, mid = "white", high = orange, midpoint = 0, limits = c(-limit, limit),
                       na.value = "grey80", name = "log2FC\n(MS vs control)") +
  labs(title = paste("Bulk EAE program genes in", tissue), x = NULL, y = NULL,
       caption = paste0("40 of ", length(detected), " detected program genes with the largest mean |Wald statistic|. * padj < 0.05. Grey: gene removed by the count filter.\n", caption_base)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 9), panel.grid = element_blank()) +
  caption_theme
save_fig(p, "program_gene_effects", 10, 10)

# Figure 2: shift of the whole program
p <- plot_shift(shift, paste("Bulk EAE program in", tissue),
                paste0(length(detected), " of ", length(program), " program genes are detected.\n", caption_base))
save_fig(p, "program_shift", 8, 6)
