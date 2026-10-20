# Enrichment of each TF's targets in MS vs control neurons, one figure per dataset (run 03_tf_go.R first)
source("common.R")
library(fgsea)

neurons <- subtype_order
d <-de[de$subtype %in% neurons & !is.na(de$stat), ]

# Wald statistic averaged over neuron groups, for genes tested in at least 3 groups
n_groups <- table(d$gene)
pooled <- tapply(d$stat, d$gene, mean)
pooled <- pooled[names(pooled) %in% names(n_groups)[n_groups >= 3]]

tf_sets <- lapply(network$tf_targets, human)
tf_sets <- lapply(tf_sets, intersect, names(pooled))
tf_sets <- tf_sets[lengths(tf_sets) >= 5]
tf <- fgsea(tf_sets, pooled, minSize = 5, maxSize = 1000)
tf <- data.frame(TF = tf$pathway, NES = tf$NES, padj = tf$padj, n_targets = tf$size)
write.csv(tf, res("tf_neuron_enrichment.csv"), row.names = FALSE)

per_group <- read.csv(res("tf_target_nes.csv"))
per_group <- per_group[per_group$subtype %in% neurons & per_group$TF %in% tf$TF, ]

tf_order <- tf$TF[order(tf$NES)]
tf$TF <- factor(tf$TF, levels = tf_order)
per_group$TF <- factor(per_group$TF, levels = tf_order)
tf$significant <- tf$padj < 0.05
p <- ggplot(tf, aes(NES, TF)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(data = per_group, aes(NES), color = "grey70", size = 1.2) +
  geom_point(aes(size = n_targets, fill = ifelse(significant, "filled", "open"), color = NES > 0), shape = 21, stroke = 1) +
  scale_color_manual(values = c("TRUE" = orange, "FALSE" = blue), guide = "none") +
  scale_fill_manual(values = c(filled = "black", open = "white"), guide = "none") +
  scale_size(range = c(2, 6), name = "Detected\ntargets") +
  labs(title = paste("TF target enrichment in", dataset, "neurons"), x = "NES (positive = targets higher in MS)", y = NULL,
       caption = paste0(wrap("Coloured points: fgsea of each TF's targets on the Wald statistic averaged over", length(neurons), "neuron groups (",
                             paste(neurons, collapse = ", "), "). Filled: BH padj < 0.05. Grey points: NES in each single group."), "\n", caption_base)) +
  theme(axis.text.y = element_text(size = 9)) +
  caption_theme
save_fig(p, "tf_neuron_enrichment", 9, 10)
