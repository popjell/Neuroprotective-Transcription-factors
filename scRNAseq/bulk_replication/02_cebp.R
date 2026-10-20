# Cebp family (Cebpa, Cebpb, Cebpd) and their iRegulon targets in MS cortex neurons
source("common.R")

tfs <- c("CEBPA", "CEBPB", "CEBPD")
target_set <- function(pattern) human(unlist(network$tf_targets[grepl(pattern, names(network$tf_targets))]))
sets <- list(Cebp = target_set("^Cebp"), Stat = target_set("^Stat"), Irf = target_set("^Irf"))
cat("targets:", paste(names(sets), lengths(sets), collapse = ", "), "\n")

# Figure 3: the Cebp TFs themselves
tf_effects <- de[de$gene %in% tfs, ]
tf_effects$lower <- tf_effects$log2FoldChange - 1.96 * tf_effects$lfcSE
tf_effects$upper <- tf_effects$log2FoldChange + 1.96 * tf_effects$lfcSE
tf_effects$significant <- !is.na(tf_effects$padj) & tf_effects$padj < 0.05
write.csv(tf_effects, res("cebp_tf_effects.csv"), row.names = FALSE)
p <- ggplot(tf_effects, aes(log2FoldChange, subtype, color = log2FoldChange > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.25) +
  geom_point(aes(fill = ifelse(significant, "filled", "open")), shape = 21, size = 3) +
  scale_color_manual(values = c("TRUE" = orange, "FALSE" = blue), guide = "none") +
  scale_fill_manual(values = c(filled = "black", open = "white"), guide = "none") +
  scale_y_discrete(limits = rev(subtype_order)) +
  facet_wrap(~gene, nrow = 1) +
  labs(title = paste("Cebp transcription factors in", tissue), x = "log2 fold change (MS vs control), 95% interval", y = NULL,
       caption = paste0("Filled points: padj < 0.05. Open points: padj >= 0.05. Missing rows: gene removed by the count filter.\n", caption_base)) +
  caption_theme
save_fig(p, "cebp_tf_effects", 10, 6)

# Figure 4: do the Cebp targets move, compared with Stat and Irf targets
shift <- do.call(rbind, lapply(names(sets), function(n) set_shift(sets[[n]], paste(n, "targets"))))
shift$set <- factor(shift$set, levels = paste(names(sets), "targets"))
write.csv(shift, res("cebp_target_shift.csv"), row.names = FALSE)
p <- plot_shift(shift, paste("Cebp, Stat and Irf target sets in", tissue),
                paste0("Target sets overlap the bulk program by construction, so this tests whether the bulk regulon direction carries over to neurons.\n", caption_base))
save_fig(p, "cebp_target_shift", 12, 6)
