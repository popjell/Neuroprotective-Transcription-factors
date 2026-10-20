# TF target sets x subtype (fgsea) and TF x bulk GO community in MS cortex neurons
source("common.R")
library(fgsea)
library(igraph)

tf_sets <- lapply(network$tf_targets, human)
tf_sets <- lapply(tf_sets, intersect, de$gene)
tf_sets <- tf_sets[lengths(tf_sets) >= 5]
cat("TFs with at least 5 detected targets:", length(tf_sets), "of", length(network$tf_targets), "\n")

# Figure 5: NES of each TF's targets per subtype
enrich <- do.call(rbind, lapply(subtype_order, function(s) {
  d <- de[de$subtype == s & !is.na(de$stat), ]
  res <- fgsea(tf_sets, setNames(d$stat, d$gene), minSize = 5, maxSize = 1000)
  data.frame(subtype = s, TF = res$pathway, NES = res$NES, padj = res$padj, size = res$size)
}))
write.csv(enrich, res("tf_target_nes.csv"), row.names = FALSE)

tf_order <- names(sort(tapply(enrich$NES, enrich$TF, mean, na.rm = TRUE)))
enrich$TF <- factor(enrich$TF, levels = tf_order)
enrich$subtype <- factor(enrich$subtype, levels = subtype_order)
limit <- max(abs(enrich$NES), na.rm = TRUE)
p <- ggplot(enrich, aes(subtype, TF, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(!is.na(padj) & padj < 0.05, "*", "")), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = blue, mid = "white", high = orange, midpoint = 0, limits = c(-limit, limit),
                       na.value = "grey80", name = "NES\n(MS vs control)") +
  labs(title = paste("TF target sets in", tissue), x = NULL, y = NULL,
       caption = paste0("fgsea of each TF's iRegulon targets on the subtype Wald statistic. * BH padj < 0.05 within subtype. TFs need at least 5 detected targets.\n", caption_base)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 9), panel.grid = element_blank()) +
  caption_theme
save_fig(p, "tf_target_nes", 10, 11)

# Figure 6: TF x bulk GO community. Communities are rebuilt as in community_detection/run_community.R (Leiden, seed 42, resolution 1.5)
edges <- network$edges
graph <- as_undirected(graph_from_data_frame(edges[, c("Regulator Gene", "Target Gene")], directed = TRUE), mode = "collapse")
set.seed(42)
membership <- igraph::membership(cluster_leiden(graph, objective_function = "modularity", resolution = 1.5))
labels <- unique(read.csv("../../community_detection/tf_community_enrichment.csv")[, c("Community", "label")])
stopifnot(length(unique(membership)) == nrow(labels))
labels <- labels[grepl(":", labels$label), ]

mean_stat <- tapply(de$stat, de$gene, mean, na.rm = TRUE)
rows <- do.call(rbind, lapply(names(network$tf_targets), function(tf) {
  targets <- network$tf_targets[[tf]]
  do.call(rbind, lapply(labels$Community, function(cl) {
    in_cl <- targets[membership[targets] == cl]
    detected <- intersect(human(in_cl), names(mean_stat))
    if (length(detected) < 3) return(NULL)
    data.frame(TF = tf, community = labels$label[labels$Community == cl],
               pct_of_tf_targets = 100 * length(in_cl) / length(targets), mean_stat = mean(mean_stat[detected]), n_detected = length(detected))
  }))
}))
write.csv(rows, res("tf_community_neurons.csv"), row.names = FALSE)

rows$TF <- factor(rows$TF, levels = names(sort(tapply(rows$mean_stat, rows$TF, mean))))
rows$community <- factor(rows$community, levels = labels$label[order(as.integer(sub("C([0-9]+):.*", "\\1", labels$label)))])
limit <- max(abs(rows$mean_stat))
p <- ggplot(rows, aes(community, TF, size = pct_of_tf_targets, color = mean_stat)) +
  geom_point() +
  scale_color_gradient2(low = blue, mid = "grey85", high = orange, midpoint = 0, limits = c(-limit, limit),
                        name = "Mean Wald stat\n(MS vs control)") +
  scale_size(range = c(1, 7), name = "% of TF targets\nin community") +
  labs(title = paste("TF targets by bulk GO community in", tissue), x = NULL, y = NULL,
       caption = paste0("Communities and labels from the bulk Leiden analysis. Colour: Wald statistic averaged over the TF's detected targets in the community and over subtypes. At least 3 detected targets.\n", caption_base)) +
  scale_x_discrete(labels = function(x) vapply(strwrap(x, 16, simplify = FALSE), paste, "", collapse = "\n")) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  caption_theme
save_fig(p, "tf_community_neurons", 14, 10)
