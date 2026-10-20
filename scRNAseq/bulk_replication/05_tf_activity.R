# TF activity in MS vs control with DoRothEA (A-C) and CollecTRI regulons, scored from the per-group DESeq2 Wald statistics (decoupleR ulm)
source("common.R")
library(dorothea)
library(decoupleR)
library(patchwork)

# CollecTRI is downloaded once from the OmniPath web service and cached (OmnipathR 3.18 fails to parse it).
# Sign: inhibition only = -1, otherwise +1; complexes are dropped.
collectri_file <- "results/collectri_human.csv"
if (!file.exists(collectri_file)) {
  options(timeout = 300)
  raw <- read.delim("https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606&license=academic&format=tsv")
  raw <- raw[!grepl("^COMPLEX", raw$source), ]
  net <- data.frame(source = raw$source_genesymbol, target = raw$target_genesymbol,
                    mor = ifelse(as.logical(raw$is_inhibition) & !as.logical(raw$is_stimulation), -1, 1))
  write.csv(net[!duplicated(net[, c("source", "target")]), ], collectri_file, row.names = FALSE)
}
networks <- list(
  dorothea = setNames(dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"), c("tf", "target", "mor")], c("source", "target", "mor")),
  collectri = read.csv(collectri_file)[, c("source", "target", "mor")]
)

bulk_tfs <- toupper(unique(network$edges$`Regulator Gene`))
cebp_tfs <- c("CEBPA", "CEBPB", "CEBPD", "CEBPG")

top_bars <- function(scores, direction, name, network_name) {
  top <- head(scores[order(direction * -scores$score), ], 20)
  top$TF <- factor(top$TF, levels = top$TF[order(top$score)])
  p <- ggplot(top, aes(score, TF)) +
    geom_col(fill = if (direction > 0) orange else blue) +
    geom_vline(xintercept = 0, color = "grey50") +
    labs(title = paste("Top 20", if (direction > 0) "up" else "down", "TFs,", network_name, "in", tissue),
         x = "Mean TF activity score across groups", y = NULL) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p, name, 7, 5)
  p
}

for (net in names(networks)) {
  regulons <- networks[[net]]
  activity <- do.call(rbind, lapply(levels(de$subtype), function(s) {
    d <- de[de$subtype == s & !is.na(de$stat), ]
    mat <- matrix(d$stat, ncol = 1, dimnames = list(d$gene, s))
    out <- as.data.frame(run_ulm(mat, regulons, .source = "source", .target = "target", .mor = "mor", minsize = 5))
    data.frame(subtype = s, TF = out$source, score = out$score, p = out$p_value)
  }))
  activity$padj <- ave(activity$p, activity$subtype, FUN = function(x) p.adjust(x, "BH"))
  write.csv(activity, res(paste0(net, "_activity.csv")), row.names = FALSE)
  activity$subtype <- factor(activity$subtype, levels = subtype_order)
  activity$significant <- activity$padj < 0.05
  label <- c(dorothea = "DoRothEA", collectri = "CollecTRI")[[net]]

  # Bulk iRegulon TFs covered by the network
  bulk <- activity[activity$TF %in% bulk_tfs, ]
  bulk$TF <- factor(bulk$TF, levels = names(sort(tapply(bulk$score, bulk$TF, mean))))
  limit <- max(abs(bulk$score))
  p <- ggplot(bulk, aes(subtype, TF, fill = score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(significant, "*", "")), size = 6, vjust = 0.75) +
    scale_fill_gradient2(low = blue, mid = "white", high = orange, midpoint = 0, limits = c(-limit, limit), name = "TF activity\n(MS vs control)") +
    labs(title = paste(label, "activity of the TFs found in EAE in", tissue), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 9), panel.grid = element_blank())
  save_fig(p, paste0(net, "_bulk_tfs"), 10, 9)
  p_bulk <- p

  # CEBP family
  cebp <- activity[activity$TF %in% cebp_tfs, ]
  p <- ggplot(cebp, aes(score, subtype, color = score > 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = 0, xend = score, yend = subtype), linewidth = 0.6) +
    geom_point(aes(fill = ifelse(significant, "filled", "open")), shape = 21, size = 3) +
    scale_color_manual(values = c("TRUE" = orange, "FALSE" = blue), guide = "none") +
    scale_fill_manual(values = c(filled = "black", open = "white"), guide = "none") +
    scale_y_discrete(limits = rev(subtype_order)) +
    facet_wrap(~TF, nrow = 1) +
    labs(title = paste("CEBP family activity,", label, "in", tissue), x = "TF activity score (filled: BH padj < 0.05)", y = NULL)
  save_fig(p, paste0(net, "_cebp"), 12, 6)
  p_cebp <- p

  # Top 50 up and down by mean activity across groups
  mean_score <- aggregate(score ~ TF, activity, mean)
  p_up <- top_bars(mean_score, 1, paste0(net, "_top20_up"), label)
  p_down <- top_bars(mean_score, -1, paste0(net, "_top20_down"), label)

  # Four-panel summary
  panel <- (p_bulk | p_cebp) / (p_up | p_down) + plot_layout(heights = c(1.6, 1)) + plot_annotation(tag_levels = "A")
  save_fig(panel, paste0(net, "_panel"), 18, 18)
}
