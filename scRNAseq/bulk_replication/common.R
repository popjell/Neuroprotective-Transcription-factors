# Shared setup for the bulk-replication scripts. Run them from scRNAseq/bulk_replication/.
# DATASET=Trobisch (default) or DATASET=Schirmer picks the single-nucleus data; outputs go to results/<dataset>/ and figures/<dataset>/.
library(ggplot2)

set.seed(1010)
dataset <- Sys.getenv("DATASET", "Trobisch")
stopifnot(dataset %in% c("Trobisch", "Schirmer"))
res_dir <- file.path("results", dataset)
fig_dir <- file.path("figures", dataset)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res <- function(file) file.path(res_dir, file)

theme_set(theme_bw(base_size = 12))
orange <- "#D9702E"
blue <- "#2F78A8"
caption_theme <- theme(plot.caption = element_blank())  # figures carry no bottom text; details live in the scripts and results csvs
wrap <- function(...) paste(strwrap(paste(...), width = 110), collapse = "\n")
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

source("../../utils/helpers.R")
network <- load_tf_network("../../reformatted_edges.csv")

# One DESeq2 fit per fitted group (Trobisch cortex neuron subtype or Schirmer cell type), MS vs control
if (dataset == "Trobisch") {
  status <- read.csv("../results/fit_status.csv")
  status <- status[status$region == "Cortex" & status$status == "fitted", ]
  subtype_order <- c(sort(status$subtype[grepl("EN$", status$subtype)]), sort(status$subtype[!grepl("EN$", status$subtype)]))
  de_path <- function(s) file.path("../results/Cortex", s, "deseq2.csv")
  tissue <- "MS cortex neurons"
  caption_base <- wrap("Trobisch cortex neurons, MS vs control (4 control / 3 MS donors, 3 / 3 in L5-L6-ADAMTSL1-EN),",
                       "donor-level pseudobulk DESeq2. Gene sets come from the mouse EAE bulk analysis (iRegulon, upper-cased to human symbols).",
                       "Set tests are exploratory; BH is within one gene set across groups.")
} else {
  status <- read.csv("results_schirmer/fit_status.csv")
  status <- status[status$status == "fitted", ]
  neuron_types <- status$cell_type[grepl("^(EN|IN)-", status$cell_type)]
  subtype_order <- c(sort(neuron_types[grepl("^EN", neuron_types)]), sort(neuron_types[grepl("^IN", neuron_types)]))
  de_path <- function(s) file.path("results_schirmer", s, "deseq2.csv")
  tissue <- "MS cortex neurons"
  caption_base <- wrap("Schirmer 2019 cortex snRNA-seq, MS vs control, donor-level pseudobulk DESeq2 (donors with at least 5 cells per cell type,",
                       "at least 3 per condition; microglia have only 3 control donors). Gene sets come from the mouse EAE bulk analysis",
                       "(iRegulon, upper-cased to human symbols). Set tests are exploratory; BH is within one gene set across cell types.")
}
de <- do.call(rbind, lapply(subtype_order, function(s) {
  d <- read.csv(de_path(s))
  d$subtype <- s
  d
}))
de$subtype <- factor(de$subtype, levels = subtype_order)

# Mouse symbols from the bulk analysis to human symbols (upper case; unmatched genes are reported by the caller)
human <- function(x) unique(toupper(x))
mouse_only <- function(x) setdiff(human(x), de$gene)

# Shift in Wald statistic of a gene set vs all other tested genes, per group (positive = higher in MS)
set_shift <- function(set, label) {
  out <- do.call(rbind, lapply(levels(de$subtype), function(s) {
    d <- de[de$subtype == s & !is.na(de$stat), ]
    inset <- d$gene %in% set
    if (sum(inset) < 5) return(NULL)
    tt <- t.test(d$stat[inset], d$stat[!inset])
    data.frame(set = label, subtype = s, n_set = sum(inset), shift = unname(tt$estimate[1] - tt$estimate[2]),
               lower = tt$conf.int[1], upper = tt$conf.int[2],
               p = wilcox.test(d$stat[inset], d$stat[!inset])$p.value)
  }))
  out$padj <- p.adjust(out$p, "BH")
  out
}

plot_shift <- function(shift, title, caption) {
  shift$subtype <- factor(shift$subtype, levels = subtype_order)
  ggplot(shift, aes(shift, subtype, color = shift > 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.25) +
    geom_point(aes(fill = ifelse(padj < 0.05, "filled", "open")), shape = 21, size = 3) +
    scale_color_manual(values = c("TRUE" = orange, "FALSE" = blue), guide = "none") +
    scale_fill_manual(values = c(filled = "black", open = "white"), guide = "none") +
    scale_y_discrete(limits = rev(subtype_order)) +
    facet_wrap(~set, nrow = 1) +
    labs(title = title, x = "Shift in Wald statistic vs other genes, 95% interval (positive = higher in MS)", y = NULL,
         caption = paste0("Filled points: BH padj < 0.05 (Wilcoxon). Open points: padj >= 0.05.\n", caption)) +
    caption_theme
}
