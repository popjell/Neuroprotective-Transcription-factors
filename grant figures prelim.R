# Grant figures prelim
# Run from the repository root: Rscript "grant figures prelim.R" [optional_output_directory]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(readxl)
  library(GSVA)
})

repo <- getwd()
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else file.path(repo, "grant_figures_prelim_figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ink <- "#202838"
blue <- "#315C86"
coral <- "#C65A5A"
grey <- "#7B8490"
theme_set(theme_classic(base_size = 12) +
  theme(text = element_text(colour = ink), plot.title = element_text(face = "bold"),
        plot.title.position = "plot", strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "top"))
save_fig <- function(plot, name, w = 7.5, h = 5.2, pdf = FALSE) {
  ggsave(file.path(out, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 320, bg = "white")
  if (pdf) ggsave(file.path(out, paste0(name, ".pdf")), plot,
                  width = w, height = h, bg = "white", device = cairo_pdf)
}
pfmt <- function(p) format.pval(p, digits = 2, eps = 0.001)

# 1. Load the three lists and rebuild the shared upregulated gene set.
expr <- read.delim(file.path(repo, "Expression_data.txt"), check.names = FALSE)
rank <- read.csv(file.path(repo, "rank_genes.csv"))
go <- read.csv(file.path(repo, "filtered_genelists/GO_metabolism_genes_expression_data.csv"))
book <- file.path(repo, "Gene_Lists_Rank1_and_Rank2_Pathways.xlsx")
labels <- c(
  nadph = "NADPH supply–demand",
  lysosomal = "PSAP-centered lysosomal lipid/cholesterol",
  metabolism = "GO metabolism"
)
sets <- list(
  nadph = unique(rank$mouse_ensembl[rank$rank == "Rank 1"]),
  lysosomal = unique(rank$mouse_ensembl[rank$rank == "Rank 2"]),
  metabolism = unique(go$ensembl_id)
)
for (i in 1:2) {
  sheet <- read_xlsx(book, sheet = i, col_names = FALSE, .name_repair = "minimal")[[1]]
  book_genes <- sheet[(match("Gene", sheet) + 1):length(sheet)]
  csv_genes <- unique(toupper(rank$mouse_symbol[rank$rank == paste("Rank", i)]))
  stopifnot(setequal(book_genes, csv_genes))
}
stopifnot(!anyDuplicated(expr$ensembl_id))

up <- lapply(c("MN", "RGC"), function(tissue) {
  p <- expr[[paste0("pvalue_", tissue)]]
  fc <- expr[[paste0("log2FoldChange_", tissue)]]
  unique(expr$ensembl_id[!is.na(p) & p < 0.05 & !is.na(fc) & fc > 0])
})
names(up) <- c("MN", "RGC")
shared <- intersect(up$MN, up$RGC)
stopifnot(setequal(shared, scan(file.path(repo, "intersected_upregulated_ensembl_ids.txt"),
                            what = "", quiet = TRUE)))
write.csv(data.frame(tissue = c("Motor neurons", "RGCs", "Shared"),
                     genes = c(length(up$MN), length(up$RGC), length(shared))),
          file.path(out, "venn_counts.csv"), row.names = FALSE)

circle <- function(x) data.frame(x = x + 1.18 * cos(seq(0, 2*pi, length.out = 500)),
                                  y = 1.18 * sin(seq(0, 2*pi, length.out = 500)))
venn <- ggplot() +
  geom_polygon(data = circle(-0.65), aes(x, y), fill = blue, alpha = 0.58,
               colour = blue, linewidth = 0.8) +
  geom_polygon(data = circle(0.65), aes(x, y), fill = coral, alpha = 0.58,
               colour = coral, linewidth = 0.8) +
  annotate("text", x = c(-1.18, 0, 1.18), y = 0,
           label = c(length(setdiff(up$MN, up$RGC)), length(shared),
                     length(setdiff(up$RGC, up$MN))),
           fontface = "bold", size = 7, colour = ink) +
  annotate("text", x = c(-1.18, 1.18), y = 1.48,
           label = c("Motor neurons", "RGCs"), fontface = "bold", size = 4.4) +
  coord_equal(xlim = c(-2.3, 2.3), ylim = c(-1.3, 1.85), expand = FALSE) +
  labs(title = "Genes upregulated in EAE across both neuronal datasets",
       subtitle = "p < 0.05 and log2 fold-change > 0 in each dataset") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = ink),
        plot.subtitle = element_text(colour = grey))
save_fig(venn, "venn_overlap", 8.2, 5.5, pdf = TRUE)

process <- ggplot() +
  annotate("rect", xmin = c(0, 0, 0), xmax = c(3.4, 3.4, 3.4),
           ymin = c(2.1, 1.1, 0.1), ymax = c(2.85, 1.85, 0.85),
           fill = c("#EDF2F7", "#EDF2F7", "#F8ECEB"), colour = NA) +
  annotate("text", x = 0.18, y = c(2.48, 1.48, 0.48), hjust = 0,
           label = c(paste0(format(nrow(expr), big.mark = ","), " genes measured in both datasets"),
                     "EAE versus healthy: p < 0.05, log2FC > 0",
                     paste(length(shared), "genes shared across datasets")),
           size = 4, fontface = c("plain", "plain", "bold"), colour = ink) +
  coord_cartesian(xlim = c(0, 3.5), ylim = c(0, 3), expand = FALSE) +
  theme_void() + labs(title = "From expression data to shared genes") +
  theme(plot.title = element_text(face = "bold", colour = ink))
save_fig(process + venn + plot_layout(widths = c(0.9, 1.3)),
         "venn_with_processing", 13, 5.6, pdf = TRUE)

# 2. Normalize full count matrices, then scale each gene across samples.
read_mn <- function(path) {
  x <- read.delim(gzfile(path), header = FALSE,
                  col.names = c("id", "count"))
  setNames(x$count[!grepl("^__", x$id)], x$id[!grepl("^__", x$id)])
}
mn_files <- list.files(file.path(repo, "GSE104897/extracted_counts"),
                       pattern = "_raw\\.[CD]_[1-5]\\.txt\\.gz$", full.names = TRUE)
stopifnot(length(mn_files) == 10)
mn_raw <- lapply(mn_files, read_mn)
stopifnot(all(vapply(mn_raw, function(x) identical(names(x), names(mn_raw[[1]])), TRUE)))
mn <- do.call(cbind, mn_raw)
rownames(mn) <- names(mn_raw[[1]])
colnames(mn) <- ifelse(grepl("_raw\\.C_", mn_files),
                        paste0("Healthy ", sub(".*\\.C_([1-5]).*", "\\1", mn_files)),
                        paste0("EAE ", sub(".*\\.D_([1-5]).*", "\\1", mn_files)))
mn_group <- rep(c("Healthy", "EAE"), each = 5)

rgc <- read_xlsx(file.path(repo, "RGC analysis/GSE247173/GSE247173_RNA_Metadata_EAE_RGC.xlsx"),
                 .name_repair = "minimal")
rgc_samples <- c("Sample_naive_n1", "Sample_naive_n2", "Sample_naive_n4", "Sample_naive_n5",
                 "Sample_onset_n1", "Sample_onset_n2", "Sample_onset_n3", "Sample_onset_n4",
                 "Sample_peak_n2")
rgc <- rgc[!grepl("^__", rgc[[1]]), c(1, match(rgc_samples, names(rgc)))]
rgc_mat <- as.matrix(rgc[, -1])
storage.mode(rgc_mat) <- "numeric"
rownames(rgc_mat) <- rgc[[1]]
stopifnot(!anyDuplicated(rownames(rgc_mat)), !anyNA(rgc_mat))
colnames(rgc_mat) <- sub("Sample_", "", rgc_samples)
rgc_group <- c(rep("Healthy", 4), rep("EAE", 5))

logcpm <- function(mat) log2(sweep(mat, 2, colSums(mat), "/") * 1e6 + 1)
matrices <- list(MN = logcpm(mn), RGC = logcpm(rgc_mat))
groups <- list(MN = mn_group, RGC = rgc_group)
symbol_of <- setNames(expr$symbol, expr$ensembl_id)
symbol_of[rank$mouse_ensembl] <- rank$mouse_symbol

rows_for <- function(ids, tissue) {
  if (tissue == "MN") return(intersect(ids, rownames(matrices$MN)))
  symbols <- unique(unname(symbol_of[ids]))
  intersect(symbols[!is.na(symbols) & nzchar(symbols)], rownames(matrices$RGC))
}
zscale <- function(x) {
  z <- t(scale(t(x)))
  z[apply(z, 1, function(row) all(is.finite(row))), , drop = FALSE]
}
coverage <- list()
for (key in names(sets)) for (tissue in names(matrices)) {
  rows <- rows_for(sets[[key]], tissue)
  mat <- matrices[[tissue]][rows, , drop = FALSE]
  z <- zscale(mat)
  if (tissue == "MN") rownames(z) <- make.unique(ifelse(
    is.na(symbol_of[rownames(z)]) | !nzchar(symbol_of[rownames(z)]),
    rownames(z), symbol_of[rownames(z)]))
  ann <- data.frame(Condition = factor(groups[[tissue]], levels = c("Healthy", "EAE")),
                    row.names = colnames(z))
  pheatmap(z, cluster_cols = FALSE, cluster_rows = TRUE,
           annotation_col = ann,
           annotation_colors = list(Condition = c(Healthy = blue, EAE = coral)),
           color = colorRampPalette(c("#315C86", "white", "#B75657"))(101),
           breaks = seq(-2.5, 2.5, length.out = 102), border_color = NA,
           fontsize_row = if (nrow(z) > 60) 5.8 else 8,
           fontsize_col = 9, angle_col = 45,
           main = paste0(labels[[key]], " · ", ifelse(tissue == "MN", "Motor neurons", "RGCs"),
                         " (", nrow(z), " genes)"),
           filename = file.path(out, paste0("heatmap_", key, "_", tolower(tissue), ".png")),
           width = 9.5, height = max(5, 2.8 + 0.23 * nrow(z)))
  coverage[[length(coverage) + 1]] <- data.frame(list = labels[[key]], tissue = tissue,
    input_genes = length(sets[[key]]), matched = length(rows), plotted = nrow(z))
}
write.csv(do.call(rbind, coverage), file.path(out, "heatmap_coverage.csv"), row.names = FALSE)

# 3. Compare each list with the other upregulated genes in the same dataset.
plot_long <- list()
sample_long <- list()
gsva_long <- list()
stats <- list()
for (tissue in names(matrices)) {
  fc <- paste0("log2FoldChange_", tissue)
  bg <- up[[tissue]]
  for (key in names(sets)) {
    targets <- intersect(sets[[key]], bg)
    others <- setdiff(bg, sets[[key]])
    g <- expr[match(c(targets, others), expr$ensembl_id), ]
    effects <- data.frame(list = labels[[key]], tissue = tissue,
      category = rep(c("Gene list", "Other upregulated"), c(length(targets), length(others))),
      log2FC = g[[fc]])
    plot_long[[length(plot_long) + 1]] <- effects

    target_rows <- rows_for(sets[[key]], tissue)
    bg_rows <- setdiff(rows_for(bg, tissue), target_rows)
    z <- zscale(matrices[[tissue]][unique(c(target_rows, bg_rows)), , drop = FALSE])
    target_rows <- intersect(target_rows, rownames(z))
    bg_rows <- intersect(bg_rows, rownames(z))
    sample_diff <- colMeans(z[target_rows, , drop = FALSE]) -
                   colMeans(z[bg_rows, , drop = FALSE])
    sample_long[[length(sample_long) + 1]] <- data.frame(
      list = labels[[key]], tissue = tissue, sample = names(sample_diff),
      condition = factor(groups[[tissue]], levels = c("Healthy", "EAE")),
      difference = as.numeric(sample_diff))

    gene_sets <- list(Target = target_rows, Background = bg_rows)
    gs <- GSVA::gsva(GSVA::gsvaParam(matrices[[tissue]], gene_sets,
                                     kcdf = "Gaussian", minSize = 5, verbose = FALSE),
                     verbose = FALSE)
    gs_diff <- gs["Target", ] - gs["Background", ]
    gsva_long[[length(gsva_long) + 1]] <- data.frame(
      list = labels[[key]], tissue = tissue, sample = names(gs_diff),
      condition = factor(groups[[tissue]], levels = c("Healthy", "EAE")),
      difference = as.numeric(gs_diff))

    stats[[length(stats) + 1]] <- data.frame(list = labels[[key]], tissue = tissue,
      n_list = length(targets), n_background = length(others),
      wilcoxon_p = wilcox.test(effects$log2FC[effects$category == "Gene list"],
                              effects$log2FC[effects$category == "Other upregulated"],
                              alternative = "greater")$p.value,
      ks_p = suppressWarnings(ks.test(effects$log2FC[effects$category == "Gene list"],
                                      effects$log2FC[effects$category == "Other upregulated"])$p.value),
      sample_p = t.test(sample_diff ~ groups[[tissue]])$p.value,
      gsva_p = t.test(gs_diff ~ groups[[tissue]])$p.value)
  }
}
effects <- do.call(rbind, plot_long)
samples <- do.call(rbind, sample_long)
gsva_scores <- do.call(rbind, gsva_long)
tests <- do.call(rbind, stats)
write.csv(tests, file.path(out, "comparison_tests.csv"), row.names = FALSE)

palette <- c("Gene list" = coral, "Other upregulated" = grey)
make_plot <- function(key, tissue, method) {
  label <- labels[[key]]
  st <- tests[tests$list == label & tests$tissue == tissue, ]
  if (method %in% c("foldchange", "ecdf")) {
    d <- effects[effects$list == label & effects$tissue == tissue, ]
    if (method == "foldchange") {
      p <- ggplot(d, aes(category, log2FC, fill = category)) +
        geom_violin(width = 0.8, alpha = 0.25, colour = NA, trim = TRUE) +
        geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.5) +
        scale_fill_manual(values = palette) +
        labs(y = "EAE vs healthy log2 fold-change", x = NULL,
             subtitle = paste0("Gene list n = ", st$n_list, "; other upregulated n = ",
                               st$n_background, "; p = ", pfmt(st$wilcoxon_p))) +
        theme(legend.position = "none")
    } else {
      p <- ggplot(d, aes(log2FC, colour = category)) +
        stat_ecdf(linewidth = 1) + scale_colour_manual(values = palette) +
        labs(x = "EAE vs healthy log2 fold-change", y = "Cumulative fraction of genes",
             colour = NULL, subtitle = paste0("KS p = ", pfmt(st$ks_p)))
    }
  } else {
    d <- if (method == "sample") samples else gsva_scores
    d <- d[d$list == label & d$tissue == tissue, ]
    p <- ggplot(d, aes(condition, difference, colour = condition)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
      geom_boxplot(width = 0.45, outlier.shape = NA, linewidth = 0.6) +
      geom_jitter(width = 0.07, size = 2.8) +
      scale_colour_manual(values = c(Healthy = blue, EAE = coral)) +
      labs(x = NULL,
        y = if (method == "sample") "Mean gene z-score: list − background"
            else "GSVA score: list − background",
        subtitle = paste0("EAE vs healthy p = ",
                          pfmt(if (method == "sample") st$sample_p else st$gsva_p))) +
      theme(legend.position = "none")
  }
  p + labs(title = label, caption = if (tissue == "MN") "Motor neurons" else "RGCs") +
    theme(plot.caption = element_text(hjust = 0, colour = grey))
}

for (method in c("foldchange", "ecdf", "sample", "gsva")) {
  for (tissue in names(matrices)) {
    panels <- lapply(names(sets), function(key) {
      p <- make_plot(key, tissue, method)
      save_fig(p, paste(method, key, tolower(tissue), sep = "_"))
      p
    })
    title <- switch(method,
      foldchange = "EAE fold-change among upregulated genes",
      ecdf = "Distribution of EAE fold-change",
      sample = "Gene-list expression relative to other upregulated genes",
      gsva = "GSVA enrichment relative to other upregulated genes")
    summary <- (panels[[1]] | panels[[2]] | panels[[3]]) +
      plot_annotation(title = paste(title, "·", if (tissue == "MN") "Motor neurons" else "RGCs"),
                      theme = theme(plot.title = element_text(face = "bold", size = 17)))
    save_fig(summary, paste("summary", method, tolower(tissue), sep = "_"),
             w = 16.5, h = 5.5, pdf = TRUE)
  }
}

cat("Figures saved to:", normalizePath(out), "\n")
print(do.call(rbind, coverage))
