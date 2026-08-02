# TF clustering by Jaccard similarity of target genes
# Groups TFs by shared target gene sets, labels each cluster via GO,
# computes enrichment scores, and saves a dendrogram + dot plot.
#
# Usage: Rscript tf_clustering/run_clustering.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(cluster)
library(ggplotify)
library(cowplot)

# ==============================================================================
# CLEAN OLD OUTPUT FILES
# ==============================================================================

old_pdfs <- list.files("tf_clustering", pattern = "^dendrograms_labeled\\.pdf$", full.names = TRUE)
old_pngs <- list.files("tf_clustering", pattern = "^(dendrogram_h|dotplot_h|combined_h).*\\.png$", full.names = TRUE)
old_csvs <- list.files("tf_clustering", pattern = "^tf_cluster_enrichment_h.*\\.csv$", full.names = TRUE)
old_files <- c(old_pdfs, old_pngs, old_csvs)
if (length(old_files) > 0) {
  file.remove(old_files)
  cat("Removed", length(old_files), "old output files\n")
}

# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

root_dir <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    script_arg <- grep("^--file=", args, value = TRUE)
    if (length(script_arg) > 0) {
      normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), ".."))
    } else {
      getwd()
    }
  }
)

edges_path <- file.path(root_dir, "reformatted_edges.csv")
exp_path <- file.path(root_dir, "Expression_data.txt")

if (!file.exists(edges_path)) {
  stop("reformatted_edges.csv not found in project root. Run the main pipeline first.")
}
if (!file.exists(exp_path)) {
  stop("Expression_data.txt not found in project root.")
}

source(file.path(root_dir, "utils", "helpers.R"))

edges <- read.csv(edges_path, stringsAsFactors = FALSE)

expr <- load_expression(exp_path)
exp_map <- expr$exp_map
mn_ranked <- expr$mn_ranked
rgc_ranked <- expr$rgc_ranked

all_tfs <- unique(edges$Regulator.Gene)
tf_targets <- list()
for (tf in all_tfs) tf_targets[[tf]] <- unique(edges$Target.Gene[edges$Regulator.Gene == tf])

cat("Loaded", length(all_tfs), "TFs and", length(unique(edges$Target.Gene)), "unique targets\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

label_clusters <- function(tf_cluster, tf_targets) {
  cls <- sort(unique(tf_cluster))
  n <- length(cls)
  cols <- rainbow(n)
  labels <- data.frame(cluster = cls, label = character(n),
                       color = cols, stringsAsFactors = FALSE)

  for (i in seq_along(cls)) {
    cl <- cls[i]
    cluster_tfs <- names(tf_cluster[tf_cluster == cl])
    cluster_targets <- unique(unlist(tf_targets[cluster_tfs]))

    if (length(cluster_targets) < 5) {
      lbl <- paste0("C", cl, " (", length(cluster_targets), " genes)")
    } else {
      ego <- tryCatch(enrichGO(cluster_targets, org.Mm.eg.db, keyType = "SYMBOL",
                               ont = "BP", pvalueCutoff = 0.05, minGSSize = 5),
                      error = function(e) NULL)
      lbl <- if (!is.null(ego) && nrow(ego) > 0) as.data.frame(ego)$Description[1]
             else paste0("C", cl, " (", length(cluster_targets), " genes)")
    }
    labels$label[i] <- ifelse(grepl("^C[0-9]", lbl), lbl, paste0("C", cl, ": ", lbl))
  }
  labels
}

# ==============================================================================
# SECTION 2: JACCARD SIMILARITY AND CLUSTERING
# ==============================================================================

all_genes <- unique(unlist(tf_targets))
gene_idx <- setNames(seq_along(all_genes), all_genes)

n_tfs <- length(all_tfs)
n_genes <- length(all_genes)
bin_mat <- matrix(0L, n_tfs, n_genes, dimnames = list(all_tfs, all_genes))
for (tf in all_tfs) {
  bin_mat[tf, gene_idx[tf_targets[[tf]]]] <- 1L
}

inter_mat <- bin_mat %*% t(bin_mat)
row_sums <- rowSums(bin_mat)
union_mat <- outer(row_sums, row_sums, "+") - inter_mat
jac_mat <- inter_mat / union_mat
diag(jac_mat) <- 1

dist_mat <- as.dist(1 - jac_mat)
hc <- hclust(dist_mat, method = "ward.D2")

cat("\nSilhouette scores by cluster count:\n")
for (k in 3:15) {
  si <- mean(silhouette(cutree(hc, k), dist_mat)[, 3], na.rm = TRUE)
  cat("  k =", k, "->", round(si, 3), "\n")
}

# ==============================================================================
# SECTION 3-6: PER-HEIGHT ANALYSIS (dendrogram + dot plot + summary)
# ==============================================================================

cut_heights <- c(0.5,1.0,1.5)

for (h in cut_heights) {
  cat(sprintf("\n=== Cut height h = %.1f ===\n", h))

  cl <- cutree(hc, h = h)
  n_clust <- length(unique(cl))
  lab_df <- label_clusters(cl, tf_targets)

  cat(sprintf("  %d clusters\n", n_clust))
  for (i in 1:nrow(lab_df)) cat("   ", lab_df$label[i], "\n")

  # --- Coverage ---
  jac_df <- expand.grid(TF = all_tfs, cluster = sort(unique(cl)), stringsAsFactors = FALSE)
  jac_df$coverage <- 0

  for (cs in sort(unique(cl))) {
    cluster_tfs <- names(cl[cl == cs])
    cluster_genes <- unique(unlist(tf_targets[cluster_tfs]))
    cluster_size <- length(cluster_genes)
    if (cluster_size == 0) next
    overlap <- rowSums(bin_mat[, gene_idx[cluster_genes], drop = FALSE])
    jac_df$coverage[jac_df$cluster == cs] <- overlap[all_tfs] / cluster_size
  }

  # --- Enrichment scores ---
  jac_df$enrichment_score <- NA
  jac_df$is_member <- (cl[jac_df$TF] == jac_df$cluster)

  for (cs in sort(unique(cl))) {
    cluster_tfs <- names(cl[cl == cs])
    cluster_targets <- unique(unlist(tf_targets[cluster_tfs]))
    pathways <- list()
    for (tf in all_tfs) {
      targs <- intersect(tf_targets[[tf]], cluster_targets)
      if (length(targs) >= 2) pathways[[tf]] <- targs
    }
    if (length(pathways) == 0) next
    mn_res <- fgsea(pathways = pathways, stats = mn_ranked,
                     minSize = 1, maxSize = 500, scoreType = "std")
    rgc_res <- fgsea(pathways = pathways, stats = rgc_ranked,
                     minSize = 1, maxSize = 500, scoreType = "std")
    for (tf in names(pathways)) {
      es <- mean(c(mn_res$NES[mn_res$pathway == tf],
                    rgc_res$NES[rgc_res$pathway == tf]), na.rm = TRUE)
      if (!is.na(es)) {
        jac_df$enrichment_score[jac_df$TF == tf & jac_df$cluster == cs] <- es
      }
    }
  }

  jac_df <- jac_df %>% left_join(lab_df, by = "cluster")
  write.csv(jac_df, sprintf("tf_clustering/tf_cluster_enrichment_h%.1f.csv", h),
            row.names = FALSE)

  # --- Dot plot ---
  cl_order <- lab_df$label
  tf_order <- jac_df %>%
    filter(is_member) %>%
    arrange(cluster, desc(enrichment_score)) %>%
    pull(TF) %>% unique()

  jac_df <- jac_df %>%
    mutate(label = factor(label, levels = cl_order),
           TF = factor(TF, levels = tf_order))

  p_dot <- ggplot(jac_df, aes(x = label, y = TF,
                               size = coverage, color = enrichment_score)) +
    geom_point(alpha = 0.85) +
    scale_color_viridis_c(option = "plasma", limits = c(0, 5), oob = scales::squish) +
    scale_size_continuous(range = c(0.5, 8)) +
    theme_bw(base_size = 11) +
    labs(title = sprintf("TF Clusters (h = %.1f, %d clusters)", h, n_clust),
         subtitle = "Dot size = fraction of cluster's target pool regulated by TF",
         x = "TF Cluster (GO label)", y = "Transcription Factor",
         size = "Cluster\ntarget\ncoverage", color = "Enrichment\nScore (NES)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
          axis.text.y = element_text(color = "black", face = "bold", size = 8),
          panel.grid.minor = element_blank())

  # --- Dendrogram with legend ---
  p_dend <- as.ggplot(function() {
    par(mar = c(5, 4, 3, 14))
    plot(hc, main = "", xlab = "TF", ylab = "Distance (Ward D2)", cex = 0.6, hang = -1)
    rect.hclust(hc, h = h, border = lab_df$color)
    title(main = sprintf("Dendrogram (h = %.1f, %d clusters)", h, n_clust),
          line = 0.5, cex.main = 1.1)
    legend("topright", legend = paste0("C", lab_df$cluster, ": ", lab_df$label),
           fill = lab_df$color, border = lab_df$color, bty = "n", cex = 0.75,
           inset = c(-0.35, 0), xpd = TRUE)
  })

  # --- PDF pages: dendrogram then dot plot ---
  print(p_dend)
  print(p_dot)

  # --- Individual PNGs ---
  ggsave(sprintf("tf_clustering/dendrogram_h%.1f.png", h), p_dend,
         width = 14, height = 8, dpi = 300)
  ggsave(sprintf("tf_clustering/dotplot_h%.1f.png", h), p_dot,
         width = 12, height = 10, dpi = 300)
  cat(sprintf("  Saved dendrogram_h%.1f.png + dotplot_h%.1f.png\n", h, h))

  # --- Summary ---
  cat(sprintf("  Cluster summary (h = %.1f):\n", h))
  for (cs in sort(unique(cl))) {
    members <- names(cl[cl == cs])
    lbl <- lab_df$label[lab_df$cluster == cs]
    cat(sprintf("    %s (%d TFs): %s\n", lbl, length(members), paste(members, collapse = ", ")))
  }
}

# --- Multi-page PDF (dendrogram then dot plot per cut height) ---
pdf("tf_clustering/dendrograms_labeled.pdf", width = 14, height = 8)
print(p_dend)
print(p_dot)
dev.off()
cat("\nSaved tf_clustering/dendrograms_labeled.pdf\n")

cat("\n=== DONE ===\n")
