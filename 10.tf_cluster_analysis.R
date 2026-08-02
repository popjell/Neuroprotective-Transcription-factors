# Script 10: TF clustering by Jaccard similarity + Cluster × GO Family dot plot
#
# Performs clustering of TFs based on either:
#   1. Jaccard similarity of target gene sets (hierarchical clustering with dendrogram)
#   2. JASPAR TF family grouping (exact family assignment, no dendrogram)
#
# Clusters are labeled by their dominant TF (Jaccard) or family name (JASPAR).
# Produces a dot plot of TF clusters vs GO parent families.

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(cluster)
library(ggplotify)
library(cowplot)
if (requireNamespace("TFBSTools", quietly = TRUE) && requireNamespace("JASPAR2024", quietly = TRUE)) {
  library(TFBSTools)
  library(JASPAR2024)
}

# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

source("utils/helpers.R")

net <- load_tf_network("reformatted_edges.csv")
all_tfs <- net$all_tfs
tf_targets <- net$tf_targets

expr <- load_expression("Expression_data.txt")
mn_ranked  <- expr$mn_ranked
rgc_ranked <- expr$rgc_ranked

if (!exists("reducedTerms") || is.null(reducedTerms)) {
  stop("'reducedTerms' not found in global environment. Run script 6 first.")
}
if (!exists("go_analysis") || is.null(go_analysis)) {
  stop("'go_analysis' not found in global environment. Run script 6 first.")
}

tf_enrich_file <- "tf_target_enrichment.csv"
if (file.exists(tf_enrich_file)) {
  tf_enrichment <- read.csv(tf_enrich_file, stringsAsFactors = FALSE)
} else {
  tf_enrichment <- NULL
  message("WARNING: tf_target_enrichment.csv not found. Enrichment scores will be computed on-the-fly.")
}

cat("Loaded", length(all_tfs), "TFs and", length(unique(unlist(tf_targets))), "unique targets\n")

# ==============================================================================
# SECTION 2: BUILD GO FAMILY → GENE MAPPING
# ==============================================================================

# Map GO terms to parent categories and extract genes per parent
go_genes_long <- go_analysis %>%
  dplyr::filter(source == "GO:BP") %>%
  dplyr::select(term_id, intersections) %>%
  tidyr::separate_rows(intersections, sep = ",\\s*|\\s+") %>%
  dplyr::rename(gene = intersections) %>%
  dplyr::filter(!is.na(gene), nchar(gene) > 0)

genes_with_parents <- reducedTerms %>%
  dplyr::select(term = go, parentTerm) %>%
  dplyr::inner_join(go_genes_long, by = c("term" = "term_id"))

parent_gene_sets_raw <- genes_with_parents %>%
  dplyr::group_by(parentTerm) %>%
  dplyr::summarise(genes = list(unique(gene)), .groups = "drop") %>%
  tibble::deframe()

# Convert Ensembl IDs to gene symbols if needed
sample_genes <- unlist(parent_gene_sets_raw)[1:5]
is_ensembl <- any(grepl("^ENS", sample_genes))
if (is_ensembl) {
  cat("Converting Ensembl IDs to gene symbols using expression data...\n")
  exp_map <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  id_map <- stats::setNames(exp_map$symbol, exp_map$ensembl_id)
  parent_gene_sets <- lapply(parent_gene_sets_raw, function(genes) {
    symbols <- unique(id_map[genes])
    symbols[!is.na(symbols)]
  })
  parent_gene_sets <- parent_gene_sets[sapply(parent_gene_sets, length) > 0]
  cat("  Mapped", length(parent_gene_sets_raw), "->", length(parent_gene_sets), "families after conversion\n")
} else {
  parent_gene_sets <- parent_gene_sets_raw
}

cat("Found", length(parent_gene_sets), "GO parent families\n")
cat("  Sample parents:", paste(head(names(parent_gene_sets), 3), collapse=", "), "\n")
cat("  Sample gene counts:", paste(sapply(head(parent_gene_sets, 3), length), collapse=", "), "\n")

# ==============================================================================
# SECTION 3: CHOOSE CLUSTERING METHOD AND BUILD DISTANCE MATRIX
# ==============================================================================

if (!exists("skip_prompts") || !skip_prompts) {
  beep(10)
  cat("\nSelect clustering method:\n")
  cat("  1 = Jaccard (target gene overlap)\n")
  cat("  2 = JASPAR (TF family grouping)\n")
  method_input <- readline(prompt = "Enter 1 or 2: ")
  clustering_method <- trimws(method_input)
} else {
  if (!exists("clustering_method")) clustering_method <- "1"
  cat(sprintf("Using pre-set clustering method: %s\n", clustering_method))
}

if (clustering_method == "2") {
  # ---- JASPAR family-based grouping ----
  if (!requireNamespace("TFBSTools", quietly = TRUE) || !requireNamespace("JASPAR2024", quietly = TRUE)) {
    stop("TFBSTools and JASPAR2024 packages are required for JASPAR clustering.\n",
         "Install with: BiocManager::install(c('TFBSTools', 'JASPAR2024'))")
  }

  cat("\nFetching JASPAR2024 mouse CORE motifs...\n")
  jdb <- JASPAR2024::JASPAR2024()
  db_con <- JASPAR2024::db(jdb)
  all_jaspar <- TFBSTools::getMatrixSet(db_con, opts = list(tax_group = "vertebrates", collection = "CORE"))
  jaspar_names <- tolower(sapply(all_jaspar, TFBSTools::name))

  # Build name -> family mapping using tags()
  family_map <- setNames(
    sapply(all_jaspar, function(m) {
      fam <- TFBSTools::tags(m)$family
      if (length(fam) == 0) return("Unknown")
      fam <- paste(fam, collapse = " / ")
      if (is.na(fam) || nchar(fam) == 0) "Unknown" else fam
    }),
    jaspar_names
  )

  # Map our TFs to JASPAR families (exact match first, then partial)
  tf_to_family <- list()
  for (tf in all_tfs) {
    idx <- which(jaspar_names == tolower(tf))
    if (length(idx) == 0) idx <- grep(paste0("^", tolower(tf), "$"), jaspar_names)
    if (length(idx) >= 1) {
      fam <- family_map[jaspar_names[idx[1]]]
      if (!is.na(fam) && nchar(fam) > 0) {
        tf_to_family[[tf]] <- fam
      }
    }
  }

  cat("  Mapped", length(tf_to_family), "/", length(all_tfs), "TFs to JASPAR families\n")
  unmapped <- setdiff(all_tfs, names(tf_to_family))
  if (length(unmapped) > 0) {
    cat("  Unmapped (excluded from clustering):", paste(unmapped, collapse = ", "), "\n")
  }

  if (length(tf_to_family) < 2) stop("Fewer than 2 TFs mapped to JASPAR families. Cannot cluster.")

  # Update all_tfs to only mapped TFs
  all_tfs <- names(tf_to_family)
  n_tfs <- length(all_tfs)

  # Assign cluster IDs by exact family
  families <- sort(unique(unlist(tf_to_family)))
  family_to_cluster <- setNames(seq_along(families), families)
  cl <- setNames(sapply(all_tfs, function(tf) family_to_cluster[tf_to_family[[tf]]]), all_tfs)

  n_clust <- length(families)
  hc <- NULL  # No dendrogram for exact family grouping
  clustering_method_name <- "JASPAR"
  cat(sprintf("  Grouped into %d JASPAR families\n", n_clust))

} else {
  # ---- Jaccard target-based clustering (default) ----
  n_tfs <- length(all_tfs)
  all_genes <- unique(unlist(tf_targets))
  gene_idx  <- setNames(seq_along(all_genes), all_genes)

  n_genes <- length(all_genes)
  bin_mat <- matrix(0L, n_tfs, n_genes, dimnames = list(all_tfs, all_genes))
  for (tf in all_tfs) {
    bin_mat[tf, gene_idx[tf_targets[[tf]]]] <- 1L
  }

  inter_mat <- bin_mat %*% t(bin_mat)
  row_sums  <- rowSums(bin_mat)
  union_mat <- outer(row_sums, row_sums, "+") - inter_mat
  jac_mat   <- inter_mat / union_mat
  diag(jac_mat) <- 1

  dist_mat <- as.dist(1 - jac_mat)
  hc <- hclust(dist_mat, method = "ward.D2")
  clustering_method_name <- "Jaccard"
}

cat("\nClustering method:", clustering_method_name, "\n")
cat("TFs clustered:", n_tfs, "\n")

if (!is.null(hc)) {
  cat("\nSilhouette scores by cluster count:\n")
  for (k in 3:min(15, n_tfs - 1)) {
    si <- mean(silhouette(cutree(hc, k), dist_mat)[, 3], na.rm = TRUE)
    cat("  k =", k, "->", round(si, 3), "\n")
  }
}

# ==============================================================================
# SECTION 4: CLUSTER LABELING
# ==============================================================================

if (clustering_method == "2") {
  # JASPAR family grouping: labels are family names, no cut height needed
  cut_height <- NA
  cluster_ids <- sort(unique(cl))
  n_clust <- length(cluster_ids)

  dominant_tf_per_cluster <- data.frame(
    cluster = cluster_ids,
    dominant_tf = NA_character_,
    label = NA_character_,
    color = rainbow(n_clust),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(dominant_tf_per_cluster))) {
    cs <- dominant_tf_per_cluster$cluster[i]
    members <- names(cl[cl == cs])
    fam_name <- tf_to_family[[members[1]]]

    # Dominant TF = member with the most targets overall
    tf_target_counts <- sapply(members, function(tf) length(tf_targets[[tf]]))
    dom_tf <- names(which.max(tf_target_counts))

    dominant_tf_per_cluster$dominant_tf[i] <- dom_tf
    dominant_tf_per_cluster$label[i] <- fam_name
  }

} else {
  # Jaccard clustering: prompt for cut height and label by dominant TF
  if (!exists("skip_prompts") || !skip_prompts) {
    beep(10)
    user_h <- readline(prompt = "Enter cut height for the dendrogram (e.g. 1.0): ")
    cut_height <- as.numeric(trimws(user_h))
    if (is.na(cut_height) || cut_height <= 0) {
      stop("Invalid cut height. Please enter a positive number.")
    }
  } else {
    if (!exists("cut_height")) cut_height <- 1.0
    cat(sprintf("Using pre-set cut height: %.2f\n", cut_height))
  }

  cl <- cutree(hc, h = cut_height)
  n_clust <- length(unique(cl))
  cluster_ids <- sort(unique(cl))
  cat(sprintf("\nCut height = %.2f → %d clusters\n\n", cut_height, n_clust))

  # Label each cluster with its dominant TF (most targets overall)
  dominant_tf_per_cluster <- data.frame(
    cluster = sort(unique(cl)),
    dominant_tf = NA_character_,
    label = NA_character_,
    color = rainbow(n_clust),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(dominant_tf_per_cluster))) {
    cs <- dominant_tf_per_cluster$cluster[i]
    members <- names(cl[cl == cs])

    # Dominant TF = member with the most targets overall
    tf_target_counts <- sapply(members, function(tf) length(tf_targets[[tf]]))
    dom_tf <- names(which.max(tf_target_counts))

    dominant_tf_per_cluster$dominant_tf[i] <- dom_tf
    dominant_tf_per_cluster$label[i] <- paste0(dom_tf, "_cluster")
  }
}

cat("Cluster assignments:\n")
for (i in seq_len(nrow(dominant_tf_per_cluster))) {
  cs <- dominant_tf_per_cluster$cluster[i]
  members <- names(cl[cl == cs])
  cat(sprintf("  %s (%d TFs): %s\n",
              dominant_tf_per_cluster$label[i], length(members),
              paste(members, collapse = ", ")))
}

# ---- JASPAR family membership table plot ----
if (clustering_method == "2") {
  family_df <- data.frame(
    family = sapply(all_tfs, function(tf) tf_to_family[[tf]]),
    tf = all_tfs,
    stringsAsFactors = FALSE
  )

  family_summary <- family_df %>%
    dplyr::count(family, sort = TRUE) %>%
    dplyr::mutate(
      tfs = sapply(family, function(f) {
        paste(strwrap(paste(sort(family_df$tf[family_df$family == f]), collapse = ", "), width = 28), collapse = "\n")
      }),
      family_wrapped = sapply(family, function(x) paste(strwrap(x, width = 24), collapse = "\n")),
      row = rev(seq_len(dplyr::n()))
    )

  table_df <- data.frame(
    x = rep(1:2, each = nrow(family_summary)),
    y = rep(family_summary$row, 2),
    label = c(family_summary$family_wrapped, family_summary$tfs),
    n = rep(family_summary$n, 2),
    row_bg = rep(ifelse(family_summary$row %% 2 == 0, "grey94", "white"), 2),
    stringsAsFactors = FALSE
  )
  table_df$x <- factor(table_df$x, levels = c(1, 2))

  col_header <- data.frame(
    x = factor(c(1, 2), levels = c(1, 2)),
    y = rep(nrow(family_summary) + 1, 2),
    label = c("JASPAR Family", "Members"),
    stringsAsFactors = FALSE
  )

  p_table <- ggplot(table_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = row_bg), color = "grey80", linewidth = 0.3) +
    geom_text(aes(label = label), hjust = 0, vjust = 0.5, size = 3.3, color = "black", lineheight = 0.9) +
    scale_fill_identity() +
    geom_tile(data = col_header, aes(x = x, y = y), fill = "grey25", color = "grey25", linewidth = 0.3) +
    geom_text(data = col_header, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold", color = "white", lineheight = 0.9) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    labs(title = "TFs Grouped by JASPAR Family") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.margin = margin(10, 15, 10, 15)
    )

  print(p_table)
  ggsave("tf_clustering/jaspar_family_members.png", p_table,
         width = 12, height = max(3, nrow(family_summary) * 0.45 + 1), dpi = 300)
}

# ==============================================================================
# SECTION 4B: IDENTIFY GENES UNIQUE TO EACH CLUSTER
# ==============================================================================

cat("\n--- Identifying genes unique to each cluster ---\n")

cluster_target_sets <- list()
for (cs in cluster_ids) {
  cluster_tfs <- names(cl[cl == cs])
  cluster_target_sets[[cs]] <- unique(unlist(tf_targets[cluster_tfs]))
}

unique_genes_per_cluster <- list()
for (cs in cluster_ids) {
  cluster_genes <- cluster_target_sets[[cs]]
  other_genes <- unique(unlist(cluster_target_sets[setdiff(cluster_ids, cs)]))
  unique_genes <- setdiff(cluster_genes, other_genes)
  unique_genes_per_cluster[[cs]] <- sort(unique_genes)
}

# Print summary
for (cs in cluster_ids) {
  cs_label <- dominant_tf_per_cluster$label[dominant_tf_per_cluster$cluster == cs]
  all_count <- length(cluster_target_sets[[cs]])
  uniq_count <- length(unique_genes_per_cluster[[cs]])
  cat(sprintf("  %s: %d unique / %d total target genes\n", cs_label, uniq_count, all_count))
}

# Save unique genes to CSV
uniq_rows <- list()
for (cs in cluster_ids) {
  cs_label <- dominant_tf_per_cluster$label[dominant_tf_per_cluster$cluster == cs]
  genes <- unique_genes_per_cluster[[cs]]
  if (length(genes) > 0) {
    uniq_rows[[length(uniq_rows) + 1]] <- data.frame(
      cluster = cs,
      cluster_label = cs_label,
      gene = genes,
      stringsAsFactors = FALSE
    )
  }
}
unique_genes_df <- dplyr::bind_rows(uniq_rows)

uniq_csv_name <- if (clustering_method == "2") {
  "tf_clustering/unique_genes_per_cluster_JASPAR.csv"
} else {
  sprintf("tf_clustering/unique_genes_per_cluster_Jaccard_h%.2f.csv", cut_height)
}
if (nrow(unique_genes_df) > 0) {
  write.csv(unique_genes_df, uniq_csv_name, row.names = FALSE)
  cat(sprintf("\nSaved unique genes to %s (%d total unique genes across clusters)\n",
              basename(uniq_csv_name), nrow(unique_genes_df)))
} else {
  cat("\nNo unique genes found across any cluster.\n")
}

# ==============================================================================
# SECTION 5: BUILD CLUSTER × GO FAMILY ENRICHMENT MATRIX
# ==============================================================================

cluster_ids <- sort(unique(cl))

# Compute per-TF per-ParentTerm enrichment if not pre-computed
if (!is.null(tf_enrichment)) {
  tf_go_scores <- tf_enrichment %>%
    dplyr::select(ParentTerm, TF, enrichment_score) %>%
    dplyr::mutate(ParentTerm = gsub("_", " ", ParentTerm))
} else {
  # Compute on-the-fly using fgsea for each TF × parent term
  cat("Computing TF × GO parent enrichment scores on-the-fly...\n")
  tf_go_list <- list()
  for (pt in names(parent_gene_sets)) {
    pt_safe <- gsub("[^[:alnum:]_]", "_", pt)
    pathways <- list()
    for (tf in all_tfs) {
      targs <- intersect(tf_targets[[tf]], parent_gene_sets[[pt]])
      if (length(targs) >= 2) pathways[[tf]] <- targs
    }
    if (length(pathways) == 0) next
    mn_res <- fgsea::fgsea(pathways = pathways, stats = mn_ranked,
                           minSize = 1, maxSize = 500, scoreType = "std")
    rgc_res <- fgsea::fgsea(pathways = pathways, stats = rgc_ranked,
                           minSize = 1, maxSize = 500, scoreType = "std")
    for (tf in names(pathways)) {
      nes_mn <- mn_res$NES[mn_res$pathway == tf]
      nes_rgc <- rgc_res$NES[rgc_res$pathway == tf]
      combined <- mean(c(nes_mn, nes_rgc), na.rm = TRUE)
      if (!is.na(combined)) {
        tf_go_list[[length(tf_go_list) + 1]] <- data.frame(
          ParentTerm = pt, TF = tf, enrichment_score = combined,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tf_go_scores <- dplyr::bind_rows(tf_go_list)
}

# Build the cluster × GO parent matrix
plot_rows <- list()

cat("\nDEBUG: cluster_ids:", length(cluster_ids), "\n")
cat("DEBUG: parent_gene_sets:", length(parent_gene_sets), "\n")
cat("DEBUG: go_scores rows:", nrow(tf_go_scores), "\n")
if (nrow(tf_go_scores) > 0) {
  cat("  Sample TFs in scores:", paste(head(unique(tf_go_scores$TF), 5), collapse=", "), "\n")
  cat("  Sample ParentTerms in scores:", paste(head(unique(tf_go_scores$ParentTerm), 5), collapse=", "), "\n")
}
cat("DEBUG: first cluster TFs:", paste(names(cl[cl == cluster_ids[1]]), collapse=", "), "\n")

for (cs in cluster_ids) {
  cluster_tfs <- names(cl[cl == cs])
  cluster_label <- dominant_tf_per_cluster$label[dominant_tf_per_cluster$cluster == cs]
  cluster_all_targets <- unique(unlist(tf_targets[cluster_tfs]))
  cluster_size <- length(cluster_all_targets)

  for (pt in names(parent_gene_sets)) {
    go_genes <- parent_gene_sets[[pt]]
    overlap <- length(intersect(cluster_all_targets, go_genes))
    pct_coverage <- if (cluster_size > 0) overlap / cluster_size * 100 else 0

    # Average enrichment score of cluster TFs for this GO parent
    cluster_tf_scores <- tf_go_scores %>%
      dplyr::filter(TF %in% cluster_tfs, ParentTerm == pt)

    avg_es <- if (nrow(cluster_tf_scores) > 0) {
      mean(cluster_tf_scores$enrichment_score, na.rm = TRUE)
    } else {
      0
    }

    plot_rows[[length(plot_rows) + 1]] <- data.frame(
      cluster = cs,
      cluster_label = cluster_label,
      ParentTerm = pt,
      pct_coverage = pct_coverage,
      enrichment_score = avg_es,
      stringsAsFactors = FALSE
    )
  }
}

plot_df <- dplyr::bind_rows(plot_rows)
cat("DEBUG: plot_df rows after bind:", nrow(plot_df), "\n")
if (nrow(plot_df) > 0) {
  cat("  pct_coverage range:", range(plot_df$pct_coverage), "\n")
  cat("  enrichment_score range:", range(plot_df$enrichment_score, na.rm = TRUE), "\n")
} else {
  cat("  ERROR: No rows generated. Checking individual cluster/parent combos...\n")
  cs <- cluster_ids[1]
  pt <- names(parent_gene_sets)[1]
  cluster_tfs <- names(cl[cl == cs])
  cluster_all_targets <- unique(unlist(tf_targets[cluster_tfs]))
  go_genes <- parent_gene_sets[[pt]]
  overlap <- length(intersect(cluster_all_targets, go_genes))
  cat("  Test cluster:", cs, "- TFs:", paste(cluster_tfs, collapse=", "), "\n")
  cat("  Test parent:", pt, "- genes:", length(go_genes), "\n")
  cat("  Cluster all targets:", length(cluster_all_targets), "\n")
  cat("  Overlap:", overlap, "\n")
  cat("  TFs in tf_go_scores:", sum(cluster_tfs %in% tf_go_scores$TF), "\n")
}

# Drop GO families with zero coverage across all clusters
active_parents <- plot_df %>%
  dplyr::group_by(ParentTerm) %>%
  dplyr::summarise(max_cov = max(pct_coverage), .groups = "drop") %>%
  dplyr::filter(max_cov > 0) %>%
  dplyr::pull(ParentTerm)

plot_df <- plot_df %>% dplyr::filter(ParentTerm %in% active_parents)
cat("DEBUG: plot_df rows after filtering:", nrow(plot_df), "\n")
if (nrow(plot_df) == 0 && length(active_parents) == 0) {
  cat("  All parents had max coverage of 0\n")
}

# ==============================================================================
# SECTION 6: ORDER AXES
# ==============================================================================
# Change these variables to control how axes are ordered.
# Set to NULL to keep alphabetical/default order.
# Options: "enrichment_score", "pct_coverage", or a manual character vector.
Y_AXIS_ORDER_BY <- "enrichment_score"   # Y-axis: TF clusters
X_AXIS_ORDER_BY <- "pct_coverage"   # X-axis: GO parent families

# --- Y-axis ordering ---
if (!is.null(Y_AXIS_ORDER_BY)) {
  if (is.character(Y_AXIS_ORDER_BY) && length(Y_AXIS_ORDER_BY) > 1) {
    cluster_order <- Y_AXIS_ORDER_BY
  } else if (is.character(Y_AXIS_ORDER_BY) && Y_AXIS_ORDER_BY %in% colnames(plot_df)) {
    cluster_order <- plot_df %>%
      dplyr::group_by(cluster_label) %>%
      dplyr::summarise(order_var = mean(.data[[Y_AXIS_ORDER_BY]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(order_var)) %>%
      dplyr::pull(cluster_label)
  } else {
    cluster_order <- sort(unique(plot_df$cluster_label))
  }
} else {
  cluster_order <- sort(unique(plot_df$cluster_label))
}

# --- X-axis ordering ---
if (!is.null(X_AXIS_ORDER_BY)) {
  if (is.character(X_AXIS_ORDER_BY) && length(X_AXIS_ORDER_BY) > 1) {
    parent_order <- X_AXIS_ORDER_BY
  } else if (is.character(X_AXIS_ORDER_BY) && X_AXIS_ORDER_BY %in% colnames(plot_df)) {
    parent_order <- plot_df %>%
      dplyr::group_by(ParentTerm) %>%
      dplyr::summarise(order_var = mean(.data[[X_AXIS_ORDER_BY]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(order_var) %>%
      dplyr::pull(ParentTerm)
  } else {
    parent_order <- sort(unique(plot_df$ParentTerm))
  }
} else {
  parent_order <- sort(unique(plot_df$ParentTerm))
}

plot_df <- plot_df %>%
  dplyr::mutate(
    cluster_label = factor(cluster_label, levels = cluster_order),
    ParentTerm = factor(ParentTerm, levels = parent_order)
  )

# ==============================================================================
# SECTION 7: DOT PLOT — TF CLUSTERS × GO FAMILIES
# ==============================================================================

cluster_go_dotplot <- ggplot(plot_df, aes(
  x = ParentTerm,
  y = cluster_label,
  size = pct_coverage,
  color = enrichment_score
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", limits = c(0, 10), oob = scales::squish) +
  scale_size_continuous(range = c(2, 10)) +
  theme_bw(base_size = 11) +
  labs(
    title = sprintf("TF Clusters vs GO Parent Families (%s, %d families)", clustering_method_name, n_clust),
    subtitle = "Dot size = % of cluster targets in GO family; Color = enrichment (avg NES)",
    y = "TF Family",
    x = "GO Parent Term",
    size = "% of Cluster\nTargets in\nGO Family",
    color = "Enrichment\nScore (NES)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", face = "bold"),
    panel.grid.minor = element_blank()
  )

print(cluster_go_dotplot)

# ==============================================================================
# SECTION 8: DENDROGRAM WITH DOMINANT-TF CLUSTER LABELS
# ==============================================================================

if (!is.null(hc)) {
  p_dend <- as.ggplot(function() {
    par(mar = c(5, 4, 3, 16))
    plot(hc, main = "", xlab = "TF", ylab = "Distance (Ward D2)", cex = 0.6, hang = -1)
    rect.hclust(hc, h = cut_height, border = dominant_tf_per_cluster$color)
    title(main = sprintf("Dendrogram - %s (h = %.2f, %d clusters)", clustering_method_name, cut_height, n_clust),
          line = 0.5, cex.main = 1.1)
    legend("topright",
           legend = dominant_tf_per_cluster$label,
           fill = dominant_tf_per_cluster$color,
           border = dominant_tf_per_cluster$color,
           bty = "n", cex = 0.75, inset = c(-0.35, 0), xpd = TRUE)
  })

  print(p_dend)
} else {
  p_dend <- NULL
}

# ==============================================================================
# SECTION 9: SAVE OUTPUTS
# ==============================================================================

dir.create("tf_clustering", showWarnings = FALSE)

# PNGs
ggsave(paste0("tf_clustering/cluster_go_dotplot_", tolower(clustering_method_name), ".png"), cluster_go_dotplot,
       width = 14, height = 8, dpi = 300)

# Summary CSV
csv_name <- if (clustering_method == "2") {
  sprintf("tf_clustering/cluster_go_summary_%s.csv", tolower(clustering_method_name))
} else {
  sprintf("tf_clustering/cluster_go_summary_%s_h%.2f.csv", tolower(clustering_method_name), cut_height)
}
write.csv(plot_df, csv_name, row.names = FALSE)

# Dendrogram PNG (Jaccard only)
if (!is.null(p_dend)) {
  ggsave(paste0("tf_clustering/cluster_dendrogram_", tolower(clustering_method_name), ".png"), p_dend,
         width = 14, height = 8, dpi = 300)
}

# PDF (Jaccard: dendrogram + dot plot; JASPAR: dot plot only)
pdf_name <- if (clustering_method == "2") {
  sprintf("tf_clustering/cluster_report_%s.pdf", tolower(clustering_method_name))
} else {
  sprintf("tf_clustering/cluster_report_%s_h%.2f.pdf", tolower(clustering_method_name), cut_height)
}
pdf(pdf_name, width = 14, height = 8)
if (!is.null(p_dend)) print(p_dend)
print(cluster_go_dotplot)
dev.off()

cat(sprintf("\nSaved outputs to tf_clustering/:\n"))
cat(sprintf("  cluster_go_dotplot_%s.png\n", tolower(clustering_method_name)))
if (!is.null(p_dend)) {
  cat(sprintf("  cluster_dendrogram_%s.png\n", tolower(clustering_method_name)))
}
cat(sprintf("  %s\n", basename(csv_name)))
cat(sprintf("  %s\n", basename(pdf_name)))
if (nrow(unique_genes_df) > 0) {
  cat(sprintf("  %s\n", basename(uniq_csv_name)))
}

cat("\n=== Script 10 complete ===\n")
