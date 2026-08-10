# Script 14: TF cluster/family × gene-list bubble plots
#
# Iterates over the 4 filtered genelists (GO_cell_death, GO_metabolism,
# LLM_cell_death, LLM_metabolism) exactly like script 13, reusing the
# <name>_cytoscape_output objects script 13 made plus the cluster assignments
# from script 10.
#
# For each (TF cluster × gene list) pair it collects the target genes of the
# cluster TFs that appear in that gene list, then plots a bubble chart with:
#   x    = gene list
#   y    = TF cluster / family
#   size = % out-degree (targets in the gene list / total cluster targets)
#   color= NES (fgsea enrichment of those targets, averaged MN + RGC)
#   label= number of targets in the bubble
#
# A per-TF version of the same plot is also produced.
#
# Outputs:
#   filtered_genelists/bubble_plots/tf_families_by_genelist.png/.csv
#   filtered_genelists/bubble_plots/tfs_by_genelist.png/.csv
#   (CSVs include the full target-gene lists per bubble)

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)

source("utils/helpers.R")

# ==============================================================================
# SECTION 1: OBJECTS FROM EARLIER SCRIPTS
# ==============================================================================

if (!exists("cl")) {
  stop("'cl' (cluster assignments) not found. Run script 10 first.")
}
if (!exists("clustering_method")) {
  clustering_method <- "1"
}

net <- load_tf_network("reformatted_edges.csv")
tf_targets <- net$tf_targets

expr <- load_expression("Expression_data.txt")
mn_ranked  <- expr$mn_ranked
rgc_ranked <- expr$rgc_ranked

out_dir <- "filtered_genelists/bubble_plots"
dir.create(out_dir, showWarnings = FALSE)

gene_list_order <- sub("\\_graph_output.csv$", "",
                       list.files("filtered_genelists/graph_output", pattern = "\\.csv$"))
gene_list_order <- sub("\\_graph_output.csv$", "", gene_list_order)

# ==============================================================================
# SECTION 2: CLUSTER LABELS (same conventions as script 10)
# ==============================================================================

cluster_label <- function(cid) {
  if (exists("dominant_tf_per_cluster") && cid %in% dominant_tf_per_cluster$cluster) {
    lbl <- dominant_tf_per_cluster$label[dominant_tf_per_cluster$cluster == cid][1]
    if (!is.na(lbl) && nzchar(lbl)) return(lbl)
  }
  if (clustering_method == "2" && exists("tf_to_family")) {
    members_all <- names(cl[cl == cid])
    fam <- if (length(members_all) > 0) tf_to_family[[members_all[1]]] else NULL
    if (!is.null(fam) && length(fam) > 0 && !is.na(fam) && nzchar(fam)) return(fam)
  }
  members_all <- names(cl[cl == cid])
  counts <- sapply(members_all, function(tf) length(tf_targets[[tf]]))
  paste0(names(which.max(counts)), "_cluster")
}

cluster_order <- vapply(sort(unique(cl)), cluster_label, character(1))
cluster_order <- unique(cluster_order)

# ==============================================================================
# SECTION 3: BUILD BUBBLE DATA (one row per unit × gene list)
# ==============================================================================

fgsea_nes <- function(targets, mn_ranked, rgc_ranked) {
  pathway <- list(unit = targets)
  mn_res  <- tryCatch(fgsea::fgsea(pathways = pathway, stats = mn_ranked,
                                   minSize = 1, maxSize = 500, scoreType = "std"),
                      error = function(e) NULL)
  rgc_res <- tryCatch(fgsea::fgsea(pathways = pathway, stats = rgc_ranked,
                                   minSize = 1, maxSize = 500, scoreType = "std"),
                      error = function(e) NULL)
  mean(c(mn_res$NES, rgc_res$NES), na.rm = TRUE)
}

build_plot_df <- function(group_by = c("cluster", "tf")) {
  group_by <- match.arg(group_by)
  rows <- list()
  for (csv in list.files("filtered_genelists/graph_output", pattern = "\\.csv$", full.names = TRUE)) {
    name <- sub("\\_graph_output.csv$", "", basename(csv))

    obj_name <- paste0(name, "_cytoscape_output")
    if (!exists(obj_name)) {
      stop(sprintf("'%s' not found in environment. Run script 13 first.", obj_name))
    }
    cyto  <- get(obj_name)
    genelist <- readLines(file.path("filtered_genelists", paste0(name, "_genes.txt")))

    regs <- intersect(unique(cyto$name[cyto$`Regulatory function` == "Regulator"]),
                      names(cl))
    if (length(regs) == 0) next

    if (group_by == "cluster") {
      for (cid in sort(unique(cl[regs]))) {
        members <- regs[cl[regs] == cid]
        total   <- unique(unlist(tf_targets[members]))
        targets <- intersect(total, genelist)
        rows[[length(rows) + 1]] <- data.frame(
          gene_list = name,
          label     = cluster_label(cid),
          n_tfs     = length(members),
          n_targets = length(targets),
          pct_outdeg = if (length(total) > 0) length(targets) / length(total) * 100 else 0,
          NES       = if (length(targets) >= 5) fgsea_nes(targets, mn_ranked, rgc_ranked) else NA_real_,
          targets   = I(list(targets)),
          stringsAsFactors = FALSE
        )
      }
    } else {
      for (tf in regs) {
        total   <- unique(tf_targets[[tf]])
        targets <- intersect(total, genelist)
        rows[[length(rows) + 1]] <- data.frame(
          gene_list = name,
          label     = tf,
          n_tfs     = 1,
          n_targets = length(targets),
          pct_outdeg = if (length(total) > 0) length(targets) / length(total) * 100 else 0,
          NES       = if (length(targets) >= 5) fgsea_nes(targets, mn_ranked, rgc_ranked) else NA_real_,
          targets   = I(list(targets)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

plot_df_cluster <- build_plot_df("cluster")
plot_df_tf      <- build_plot_df("tf")

plot_df_cluster <- plot_df_cluster %>%
  dplyr::mutate(
    gene_list = factor(gene_list, levels = gene_list_order),
    label     = factor(label, levels = intersect(cluster_order, unique(label)))
  )
plot_df_tf <- plot_df_tf %>%
  dplyr::mutate(
    gene_list = factor(gene_list, levels = gene_list_order)
  )

# ==============================================================================
# SECTION 4: BUBBLE PLOT
# ==============================================================================

make_bubble_plot <- function(plot_df, title) {
  ggplot(plot_df, aes(x = gene_list, y = label, size = pct_outdeg, color = NES)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = n_targets), color = "black", size = 3, vjust = -1.2) +
    scale_color_viridis_c(option = "plasma", na.value = "grey60", name = "NES") +
    scale_size_continuous(range = c(2, 12), name = "% out-degree") +
    labs(x = NULL, y = NULL, title = title,
         subtitle = "Bubble size = % of cluster targets in gene list; color = avg NES (MN + RGC); number = target count") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black", face = "bold"))
}

p_cluster <- make_bubble_plot(plot_df_cluster,
                              title = "TF Clusters vs Filtered Gene Lists")
print(p_cluster)
ggsave(file.path(out_dir, "tf_families_by_genelist.png"), p_cluster,
       width = 10, height = max(6, nrow(plot_df_cluster) * 0.28), dpi = 300)

p_tf <- make_bubble_plot(plot_df_tf, title = "Individual TFs vs Filtered Gene Lists")
print(p_tf)
ggsave(file.path(out_dir, "tfs_by_genelist.png"), p_tf,
       width = 10, height = max(6, nrow(plot_df_tf) * 0.28), dpi = 300)

# ==============================================================================
# SECTION 5: SAVE DATA (targets stored per bubble)
# ==============================================================================

write_units <- function(plot_df, file_name) {
  out <- plot_df %>%
    dplyr::mutate(targets = sapply(targets, paste, collapse = "; "))
  write.csv(out, file.path(out_dir, file_name), row.names = FALSE)
}

write_units(plot_df_cluster, "tf_families_by_genelist.csv")
write_units(plot_df_tf, "tfs_by_genelist.csv")

cat("\nSaved bubble plots + data to", out_dir, "\n")
cat("\n=== Script 14 complete ===\n")
