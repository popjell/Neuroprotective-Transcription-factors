# Standalone community detection analysis
# Clusters TF-target network via Leiden, labels communities via GO,
# computes TF activity, and exports a dot plot.
#
# Usage:
#   1. Run the main pipeline first (or ensure reformatted_edges.csv + Expression_data.txt exist in root)
#   2. Rscript community_detection/run_community.R
#
# Input files needed (from project root):
#   - reformatted_edges.csv  (TF -> target gene edges)
#   - Expression_data.txt    (expression data for enrichment scoring)

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)

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
colnames(edges) <- gsub("\\.", " ", colnames(edges))

expr <- load_expression(exp_path)
exp_map <- expr$exp_map
mn_ranked <- expr$mn_ranked
rgc_ranked <- expr$rgc_ranked

all_tfs <- unique(edges$`Regulator Gene`)
tf_targets <- list()
for (tf in all_tfs) {
  tf_targets[[tf]] <- unique(edges$`Target Gene`[edges$`Regulator Gene` == tf])
}

cat("\n=== COMMUNITY DETECTION ===\n")
cat("TFs:", length(all_tfs), "\n")
cat("Targets:", length(unique(edges$`Target Gene`)), "\n")
cat("Edges:", nrow(edges), "\n")

# ==============================================================================
# SECTION 2: COMMUNITY DETECTION
# ==============================================================================

g <- graph_from_data_frame(d = edges[, c("Regulator Gene", "Target Gene")], directed = TRUE)
g_undir <- as_undirected(g, mode = "collapse")

cat("\n--- Resolution sweep ---\n")
for (res in c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) {
  set.seed(42)
  result <- cluster_leiden(g_undir, objective_function = "modularity", resolution = res)
  cat("  Resolution", res, "->", length(unique(membership(result))), "communities\n")
}

cat("\nEnter resolution parameter (default 1.0): ")
chosen_res <- readline()
chosen_res <- trimws(chosen_res)
if (chosen_res == "") chosen_res <- 1.0 else chosen_res <- as.numeric(chosen_res)
cat("Using resolution =", chosen_res, "\n")

set.seed(42)
leiden_result <- cluster_leiden(g_undir, objective_function = "modularity", resolution = chosen_res)
mem <- membership(leiden_result)

comm_df <- data.frame(
  node = names(mem),
  community = as.integer(mem),
  stringsAsFactors = FALSE
)
all_targets <- unique(edges$`Target Gene`)
comm_df$is_target <- comm_df$node %in% all_targets

cat("\nCommunity sizes:\n")
print(table(comm_df$community))

# ==============================================================================
# SECTION 3: LABEL COMMUNITIES VIA GO ENRICHMENT
# ==============================================================================

cat("\n--- Labeling communities via GO enrichment ---\n")

label_df <- data.frame(community = integer(), label = character(), stringsAsFactors = FALSE)

for (cl in sort(unique(comm_df$community))) {
  genes <- comm_df$node[comm_df$community == cl & comm_df$is_target]
  n_total <- sum(comm_df$community == cl)

  if (length(genes) < 2) {
    lbl <- paste0("C", cl, " (", n_total, " genes)")
    label_df <- rbind(label_df, data.frame(community = cl, label = lbl, stringsAsFactors = FALSE))
    cat("  C", cl, ": too few targets for enrichment ->", lbl, "\n", sep = "")
    next
  }

  ego <- tryCatch({
    enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
             qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 500)
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    top_term <- as.data.frame(ego)$Description[1]
    lbl <- paste0("C", cl, ": ", top_term)
  } else {
    lbl <- paste0("C", cl, " (", length(genes), " target genes)")
  }
  label_df <- rbind(label_df, data.frame(community = cl, label = lbl, stringsAsFactors = FALSE))
  cat("  ", lbl, "\n")
}

# ==============================================================================
# SECTION 4: TF ACTIVITY PER COMMUNITY
# ==============================================================================

cat("\n--- Computing TF activity per community ---\n")

comm_rows <- list()

for (cl in sort(unique(comm_df$community))) {
  comm_genes <- comm_df$node[comm_df$community == cl]
  comm_edges <- edges %>% filter(`Target Gene` %in% comm_genes)
  active_tfs <- unique(comm_edges$`Regulator Gene`)

  pathways <- list()
  fallback <- list()
  for (tf in active_tfs) {
    targs <- unique(comm_edges$`Target Gene`[comm_edges$`Regulator Gene` == tf])
    if (length(targs) >= 2) pathways[[tf]] <- targs else fallback[[tf]] <- targs
  }

  cl_rows <- list()

  if (length(pathways) > 0) {
    mn_res <- fgsea(pathways = pathways, stats = mn_ranked,
                     minSize = 1, maxSize = 500, scoreType = "std")
    rgc_res <- fgsea(pathways = pathways, stats = rgc_ranked,
                     minSize = 1, maxSize = 500, scoreType = "std")

    for (tf in names(pathways)) {
      n_targ <- length(pathways[[tf]])
      nes_mn <- mn_res$NES[mn_res$pathway == tf]
      nes_rgc <- rgc_res$NES[rgc_res$pathway == tf]
      if (length(nes_mn) == 0) nes_mn <- NA
      if (length(nes_rgc) == 0) nes_rgc <- NA
      combined <- mean(c(nes_mn, nes_rgc), na.rm = TRUE)
      if (is.na(combined)) next

      pct <- n_targ / length(comm_genes) * 100

      cl_rows[[tf]] <- data.frame(
        Community = cl,
        TF = tf,
        enrichment_score = combined,
        NES_MN = if (length(nes_mn) > 0) nes_mn else NA,
        NES_RGC = if (length(nes_rgc) > 0) nes_rgc else NA,
        n_targets = n_targ,
        pct_community = pct,
        method = "fgsea",
        stringsAsFactors = FALSE
      )
    }
  }

  for (tf in names(fallback)) {
    targs <- fallback[[tf]]
    n_targ <- length(targs)
    if (n_targ == 0) next

    target_stats <- exp_map %>% filter(symbol %in% targs)
    mean_mn <- mean(-log10(target_stats$pvalue_MN) * sign(target_stats$log2FoldChange_MN), na.rm = TRUE)
    mean_rgc <- mean(-log10(target_stats$pvalue_RGC) * sign(target_stats$log2FoldChange_RGC), na.rm = TRUE)
    combined <- mean(c(mean_mn, mean_rgc), na.rm = TRUE)
    if (is.na(combined)) next

    pct <- n_targ / length(comm_genes) * 100
    cl_rows[[tf]] <- data.frame(
      Community = cl,
      TF = tf,
      enrichment_score = combined,
      NES_MN = mean_mn,
      NES_RGC = mean_rgc,
      n_targets = n_targ,
      pct_community = pct,
      method = "fallback",
      stringsAsFactors = FALSE
    )
  }

  comm_rows[[as.character(cl)]] <- bind_rows(cl_rows)
}

tf_community_df <- bind_rows(comm_rows)

if (nrow(tf_community_df) > 0) {
  tf_community_df <- tf_community_df %>%
    group_by(TF) %>%
    mutate(specificity = n_targets / sum(n_targets) * 100) %>%
    ungroup()
}

tf_community_df <- tf_community_df %>%
  left_join(label_df, by = c("Community" = "community"))

write.csv(tf_community_df, "tf_community_enrichment.csv", row.names = FALSE)
cat("\nExported tf_community_enrichment.csv\n")

# ==============================================================================
# SECTION 5: DOT PLOT
# ==============================================================================

cl_order <- tf_community_df %>%
  group_by(Community, label) %>%
  summarise(avg_es = mean(enrichment_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(avg_es) %>%
  pull(label)

tf_community_df <- tf_community_df %>%
  mutate(label = factor(label, levels = cl_order))

tf_order <- tf_community_df %>%
  group_by(TF) %>%
  summarise(avg_es = mean(enrichment_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_es)) %>%
  pull(TF)

tf_community_df <- tf_community_df %>%
  mutate(TF = factor(TF, levels = tf_order))

p <- ggplot(tf_community_df, aes(
  x = label, y = TF,
  size = specificity,
  color = enrichment_score
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", limits = c(0, 5), oob = scales::squish) +
  scale_size_continuous(range = c(2, 10)) +
  theme_bw(base_size = 11) +
  labs(
    title = "TF Activity Across Gene Communities",
    subtitle = paste0("Leiden clustering (resolution = ", chosen_res,
                      ")\nDot size = % of each TF's targets in that community"),
    x = "Community (top GO term)",
    y = "Transcription Factor",
    size = "% of TF targets",
    color = "Enrichment\nScore (NES)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
    axis.text.y = element_text(color = "black", face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )
p
name = paste0("community_detection/tf_community_dotplot_res:", chosen_res, ".png")

if(file.exists(name)) {
  file.remove(name)
}
ggsave(name, p, width = 12, height = 10, dpi = 300)
cat("Saved", name, "\n")

cat("\n=== DONE ===\n")
cat("Output files:\n")
cat("  tf_community_dotplot.png  - dot plot\n")
cat("  tf_community_enrichment.csv  - enrichment table\n")
