# Script 11: Cluster Deep Dive
#
# Two modes:
#   1 = Run on ALL clusters (combined targets per cluster)
#   2 = Run on ALL TFs individually (every TF across all clusters)
#
# Outputs a single master table `deep_dive_units` with columns:
#   unit_id, label, type ("cluster" or "tf"), n_members, members, n_targets, targets
# (targets is a list-column of gene symbols)
#
# Then runs subnetwork + GO enrichment on each unit.
# Results saved to:
#   Mode 1: all_clusters_analysis/<label>/
#   Mode 2: all_tf_analysis/<label>/

library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(rrvgo)
if (exists("method") && method == "2") {
  if (!requireNamespace("simplifyEnrichment", quietly = TRUE)) {
    stop("simplifyEnrichment package required. Install with BiocManager::install('simplifyEnrichment')")
  }
  library(simplifyEnrichment)
}
library(org.Mm.eg.db)
library(GO.db)
library(beepr)
library(igraph)
library(RCy3)

source("utils/helpers.R")

# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

net <- load_tf_network("reformatted_edges.csv")
all_tfs <- net$all_tfs
tf_targets <- net$tf_targets

expr <- load_expression("Expression_data.txt")
exp_map <- expr$exp_map
mn_ranked  <- expr$mn_ranked
rgc_ranked <- expr$rgc_ranked

if (!file.exists("iregulon_results.csv")) {
  stop("'iregulon_results.csv' not found. Run iRegulon reformatting first.")
}
source("iregulon_reformat.R")
master_network <- reformat_networks("iregulon_results.csv", "intersected_upregulated_IDs.txt")
master_edges <- master_network$edges
master_nodes <- master_network$nodes

cat("Loaded", length(all_tfs), "TFs and", length(unique(unlist(tf_targets))), "unique targets\n")

# ==============================================================================
# SECTION 2: CYTOSCAPE CONNECTION CHECK
# ==============================================================================

run_cytoscape <- FALSE

tryCatch({
  RCy3::cytoscapePing()
  run_cytoscape <- TRUE
  message("Cytoscape connected.")
}, error = function(e) {
  message("Cytoscape not available. Subnetworks will be saved as files only.")
})

# ==============================================================================
# SECTION 3: DETERMINE CLUSTERING METHOD AND LIST CLUSTERS
# ==============================================================================

if (!exists("clustering_method")) {
  stop("'clustering_method' not found. Run script 10 first.")
}
if (!exists("cl") || is.null(cl)) {
  stop("'cl' (cluster assignments) not found. Run script 10 first.")
}

cat(sprintf("\nUsing clustering from script 10: %s\n", ifelse(clustering_method == "2", "JASPAR", "Jaccard")))

cluster_ids <- sort(unique(cl))
cat("\nAvailable clusters:\n")
for (i in seq_along(cluster_ids)) {
  cs <- cluster_ids[i]
  members <- names(cl[cl == cs])
  if (clustering_method == "2") {
    fam_name <- if (exists("tf_to_family")) tf_to_family[[members[1]]] else paste("Cluster", cs)
    label <- fam_name
  } else {
    tf_target_counts <- sapply(members, function(tf) length(tf_targets[[tf]]))
    dom_tf <- names(which.max(tf_target_counts))
    label <- paste0(dom_tf, "_cluster")
  }
  cat(sprintf("  %d. %s (%d TFs): %s\n", i, label, length(members), paste(members, collapse = ", ")))
}

# ==============================================================================
# SECTION 4: BUILD MASTER TABLE
# ==============================================================================

deep_dive_units <- data.frame(
  unit_id          = integer(),
  label            = character(),
  type             = character(),
  n_members        = integer(),
  members          = I(list()),
  n_targets        = integer(),
  targets          = I(list()),
  n_targets_unique  = integer(),
  targets_unique   = I(list()),
  parent_terms     = I(list()),
  stringsAsFactors = FALSE
)

if (!exists("skip_prompts") || !skip_prompts) {
  beep(10)
  cat("\nAnalyze:\n")
  cat("  1 = ALL clusters (combined targets per cluster)\n")
  cat("  2 = ALL TFs individually (every TF across all clusters)\n")
  mode_input <- readline(prompt = "Select mode (1 or 2): ")
  analysis_mode <- trimws(mode_input)
} else {
  if (!exists("analysis_mode")) analysis_mode <- "1"
  cat(sprintf("Using pre-set mode: %s\n", analysis_mode))
}

if (!analysis_mode %in% c("1", "2")) {
  stop("Invalid mode. Enter 1 or 2.")
}

# Wipe old output directories
if (analysis_mode == "1") {
  unlink("all_clusters_analysis", recursive = TRUE)
} else {
  unlink("all_tf_analysis", recursive = TRUE)
}

if (analysis_mode == "1") {
  cat("\nBuilding cluster-level master table...\n")
  for (i in seq_along(cluster_ids)) {
    cs <- cluster_ids[i]
    members <- names(cl[cl == cs])
    if (clustering_method == "2") {
      lbl <- if (exists("tf_to_family")) tf_to_family[[members[1]]] else paste("Cluster", cs)
    } else {
      tf_target_counts <- sapply(members, function(tf) length(tf_targets[[tf]]))
      lbl <- paste0(names(which.max(tf_target_counts)), "_cluster")
    }
    targs <- unique(unlist(tf_targets[members]))
    deep_dive_units <- rbind(deep_dive_units, data.frame(
      unit_id = i,
      label   = lbl,
      type    = "cluster",
      n_members = length(members),
      members = I(list(members)),
      n_targets = length(targs),
      targets = I(list(targs)),
      stringsAsFactors = FALSE
    ))
  }
} else {
  cat("\nBuilding TF-level master table...\n")
  all_tfs_in_clusters <- sort(unique(unlist(sapply(cluster_ids, function(cs) names(cl[cl == cs])))))
  for (i in seq_along(all_tfs_in_clusters)) {
    tf <- all_tfs_in_clusters[i]
    targs <- tf_targets[[tf]]
    deep_dive_units <- rbind(deep_dive_units, data.frame(
      unit_id = i,
      label   = tf,
      type    = "tf",
      n_members = 1,
      members = I(list(tf)),
      n_targets = length(targs),
      targets = I(list(targs)),
      stringsAsFactors = FALSE
    ))
  }
}

# Compute targets unique to each unit (not in any other unit)
all_target_lists <- deep_dive_units$targets
for (k in seq_len(nrow(deep_dive_units))) {
  others <- unique(unlist(all_target_lists[-k]))
  uniq <- setdiff(deep_dive_units$targets[[k]], others)
  deep_dive_units$targets_unique[[k]] <- uniq
  deep_dive_units$n_targets_unique[[k]] <- length(uniq)
}
# Convert n_targets_unique from list to vector
deep_dive_units$n_targets_unique <- as.integer(deep_dive_units$n_targets_unique)

cat(sprintf("Master table: %d units\n", nrow(deep_dive_units)))

# ==============================================================================
# SECTION 5: BUILD GO:BP GENE SETS
# ==============================================================================

cat("Building GO:BP gene sets...\n")

go_annots <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = keys(org.Mm.eg.db, keytype = "SYMBOL"),
                                   columns = c("GO", "SYMBOL"),
                                   keytype = "SYMBOL") %>%
  dplyr::filter(ONTOLOGY == "BP", !is.na(GO), !is.na(SYMBOL))

go_ids <- unique(go_annots$GO)
go_term_lookup <- AnnotationDbi::select(GO.db, keys = go_ids, columns = "TERM", keytype = "GOID")
go_term_lookup <- setNames(go_term_lookup$TERM, go_term_lookup$GOID)

go_gene_sets_raw <- go_annots %>%
  dplyr::group_by(GO) %>%
  dplyr::summarise(genes = list(unique(SYMBOL)), .groups = "drop") %>%
  dplyr::filter(lengths(genes) >= 5)

go_gene_list <- setNames(go_gene_sets_raw$genes, go_gene_sets_raw$GO)
cat(sprintf("  Loaded %d GO:BP gene sets\n", length(go_gene_list)))

if (!exists("method")) method <- "1"
if (method == "1") {
  if (!exists("rrvgo_threshold")) rrvgo_threshold <- 0.7
  cat(sprintf("GO collapsing method: rrvgo (threshold=%.2f)\n", rrvgo_threshold))
} else {
  if (!exists("chosen_threshold")) chosen_threshold <- "0.80"
  cat(sprintf("GO collapsing method: simplifyEnrichment (threshold=%s)\n", chosen_threshold))
}

# ==============================================================================
# SECTION 6: HELPERS
# ==============================================================================

run_go_fgsea <- function(target_genes, mn_ranked, rgc_ranked,
                        gene_sets = go_gene_list, term_lookup = go_term_lookup) {
  if (length(target_genes) < 5) {
    message("  Fewer than 5 target genes, skipping fgsea.")
    return(NULL)
  }
  tryCatch({
    avg_ranked <- sort((mn_ranked + rgc_ranked[names(mn_ranked)]) / 2, decreasing = TRUE)
    avg_ranked <- avg_ranked[!is.na(avg_ranked)]
    pathways <- lapply(gene_sets, function(gs) intersect(gs, target_genes))
    pathways <- pathways[lengths(pathways) >= 1]
    if (length(pathways) == 0) {
      message("  No GO terms with target overlap for fgsea.")
      return(NULL)
    }
    fgsea_res <- fgsea::fgsea(pathways = pathways, stats = avg_ranked,
                              minSize = 1, maxSize = 500, scoreType = "std")
    if (nrow(fgsea_res) == 0 || all(is.na(fgsea_res$pval))) {
      message("  No significant fgsea results.")
      return(NULL)
    }
    fgsea_res %>%
      dplyr::filter(!is.na(pval), pval < 0.05) %>%
      dplyr::mutate(
        term_id = pathway,
        term_name = term_lookup[pathway],
        neg_log10_p = -log10(pval)
      ) %>%
      dplyr::arrange(pval) %>%
      dplyr::select(term_id, term_name, pval, NES, size, neg_log10_p)
  }, error = function(e) {
    message("  fgsea failed: ", e$message)
    return(NULL)
  })
}

run_go_gsea <- function(target_genes, mn_ranked, rgc_ranked,
                        gene_sets = go_gene_list, term_lookup = go_term_lookup) {
  if (length(target_genes) < 5) {
    message("  Fewer than 5 target genes, skipping GSEA.")
    return(NULL)
  }
  tryCatch({
    avg_ranked <- sort((mn_ranked + rgc_ranked[names(mn_ranked)]) / 2, decreasing = TRUE)
    avg_ranked <- avg_ranked[!is.na(avg_ranked)]
    pathways <- lapply(gene_sets, function(gs) intersect(gs, target_genes))
    pathways <- pathways[lengths(pathways) >= 2]
    if (length(pathways) == 0) {
      message("  No GO terms with sufficient target overlap for GSEA.")
      return(NULL)
    }
    fgsea_res <- fgsea::fgsea(pathways = pathways, stats = avg_ranked,
                              minSize = 2, maxSize = 500, scoreType = "std")
    if (nrow(fgsea_res) == 0 || all(is.na(fgsea_res$pval))) {
      message("  No significant GSEA results.")
      return(NULL)
    }
    fgsea_res %>%
      dplyr::filter(!is.na(pval), pval < 0.05) %>%
      dplyr::mutate(
        term_id = pathway,
        term_name = term_lookup[pathway],
        neg_log10_p = -log10(pval)
      ) %>%
      dplyr::arrange(pval) %>%
      dplyr::select(term_id, term_name, pathway, pval, NES, size, neg_log10_p)
  }, error = function(e) {
    message("  GSEA failed: ", e$message)
    return(NULL)
  })
}

reduce_go_results <- function(go_df, score_col = "neg_log10_p") {
  if (is.null(go_df) || nrow(go_df) < 3) return(go_df)
  term_id_col <- if ("term_id" %in% names(go_df)) "term_id" else "pathway"
  term_ids <- go_df[[term_id_col]]
  scores <- setNames(go_df[[score_col]], term_ids)

  if (method == "2") {
    tryCatch({
      cat(sprintf("    Computing GO similarity matrix (simplifyEnrichment) for %d terms...\n", length(term_ids)))
      mat <- GO_similarity(term_ids, db = "org.Mm.eg.db")
      hc <- hclust(as.dist(1 - mat), method = "average")
      clusters <- cutree(hc, h = as.numeric(chosen_threshold))
      cat(sprintf("    Cut at h=%s -> %d clusters\n", chosen_threshold, length(unique(clusters))))
      cluster_df <- data.frame(term_id = term_ids, cluster_id = clusters, score = scores[term_ids], stringsAsFactors = FALSE)
      parents <- cluster_df %>% dplyr::group_by(cluster_id) %>% dplyr::slice_max(order_by = score, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
      parent_map <- setNames(parents$term_id, parents$cluster_id)
      go_df$parentTerm <- go_df$term_name[match(parent_map[as.character(clusters[go_df[[term_id_col]]])], go_df[[term_id_col]])]
      go_df$parentTerm[is.na(go_df$parentTerm)] <- go_df$term_name[is.na(go_df$parentTerm)]
      go_df %>% dplyr::group_by(parentTerm) %>% dplyr::slice_max(order_by = abs(.data[[score_col]]), n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::arrange(pval)
    }, error = function(e) {
      message("    simplifyEnrichment failed: ", e$message, ". Using unreduced terms.")
      go_df
    })
  } else {
    tryCatch({
      cat(sprintf("    Computing similarity matrix (rrvgo) for %d terms...\n", length(term_ids)))
      sim_matrix <- calculateSimMatrix(term_ids, orgdb = "org.Mm.eg.db", ont = "BP", method = "Rel")
      reduced <- reduceSimMatrix(sim_matrix, scores, threshold = rrvgo_threshold, orgdb = "org.Mm.eg.db")
      parent_map <- setNames(reduced$parentTerm, reduced$go)
      go_df$parentTerm <- parent_map[go_df[[term_id_col]]]
      go_df$parentTerm[is.na(go_df$parentTerm)] <- go_df[[term_id_col]][is.na(go_df$parentTerm)]
      go_df %>% dplyr::group_by(parentTerm) %>% dplyr::slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::arrange(pval)
    }, error = function(e) {
      message("    rrvgo failed: ", e$message, ". Using unreduced terms.")
      go_df
    })
  }
}

build_subnetwork <- function(target_genes, regulators = NULL) {
  is_ensembl <- any(grepl("^ENS", target_genes))
  if (is_ensembl) {
    target_symbols <- exp_map %>% dplyr::filter(ensembl_id %in% target_genes) %>% dplyr::pull(symbol)
  } else {
    target_symbols <- target_genes
  }
  sub_edges <- master_edges %>% dplyr::filter(`Target Gene` %in% target_symbols)
  if (!is.null(regulators)) {
    sub_edges <- sub_edges %>% dplyr::filter(`Regulator Gene` %in% regulators)
  }
  if (nrow(sub_edges) == 0) return(NULL)
  active_regulators <- unique(sub_edges$`Regulator Gene`)
  active_targets <- unique(sub_edges$`Target Gene`)
  all_nodes <- unique(c(active_regulators, active_targets))
  sub_nodes <- master_nodes %>% dplyr::filter(id %in% all_nodes)
  g <- graph_from_data_frame(d = sub_edges, directed = TRUE)
  in_deg  <- degree(g, mode = "in")
  out_deg <- degree(g, mode = "out")
  bet     <- betweenness(g, directed = TRUE, normalized = FALSE)
  close   <- closeness(g, mode = "out")
  metrics <- data.frame(
    name = names(in_deg),
    InDegree = as.numeric(in_deg),
    OutDegree = as.numeric(out_deg),
    BetweennessCentrality = as.numeric(bet),
    ClosenessCentrality = as.numeric(close),
    stringsAsFactors = FALSE
  )
  list(edges = sub_edges, nodes = sub_nodes, metrics = metrics,
       regulators = active_regulators, targets = active_targets)
}

send_to_cytoscape <- function(sub, net_title) {
  cyto_edges <- sub$edges %>%
    dplyr::rename(source = `Regulator Gene`, target = `Target Gene`) %>%
    dplyr::mutate(source = as.character(source), target = as.character(target))
  cyto_nodes <- sub$nodes %>% dplyr::mutate(id = as.character(id)) %>% dplyr::select(id, everything()) %>% as.data.frame()
  captured <- tryCatch({
    capture.output(suppressMessages(RCy3::createNetworkFromDataFrames(nodes = cyto_nodes, edges = cyto_edges, title = net_title, collection = "Cluster_DeepDive", key.col.index = 1)))
    TRUE
  }, error = function(e) { message("  Cytoscape network creation failed: ", e$message); FALSE })
  if (captured) {
    tryCatch({ capture.output(suppressMessages(RCy3::loadTableData(data = sub$metrics, data.key.column = "name", table = "node", table.key.column = "name"))) }, error = function(e) NULL)
    tryCatch({ capture.output(suppressMessages(RCy3::setVisualStyle("GRN"))) }, error = function(e) warning("Could not apply GRN style to ", net_title, ": ", e$message))
  }
}

make_go_dotplot <- function(go_result, title = "GO:BP fgsea Enrichment", p_cutoff = 0.05, term_col = "term_name") {
  if (is.null(go_result) || nrow(go_result) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant GO:BP terms", size = 5) + theme_void() + labs(title = title))
  }
  plot_df <- go_result %>% dplyr::filter(pval < p_cutoff) %>% dplyr::arrange(dplyr::desc(NES)) %>% dplyr::mutate(term_label = reorder(.data[[term_col]], NES), neg_log10_p = -log10(pval))
  ggplot(plot_df, aes(x = term_label, y = NES)) +
    geom_point(aes(color = neg_log10_p), alpha = 0.85, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_viridis_c(option = "plasma", name = "-log10(p)") +
    coord_flip() +
    labs(title = paste0(title, " (", nrow(plot_df), " terms)"), x = NULL, y = "NES") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(color = "black", size = 7), axis.text.x = element_text(color = "black"), panel.grid.minor = element_blank(), legend.position = "right")
}

make_go_gsea_dotplot <- function(gsea_result, title = "GO:BP GSEA Enrichment", p_cutoff = 0.05, term_col = "term_name") {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant GSEA terms", size = 5) + theme_void() + labs(title = title))
  }
  plot_df <- gsea_result %>% dplyr::filter(pval < p_cutoff) %>% dplyr::arrange(dplyr::desc(NES)) %>% dplyr::mutate(term_label = reorder(.data[[term_col]], NES), neg_log10_p = -log10(pval))
  ggplot(plot_df, aes(x = term_label, y = NES)) +
    geom_point(aes(color = neg_log10_p), alpha = 0.85, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_viridis_c(option = "plasma", name = "-log10(p)") +
    coord_flip() +
    labs(title = paste0(title, " (", nrow(plot_df), " terms)"), x = NULL, y = "NES") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(color = "black", size = 7), axis.text.x = element_text(color = "black"), panel.grid.minor = element_blank(), legend.position = "right")
}

build_gene_process_table <- function(target_genes, fgsea_result, top_n_processes = 15, term_lookup = go_term_lookup) {
  if (is.null(fgsea_result) || nrow(fgsea_result) == 0) return(NULL)
  top_terms <- fgsea_result %>% dplyr::arrange(pval)
  process_names <- setNames(top_terms$term_name, top_terms$term_id)
  all_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
  go_annots <- AnnotationDbi::select(org.Mm.eg.db, keys = all_symbols[all_symbols %in% target_genes], columns = c("GO", "SYMBOL"), keytype = "SYMBOL") %>% dplyr::filter(ONTOLOGY == "BP", !is.na(GO))
  gene_process <- data.frame(symbol = unique(target_genes), stringsAsFactors = FALSE)
  for (i in seq_along(process_names)) {
    go_id <- names(process_names)[i]
    proc_name <- process_names[i]
    col_name <- gsub("[^[:alnum:]_]", "_", substr(proc_name, 1, 50))
    genes_in_proc <- go_annots %>% dplyr::filter(GO == go_id) %>% dplyr::pull(SYMBOL)
    gene_process[[col_name]] <- ifelse(gene_process$symbol %in% genes_in_proc, 1, 0)
  }
  gene_process$top_processes <- sapply(gene_process$symbol, function(g) {
    procs <- go_annots %>% dplyr::filter(SYMBOL == g, GO %in% names(process_names)) %>% dplyr::pull(GO)
    if (length(procs) == 0) return("")
    paste(process_names[procs], collapse = "; ")
  })
  process_cols <- names(gene_process)[!names(gene_process) %in% c("symbol", "top_processes")]
  gene_process$n_processes <- rowSums(gene_process[, process_cols, drop = FALSE])
  gene_process %>% dplyr::arrange(dplyr::desc(n_processes))
}

# ==============================================================================
# SECTION 7: RUN ANALYSIS ON EACH UNIT
# ==============================================================================

run_unit_analysis <- function(target_genes, label, out_dir, regulators = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^[:alnum:]_]", "_", label)

  cat("  Building subnetwork...\n")
  sub <- build_subnetwork(target_genes, regulators = regulators)
  if (!is.null(sub)) {
    write.table(sub$edges, file.path(out_dir, "edges.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    analyzed_nodes <- sub$nodes %>% dplyr::rename(name = id) %>% dplyr::left_join(sub$metrics, by = "name")
    write.table(analyzed_nodes, file.path(out_dir, "nodes.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    writeLines(sort(sub$regulators), file.path(out_dir, "regulators.txt"))
    cat(sprintf("  Subnetwork: %d edges, %d nodes (%d regulators, %d targets)\n",
                nrow(sub$edges), nrow(analyzed_nodes), length(sub$regulators), length(sub$targets)))
    if (run_cytoscape) {
      cat("  Sending to Cytoscape...\n")
      send_to_cytoscape(sub, paste0(safe_label, "_network"))
    }
    comp_df <- data.frame(Role = c("Regulator", "Regulated"), Count = c(length(sub$regulators), length(sub$targets)))
    p_comp <- ggplot(comp_df, aes(x = Role, y = Count, fill = Role)) +
      geom_col(width = 0.5) + geom_text(aes(label = Count), vjust = -0.3, size = 4) +
      scale_fill_manual(values = c("Regulator" = "#00BFC4", "Regulated" = "#F8766D")) +
      labs(title = paste("Network Composition —", label), y = "Number of Genes", x = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "none", panel.grid.major.x = element_blank())
    print(p_comp)
    ggsave(file.path(out_dir, "network_composition.png"), p_comp, width = 6, height = 4, dpi = 300)
  } else {
    cat("  No edges found in master network for these targets.\n")
    p_comp <- NULL
  }

  cat(sprintf("  Running GO:BP fgsea enrichment on %d target genes...\n", length(target_genes)))
  go_fgsea_raw <- run_go_fgsea(target_genes, mn_ranked, rgc_ranked)
  if (!is.null(go_fgsea_raw)) {
    write.csv(go_fgsea_raw %>% as.data.frame(), file.path(out_dir, "go_results_fgsea.csv"), row.names = FALSE)
    cat(sprintf("  fgsea: %d significant GO:BP terms\n", nrow(go_fgsea_raw)))
  }

  cat(sprintf("  Running GO:BP GSEA enrichment on %d target genes...\n", length(target_genes)))
  go_gsea_raw <- run_go_gsea(target_genes, mn_ranked, rgc_ranked)
  if (!is.null(go_gsea_raw)) {
    write.csv(go_gsea_raw %>% as.data.frame(), file.path(out_dir, "go_results_gsea.csv"), row.names = FALSE)
    cat(sprintf("  GSEA: %d significant GO:BP terms\n", nrow(go_gsea_raw)))
  }

  collapse_method <- ifelse(method == "1", sprintf("rrvgo (threshold=%.2f)", rrvgo_threshold), sprintf("simplifyEnrichment (threshold=%s)", chosen_threshold))
  cat(sprintf("  Collapsing redundant GO terms via %s...\n", collapse_method))
  go_fgsea <- reduce_go_results(go_fgsea_raw)
  go_gsea <- reduce_go_results(go_gsea_raw, score_col = "neg_log10_p")
  if (!is.null(go_fgsea)) cat(sprintf("  fgsea: %d -> %d parent terms\n", nrow(go_fgsea_raw), nrow(go_fgsea)))
  if (!is.null(go_gsea)) cat(sprintf("  GSEA: %d -> %d parent terms\n", nrow(go_gsea_raw), nrow(go_gsea)))

  p_fgsea <- make_go_dotplot(go_fgsea, title = paste("GO:BP fgsea —", label), term_col = "parentTerm")
  print(p_fgsea)
  n_fgsea <- if (!is.null(go_fgsea)) nrow(go_fgsea) else 0
  ggsave(file.path(out_dir, "go_dotplot_fgsea.png"), p_fgsea, width = 14, height = max(8, n_fgsea * 0.35), dpi = 300)

  p_gsea <- make_go_gsea_dotplot(go_gsea, title = paste("GO:BP GSEA —", label), term_col = "parentTerm")
  print(p_gsea)
  n_gsea <- if (!is.null(go_gsea)) nrow(go_gsea) else 0
  ggsave(file.path(out_dir, "go_dotplot_gsea.png"), p_gsea, width = 14, height = max(8, n_gsea * 0.35), dpi = 300)

  cat("  Building gene-process annotation table...\n")
  gene_proc <- build_gene_process_table(target_genes, go_fgsea_raw)
  if (!is.null(gene_proc)) {
    write.csv(gene_proc, file.path(out_dir, "gene_process_annotations.csv"), row.names = FALSE)
    cat(sprintf("  %d genes tagged with top GO processes\n", nrow(gene_proc)))
  }

  pdf(file.path(out_dir, "report.pdf"), width = 12, height = 8)
  if (!is.null(p_comp)) print(p_comp)
  print(p_fgsea)
  print(p_gsea)
  dev.off()

  cat(sprintf("  Outputs saved to %s/\n", out_dir))
  
  # Collect significant parent terms from both fgsea and GSEA
  all_parents <- character()
  if (!is.null(go_fgsea) && "parentTerm" %in% names(go_fgsea)) {
    all_parents <- c(all_parents, go_fgsea$parentTerm)
  }
  if (!is.null(go_gsea) && "parentTerm" %in% names(go_gsea)) {
    all_parents <- c(all_parents, go_gsea$parentTerm)
  }
  sort(unique(all_parents))
}

cat(sprintf("\n=== Running analysis on %d units ===\n", nrow(deep_dive_units)))

for (i in seq_len(nrow(deep_dive_units))) {
  unit <- deep_dive_units[i, ]
  cat(sprintf("\n--- Unit %d/%d: %s (%s, %d targets) ---\n",
              i, nrow(deep_dive_units), unit$label, unit$type, unit$n_targets))
  out_dir <- if (unit$type == "cluster") {
    file.path("all_clusters_analysis", gsub("[^[:alnum:]_]", "_", unit$label))
  } else {
    file.path("all_tf_analysis", gsub("[^[:alnum:]_]", "_", unit$label))
  }
  parent_terms <- run_unit_analysis(unit$targets[[1]], unit$label, out_dir, regulators = unit$members[[1]])
  deep_dive_units$parent_terms[[i]] <- parent_terms
}

cat(sprintf("\n=== All %d units complete ===\n", nrow(deep_dive_units)))
cat("\n=== Script 11 complete ===\n")
