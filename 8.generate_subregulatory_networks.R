#NOTE: BEFORE RUNNING THIS PART, MAKE SURE YOU HAVE:
#1. imported the intersected_upregulated_genes.txt into cytoscape as a network, ran iregulon, and exported it as a .tsv file
#2. converted the .tsv file into a .csv file and placed it in the root directory with the name "iregulon_results.csv"

library(dplyr)
library(tidyr)
library(fs)
library(ggplot2)
library(RCy3)
library(igraph)
library(fgsea)

source("utils/helpers.R")

# ==============================================================================
# SECTION 1: MAP GENES TO GO PARENT CATEGORIES
# ==============================================================================

go_genes_long <- go_analysis %>%
  dplyr::filter(source == "GO:BP") %>%
  dplyr::select(term_id = term_id, term_name, intersections) %>%
  separate_rows(intersections, sep = ",\\s*|\\s+") %>% 
  dplyr::rename(gene = intersections)

genes_with_parents <- reducedTerms %>%
  dplyr::select(term = go, parentTerm) %>%
  inner_join(go_genes_long, by = c("term" = "term_id"))

parent_gene_lists <- genes_with_parents %>%
  group_by(parentTerm) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  tibble::deframe()

# --- Optional: use only genes unique to a single parent category ---
gene_parent_count <- genes_with_parents %>%
  distinct(gene, parentTerm) %>%
  dplyr::count(gene, name = "n_parents")

unique_genes <- gene_parent_count %>%
  dplyr::filter(n_parents == 1) %>%
  pull(gene)

parent_unique_gene_lists <- parent_gene_lists
for (cat in names(parent_unique_gene_lists)) {
  parent_unique_gene_lists[[cat]] <- intersect(parent_unique_gene_lists[[cat]], unique_genes)
}
parent_unique_gene_lists <- parent_unique_gene_lists[sapply(parent_unique_gene_lists, length) > 0]

beep(10)

user_input <- readline(prompt = "Use only genes unique to a single parent category? (y/n): ")
use_unique <- tolower(user_input) %in% c("y", "yes")

if (use_unique) {
  parent_gene_lists <- parent_unique_gene_lists
  message(sprintf("Using unique-gene mode: %d categories, %d total unique genes",
                  length(parent_gene_lists), length(unique(unlist(parent_gene_lists)))))
} else {
  message(sprintf("Using all-gene mode: %d categories, %d total unique genes",
                  length(parent_gene_lists), length(unique(unlist(parent_gene_lists)))))
}
# -------------------------------------------------------------------

output_dir <- "parent_category_genes"
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
dir.create(output_dir)

for (category_name in names(parent_gene_lists)) {
  safe_name <- gsub("[^[:alnum:]_]", "_", category_name)
  file_path <- file.path(output_dir, paste0(safe_name, "_genes.txt"))
  genes_to_write <- parent_gene_lists[[category_name]]
  write(genes_to_write, file = file_path)
}

print(paste("Successfully exported", length(parent_gene_lists), "files to the folder:", output_dir))



# ==============================================================================
# SECTION 2: IMPORT THE WHOLE-GENOME MASTER NETWORKS & USER INTERACTIVE PROMPT
# ==============================================================================
if (!file.exists("iregulon_results.csv")) {
  beep(10)
  stop("'iregulon_results.csv' not found in the root directory. follow instructions in the comments of the script 7")
}


master_network <- reformat_networks("iregulon_results.csv", "intersected_upregulated_IDs.txt")
master_edges <- master_network$edges
master_nodes <- master_network$nodes

base_output_dir = "parent_subnetworks"
if (dir.exists(base_output_dir)) unlink(base_output_dir, recursive = TRUE)
dir.create(base_output_dir)

exp_map <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Build preranked gene lists for GSEA enrichment analysis
# (build_ranked is defined in utils/helpers.R, sourced above)
mn_ranked   <- build_ranked(exp_map$pvalue_MN, exp_map$log2FoldChange_MN, exp_map$symbol)
rgc_ranked  <- build_ranked(exp_map$pvalue_RGC, exp_map$log2FoldChange_RGC, exp_map$symbol)

# INTERACTIVE PROMPT: Ask the user whether to run the heavy Cytoscape steps
if (!exists("run_cytoscape")) {
  beep(10)
  user_input <- readline(prompt = "Run cytoscape analysis as well? (y/n): ")
  run_cytoscape <- lgl <- tolower(user_input) %in% c("y", "yes")
}


# If the user selected yes, make sure Cytoscape is actually open before starting
if (run_cytoscape) {
  is_connected <- tryCatch({
    RCy3::cytoscapePing()
    TRUE
  }, error = function(e) {
      FALSE
  })

  if (is_connected) {
      message("Cytoscape is connected and ready.")
  } else {
      stop("Could not connect to Cytoscape. Ensure the application is open and the CyREST port is accessible.")
  }
  if (is_connected) {
    print("Cytoscape connection established. Proceeding with the analysis.")
  } else {
    repeat {
      beep(10)
      answer <- readline(prompt = "Cytoscape is not open. Please open Cytoscape and confirm by typing 'y' to continue. Or, press 'c' to cancel the Cytoscape analysis and proceed without it: ")
      if (tolower(answer) == "c"){
        print("Cytoscape analysis will be skipped. Proceeding without Cytoscape.")
        run_cytoscape <<- FALSE
        break
      } else if (tolower(answer) == "y") {
        #check if cytoscape is open
        if (RCy3::cytoscapePing()) {
          print("Cytoscape is open. Proceeding with the analysis.")
          cytoscape_open <<- TRUE
          break
        } else {
          print("Cytoscape is still not responding. Please ensure Cytoscape is open and try again.")
        }
      } else {
        print("Invalid input. Please type 'y' to confirm Cytoscape is open or 'c' to cancel the Cytoscape analysis.")
      }
    }
  }
} else {
  print("Bypassing Cytoscape interface.")
}

chart_data_list <- list()

# Collection for TF enrichment scores across all parent terms
# Each entry: ParentTerm, TF, enrichment_score, NES_MN, NES_RGC, n_targets
tf_enrichment_list <- list()

# ==============================================================================
# SECTION 3: LOOP THROUGH EACH BROAD CATEGORY & BUILD DEDICATED SUBNETWORKS
# ==============================================================================

# Initialize the progress bar based on the total number of categories
total_categories <- length(parent_gene_lists)
pb <- txtProgressBar(min = 0, max = total_categories, style = 3)
category_counter <- 0

for (category_name in names(parent_gene_lists)) {
  
  # Increment our counter and update the progress bar at the start of every iteration
  category_counter <- category_counter + 1
  setTxtProgressBar(pb, category_counter)
  
  safe_name <- gsub("[^[:alnum:]_]", "_", category_name)
  category_dir <- file.path(base_output_dir, safe_name)
  if (!dir.exists(category_dir)) dir.create(category_dir)
  
  category_genes <- parent_gene_lists[[category_name]]
  is_ensembl <- any(grepl("^ENS", category_genes))
  
  if (is_ensembl) {
    category_symbols <- exp_map %>%
      dplyr::filter(ensembl_id %in% category_genes) %>%
      pull(symbol)
  } else {
    category_symbols <- category_genes
  }
  
  sub_edges <- master_edges %>%
    dplyr::filter(`Target Gene` %in% category_symbols)
  
  active_regulators <- unique(sub_edges$`Regulator Gene`)
  active_targets <- unique(sub_edges$`Target Gene`)
  all_active_nodes <- unique(c(active_regulators, active_targets))
  
  sub_nodes <- master_nodes %>%
    dplyr::filter(id %in% all_active_nodes)

  if (nrow(sub_edges) > 0) {
    
    cyto_edges <- sub_edges %>% 
      dplyr::rename(source = `Regulator Gene`, target = `Target Gene`) %>%
      dplyr::mutate(source = as.character(source), target = as.character(target))
    
    cyto_nodes <- sub_nodes %>%
      dplyr::mutate(id = as.character(id)) %>%
      dplyr::select(id, everything()) %>%
      as.data.frame()
    
    # ---- NETWORK CENTRALITIES CALCULATIONS (Always happens using igraph) ----
    g_temp <- graph_from_data_frame(d = cyto_edges, directed = TRUE)
    
    in_deg_v   <- degree(g_temp, mode = "in")
    out_deg_v  <- degree(g_temp, mode = "out")
    bet_v      <- betweenness(g_temp, directed = TRUE, normalized = FALSE)
    close_v    <- closeness(g_temp, mode = "out") 
    
    metrics_df <- data.frame(
      name = names(in_deg_v),
      InDegree = as.numeric(in_deg_v),
      OutDegree = as.numeric(out_deg_v),
      BetweennessCentrality = as.numeric(bet_v),
      ClosenessCentrality = as.numeric(close_v),
      stringsAsFactors = FALSE
    )
    
    # ---- TF TARGET ENRICHMENT (GSEA-style using expression data) ----
    # Batch all TFs in this parent term into two fgsea calls (MN + RGC)
    # so the permutation null is computed once per dataset.
    # Split TFs: those with >= 2 targets get full GSEA; those with < 2 get a
    # fallback (mean signed significance).
    fgsea_pathways <- list()
    fallback_list <- list()
    for (tf in active_regulators) {
      tf_targets <- sub_edges %>%
        dplyr::filter(`Regulator Gene` == tf) %>%
        dplyr::pull(`Target Gene`) %>%
        unique()
      if (length(tf_targets) >= 2) {
        fgsea_pathways[[tf]] <- tf_targets
      } else {
        fallback_list[[tf]] <- tf_targets
      }
    }
    
    category_tf_list <- list()
    
    # -- Batch GSEA for TFs with >= 2 targets --
    if (length(fgsea_pathways) > 0) {
      mn_res_all <- fgsea::fgsea(
        pathways = fgsea_pathways,
        stats = mn_ranked,
        minSize = 1,
        maxSize = 500,
        scoreType = "std"
      )
      rgc_res_all <- fgsea::fgsea(
        pathways = fgsea_pathways,
        stats = rgc_ranked,
        minSize = 1,
        maxSize = 500,
        scoreType = "std"
      )
      for (tf in names(fgsea_pathways)) {
        n_targets <- length(fgsea_pathways[[tf]])
        nes_mn <- mn_res_all$NES[mn_res_all$pathway == tf]
        nes_rgc <- rgc_res_all$NES[rgc_res_all$pathway == tf]
        # fgsea returns NA if the pathway wasn't tested; treat as missing
        if (length(nes_mn) == 0) nes_mn <- NA
        if (length(nes_rgc) == 0) nes_rgc <- NA
        combined_nes <- mean(c(nes_mn, nes_rgc), na.rm = TRUE)
        if (is.na(combined_nes)) next
        category_tf_list[[tf]] <- data.frame(
          ParentTerm = safe_name,
          TF = tf,
          enrichment_score = combined_nes,
          NES_MN = nes_mn,
          NES_RGC = nes_rgc,
          n_targets = n_targets,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # -- Fallback for TFs with < 2 targets --
    for (tf in names(fallback_list)) {
      tf_targets <- fallback_list[[tf]]
      n_targets <- length(tf_targets)
      target_stats <- exp_map %>%
        dplyr::filter(symbol %in% tf_targets)
      mean_mn <- mean(
        -log10(target_stats$pvalue_MN) * sign(target_stats$log2FoldChange_MN),
        na.rm = TRUE
      )
      mean_rgc <- mean(
        -log10(target_stats$pvalue_RGC) * sign(target_stats$log2FoldChange_RGC),
        na.rm = TRUE
      )
      combined_nes <- mean(c(mean_mn, mean_rgc), na.rm = TRUE)
      if (is.na(combined_nes)) next
      category_tf_list[[tf]] <- data.frame(
        ParentTerm = safe_name,
        TF = tf,
        enrichment_score = combined_nes,
        NES_MN = mean_mn,
        NES_RGC = mean_rgc,
        n_targets = n_targets,
        stringsAsFactors = FALSE
      )
    }
    
    # Append this category's results to the global list
    tf_enrichment_list[[safe_name]] <- bind_rows(category_tf_list)
    
    # ---- OPTIONAL CYTOSCAPE INTERACTION BLOCK ----
    if (run_cytoscape) {
      
      # capture.output and suppressMessages swallow the "Applying default style..." text
      capture.output(
        suppressMessages(
          createNetworkFromDataFrames(
            nodes = cyto_nodes, 
            edges = cyto_edges, 
            title = safe_name, 
            collection = "Parent_Subnetworks",
            key.col.index = 1
          )
        )
      )
      
      tryCatch({
        capture.output(suppressMessages(setVisualStyle("GRN")))
      }, error = function(e) {
        # Using warning() keeps the message out of the standard console stdout
        warning(paste("Could not apply style to network:", safe_name, "-", e$message))
      })
      
      capture.output(
        suppressMessages(
          loadTableData(
            data = metrics_df,
            data.key.column = "name",
            table = "node",
            table.key.column = "name"
          )
        )
      )
    }
    analyzed_nodes <- cyto_nodes %>%
      dplyr::rename(name = id) %>%
      left_join(metrics_df, by = "name")
    
    edge_name <- paste0("edges_", safe_name, ".txt")
    node_name <- paste0("metadata_", safe_name, ".txt")
    
    write.table(sub_edges, file.path(category_dir, edge_name), 
                row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(analyzed_nodes, file.path(category_dir, node_name), 
                row.names = FALSE, sep = "\t", quote = FALSE)
    writeLines(sort(active_regulators), file.path(category_dir, "regulators.txt"))
    
    counts <- analyzed_nodes %>%
      group_by(`Regulatory function`) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(parentTerm = category_name)
    
    chart_data_list[[category_name]] <- counts
  }
}

# Close the progress bar cleanly after the loop completes execution
close(pb)


# ==============================================================================
# SECTION 4: EXPORT TF TARGET ENRICHMENT RESULTS
# ==============================================================================

# Combine all per-category enrichment tables into one and save to disk
# This is used downstream by 9.visualize_subnetworks.R for the bubble plot
if (length(tf_enrichment_list) > 0) {
  all_enrichment <- bind_rows(tf_enrichment_list)
  write.csv(all_enrichment, "tf_target_enrichment.csv", row.names = FALSE)
  message("Exported TF target enrichment scores to tf_target_enrichment.csv (",
          nrow(all_enrichment), " rows)")
} else {
  message("No TF enrichment data was collected. Skipping export.")
}


# ==============================================================================
# SECTION 5: SAVE CYTOSCAPE SESSION FILE (Conditional on Choice)
# ==============================================================================

if (run_cytoscape) {
  session_file_path <- file.path(getwd(), "Parent_Subnetworks_Workspace.cys")
  saveSession(filename = session_file_path)
  print(paste("Cytoscape session saved successfully to:", session_file_path))
} else {
  print("Skipped saving .cys session file because Cytoscape interface run was disabled.")
}


# ==============================================================================
# SECTION 6: GENERATE COMPOSITION PLOT
# ==============================================================================

plot_df <- bind_rows(chart_data_list)

term_order <- plot_df %>%
  group_by(parentTerm) %>%
  summarise(total = sum(count)) %>%
  arrange(total) %>%
  pull(parentTerm)

plot_df$parentTerm <- factor(plot_df$parentTerm, levels = term_order)

network_summary_plot <- ggplot(plot_df, aes(x = parentTerm, y = count, fill = `Regulatory function`)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Regulator" = "#00BFC4", "Regulated" = "#F8766D")) +
  coord_flip() +
  labs(
    title = "Network Composition Across GO Parent Terms",
    subtitle = "Counts of unique Regulators vs Target Genes",
    x = NULL,
    y = "Number of Unique Genes",
    fill = "Role in Network"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_blank(), 
    legend.position = "bottom"
  )

print(network_summary_plot)
