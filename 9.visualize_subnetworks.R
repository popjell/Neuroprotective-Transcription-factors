# Script 9: builds TF summary plots and compiles a ranked Cytoscape network report PDF

library(dplyr)
library(tidyr)
library(fs)
library(ggplot2)
library(purrr)
library(qpdf)
library(cowplot)   
library(gridExtra) 
library(stringr)   
library(RCy3)
library(beepr)      
print("note: this requires the cytoscape analysis from script 8 to be completed first, and the user to have applied the correct layout in cytoscape before proceeding")
if (!exists("skip_prompts")) {
  cancel <- readline(prompt = "If you have not ran script 8 and applied the layout in cytoscape, type 'c' to cancel this script. Otherwise, press Enter to continue: ")
  if(tolower(cancel) == "c"){
    stop("User cancelled the script. Please run script 8 and apply the layout in Cytoscape before proceeding.")
  }
}
ping <- tryCatch({
  RCy3::cytoscapePing()
  TRUE
}, error = function(e) {
  FALSE
})
# Simple flag to prevent rerunning the layout check inside the same session
if(!exists("layout_applied")){
  layout_applied <<- FALSE
}

# Force user confirmation that Cytoscape layout is ready before exporting images
if(!layout_applied && !exists("skip_prompts")){
  repeat {
    beep(10)
    answer <- readline(prompt = "Have you applied the correct layout in cytoscape? (y/n): ")
    
    if (tolower(answer) == "y") {
      print("Proceeding.")
      layout_applied <<- TRUE
      break
    } else {
      print("Please apply the correct layout in Cytoscape before proceeding. The script will pause until you confirm.")
    }
  }
}

# Basic sanity checks so the pipeline fails early if inputs are missing
if (!dir.exists("parent_subnetworks")) {
  stop("Directory 'parent_subnetworks' not found. Please run the main generation script first.")
}

if (!file.exists("GRN_legend.gif") && !file.exists("GRN_legend.png")) {
  stop("Legend file ('GRN_legend.gif' or 'GRN_legend.png') not found in the root directory.")
}

# Pick whichever legend file exists
legend_file <- if (file.exists("GRN_legend.png")) "GRN_legend.png" else "GRN_legend.gif"


# Read all exported metadata tables and combine into one dataset
print("Reading exported local metadata sheets...")

metadata_files <- fs::dir_ls("parent_subnetworks", recurse = TRUE, glob = "*metadata_*.txt")

if (length(metadata_files) == 0) {
  stop("No metadata text files found inside 'parent_subnetworks'.")
}

# Merge all node tables and tag each row with its parent network name
all_network_nodes <- metadata_files %>%
  map_df(~ read.csv(.x, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>%
           mutate(ParentTerm = basename(dirname(.x))))

# Handle inconsistent column naming for out-degree values
outdegree_col <- if ("Outdegree" %in% as.character(colnames(all_network_nodes))) "Outdegree" else "OutDegree"

# Count total target genes (Regulated nodes) per parent term for percentage normalization
n_targets_per_parent <- all_network_nodes %>%
  dplyr::filter(`Regulatory function` == "Regulated") %>%
  dplyr::count(ParentTerm) %>%
  dplyr::rename(n_targets = n)

# Keep only regulators and compute % outdegree per subnetwork
top_regulators <- all_network_nodes %>%
  dplyr::filter(`Regulatory function` == "Regulator") %>%
  dplyr::left_join(n_targets_per_parent, by = "ParentTerm") %>%
  dplyr::mutate(pct_outdegree = (!!sym(outdegree_col)) / n_targets * 100) %>%
  dplyr::filter(pct_outdegree > 0)

# Read TF target enrichment scores computed by script 8
# Contains: ParentTerm, TF, enrichment_score (averaged NES across MN and RGC)
enrichment_file <- "tf_target_enrichment.csv"
if (file.exists(enrichment_file)) {
  tf_enrichment <- read.csv(enrichment_file, stringsAsFactors = FALSE)
  message("Loaded TF target enrichment scores from ", enrichment_file)
} else {
  message("WARNING: ", enrichment_file, " not found. Bubble plot will use % Out-Degree for color instead.")
  tf_enrichment <- NULL
}

# Merge enrichment scores into the regulator table (matching ParentTerm directory name with TF name)
if (!is.null(tf_enrichment)) {
  top_regulators <- top_regulators %>%
    dplyr::left_join(
      tf_enrichment %>% dplyr::select(ParentTerm, TF, enrichment_score, n_targets),
      by = c("ParentTerm" = "ParentTerm", "name" = "TF")
    )
} else {
  top_regulators <- top_regulators %>%
    dplyr::mutate(enrichment_score = pct_outdegree, n_targets = NA)
}

# ============== AXIS ORDERING CONFIG ==============
# Change these variables to control how axes are ordered.
# Set to NULL to keep alphabetical/default order.
# Options: "enrichment_score", "n_targets", "pct_outdegree", etc.
Y_AXIS_ORDER_BY <- "pct_outdegree"
X_AXIS_ORDER_BY <- "pct_outdegree"
# =================================================

if (!is.null(Y_AXIS_ORDER_BY) && Y_AXIS_ORDER_BY %in% colnames(top_regulators)) {
  tf_order <- top_regulators %>%
    dplyr::mutate(name = as.character(name)) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(ordering_var = mean(.data[[Y_AXIS_ORDER_BY]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(desc(ordering_var)) %>%
    dplyr::pull(name)
  top_regulators <- top_regulators %>%
    dplyr::mutate(name = as.character(name),
                  name = factor(name, levels = tf_order))
}

if (!is.null(X_AXIS_ORDER_BY) && X_AXIS_ORDER_BY %in% colnames(top_regulators)) {
  parent_order <- top_regulators %>%
    dplyr::mutate(ParentTerm = as.character(ParentTerm)) %>%
    dplyr::group_by(ParentTerm) %>%
    dplyr::summarise(ordering_var = mean(.data[[X_AXIS_ORDER_BY]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(ordering_var) %>%
    dplyr::pull(ParentTerm)
  top_regulators <- top_regulators %>%
    dplyr::mutate(ParentTerm = as.character(ParentTerm),
                  ParentTerm = factor(ParentTerm, levels = parent_order))
}

# Global overview plot showing regulator strength across GO parent terms
# Size = % outdegree (normalized by total targets in each subnetwork)
# Color = GSEA enrichment score (averaged NES from MN and RGC expression rankings)
tf_pathway_plot <- ggplot(top_regulators, aes(
  y = name,
  x = ParentTerm,
  size = pct_outdegree,
  color = enrichment_score
)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", limits = c(0, 10), oob = scales::squish) +
  theme_bw(base_size = 11) +
  labs(
    title = "Top Regulators Across GO Parent Categories",
    subtitle = "Dot size = % of target genes regulated; Color = expression enrichment (avg NES)",
    y = "Transcription Factor",
    x = "GO Parent Term",
    size = "% Out-Degree",
    color = "Enrichment\nScore (NES)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", face = "bold"),
    panel.grid.minor = element_blank()
  )

print(tf_pathway_plot)


# Start PDF report generation
output_pdf_name <- "subnetwork_report.pdf"

# Remove old report if it exists
if (file.exists(output_pdf_name)) {
  file.remove(output_pdf_name)
}

# Check if Cytoscape is reachable before trying to export images
cytoscape_active <- FALSE

tryCatch({
  cytoscapePing()
  cytoscape_active <- TRUE
}, error = function(e) {
  print("Cytoscape is not actively connected. Skipping structural map image generation phase.")
})

# Only continue if Cytoscape is available
if (cytoscape_active && !exists("skip_cytoscape_pdf")) {

  # ---- MONKEY-PATCH: Fix RCy3 commandsGET bug (v2.32.0) ----
  # When Cytoscape returns a response without <p> tags, xpathSApply returns
  # an empty list and startsWith(res.elem[1], "[") crashes. Patch it to
  # fall back to raw character parsing.
  .patched_commandsGET <- function(command, base.url = RCy3:::.defaultBaseUrl) {
    q.url <- RCy3:::.command2getQuery(command, base.url)
    res <- tryCatch(
      httr::GET(q.url),
      error = function(e) stop(e)
    )
    res.html <- XML::htmlParse(httr::content(res, as = "text", encoding = "UTF-8"), asText = TRUE)
    res.elem <- XML::xpathSApply(res.html, "//p", XML::xmlValue)
    if (length(res.elem) == 0 || !is.character(res.elem[[1]])) {
      # Fallback: grab raw text content from the response body
      raw <- httr::content(res, as = "text", encoding = "UTF-8")
      raw <- gsub("<[^>]*>", "", raw)
      res.list <- unlist(strsplit(raw, "\n\\s*"))
      res.list <- trimws(res.list)
      res.list <- res.list[!(res.list == "" | res.list == "Finished")]
      if (length(res.list) > 0) return(res.list) else return(invisible(res.list))
    }
    if (startsWith(res.elem[1], "[")) {
      res.elem[1] <- gsub("\\[|\\]|\"", "", res.elem[1])
      res.elem2 <- unlist(strsplit(res.elem[1], "\n"))[1]
      res.list <- unlist(strsplit(res.elem2, ","))
    } else {
      res.list <- unlist(strsplit(res.elem[1], "\n\\s*"))
      res.list <- res.list[!(res.list == "Finished")]
    }
    if (length(res.list) > 0) res.list else invisible(res.list)
  }
  assignInNamespace("commandsGET", .patched_commandsGET, ns = "RCy3")
  if (!identical(body(get("commandsGET", envir = asNamespace("RCy3"))), body(.patched_commandsGET))) {
    env <- as.environment("package:RCy3")
    if (bindingIsLocked("commandsGET", env)) unlockBinding("commandsGET", env)
    assign("commandsGET", .patched_commandsGET, envir = env)
  }
  # ---- END PATCH ----
  
  print("Starting programmatic layout compilation with multi-panel plots sorted by significance...")
  
  pdf_dir <- "temp_pdf_pages"
  if (!dir.exists(pdf_dir)) dir.create(pdf_dir)

  # Decide network order based on enrichment size if available
  if (exists("parent_summary") && "n_terms" %in% as.character(colnames(parent_summary))) {
    
    network_ranking <- parent_summary %>%
      dplyr::mutate(NetworkMatchName = gsub("[^[:alnum:]_]", "_", parentTerm)) %>%
      dplyr::select(NetworkMatchName, n_terms) %>%
      dplyr::arrange(desc(n_terms))
      
    ordered_networks <- network_ranking$NetworkMatchName

  } else {

    print("Warning: 'parent_summary' dataframe or 'n_terms' column not found in global env. Falling back to default list.")

    ordered_networks <- getNetworkList()

    network_ranking <- tibble(
      NetworkMatchName = ordered_networks,
      n_terms = NA
    )
  }
  
  pdf_compiled_list <- list()
  
  # First page: simple cover page with title and date
  title_main <- cowplot::ggdraw() + 
    cowplot::draw_label("Subregulatory networks", fontface = 'bold', size = 28, x = 0.1, y = 0.6, hjust = 0, lineheight = 1.2) +
    cowplot::draw_label(paste("Generated:", Sys.Date()), size = 12, x = 0.1, y = 0.4, hjust = 0, color = "gray30")
  
  title_page_path <- file.path(pdf_dir, "00_title_page.pdf")
  ggsave(title_page_path, plot = title_main, width = 8.5, height = 11)
  pdf_compiled_list <- c(pdf_compiled_list, title_page_path)

  # Second page: legend explaining visual encoding used in networks
  legend_graphic <- cowplot::ggdraw() + cowplot::draw_image(legend_file)

  legend_title <- cowplot::ggdraw() +
    cowplot::draw_label("Visual Style Mapping Key", fontface = 'bold', size = 18, x = 0.05, hjust = 0)

  legend_page <- cowplot::plot_grid(
    legend_title,
    legend_graphic,
    rel_heights = c(1, 9),
    ncol = 1
  )
  
  legend_page_path <- file.path(pdf_dir, "01_legend_page.pdf")
  ggsave(legend_page_path, plot = legend_page, width = 8.5, height = 11)
  pdf_compiled_list <- c(pdf_compiled_list, legend_page_path)
  
  # Get network list once before the loop
  cytoscape_networks <- tryCatch({
    as.character(getNetworkList())
  }, error = function(e) {
    cat("Warning: Could not retrieve network list from Cytoscape\n")
    NULL
  })
  
  cat("Compiling", length(ordered_networks), "network pages...\n")
  
  # Build one report page per network
  pages_generated <- 0
  for (rank_idx in seq_along(ordered_networks)) {

    net_name <- ordered_networks[rank_idx]
    
    # Skip networks that are not currently loaded in Cytoscape
    in_cytoscape <- is.null(cytoscape_networks) || (net_name %in% cytoscape_networks)
    if (!in_cytoscape) next
    
    # Set current network
    network_ok <- tryCatch({
      setCurrentNetwork(network = net_name)
      TRUE
    }, error = function(e) {
      warning("Could not set network '", net_name, "': ", conditionMessage(e), ". Skipping.")
      FALSE
    })
    if (!network_ok) next
    
    # Ensure clean selection state before exporting visuals
    clearSelection(network = net_name)

    # Force Cytoscape to render updated layout before screenshot
    fitContent()
    Sys.sleep(0.5)
    
    temp_png <- file.path(pdf_dir, paste0("temp_net_", rank_idx, ".png"))

    # Export current network view as image for embedding in PDF
    export_ok <- tryCatch({
      exportImage(
        filename = temp_png,
        type = "png",
        resolution = 300,
        overwriteFile = TRUE
      )
      TRUE
    }, error = function(e) {
      cat("  WARNING: exportImage failed for", net_name, ":", conditionMessage(e), "\n")
      FALSE
    })
    if (!export_ok) next
    
    # Verify the exported image exists and has non-zero size
    if (!file.exists(temp_png) || file.info(temp_png)$size == 0) {
      cat("  WARNING: Exported PNG is missing or empty for", net_name, "\n")
      cat("  Looking for alternatives:", paste(list.files(pdf_dir, pattern = paste0("temp_net_", rank_idx)), collapse = ", "), "\n")
      # Cytoscape may append extension or use different naming
      candidates <- list.files(pdf_dir, pattern = paste0("temp_net_", rank_idx), full.names = TRUE)
      if (length(candidates) > 0) {
        temp_png <- candidates[which.max(file.info(candidates)$size)]
        cat("  Using:", temp_png, "(", file.info(temp_png)$size, "bytes )\n")
      } else {
        cat("  SKIP: no valid image found\n")
        next
      }
    } else {
      cat("  Image OK:", temp_png, "(", file.info(temp_png)$size, "bytes )\n")
    }
  
    # Extract nodes belonging to this specific network
    local_nodes <- all_network_nodes %>%
      dplyr::filter(ParentTerm == net_name)
    
    # Top transcription factors ranked by outgoing connections
    local_tfs <- local_nodes %>%
      dplyr::filter(`Regulatory function` == "Regulator") %>%
      dplyr::arrange(desc(!!sym(outdegree_col))) %>%
      head(12)
    
    # If no regulators exist, show placeholder instead of empty plot
    if(nrow(local_tfs) == 0) {

      tf_bar_chart <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No local regulators found", size = 3) +
        theme_void()

    } else {

      tf_bar_chart <- ggplot(local_tfs, aes(
        x = reorder(name, !!sym(outdegree_col)),
        y = !!sym(outdegree_col)
      )) +
        geom_col(fill = "steelblue", width = 0.7) +
        coord_flip() +
        theme_minimal(base_size = 8) +
        labs(title = "Top TFs by Target Count (Out-Degree)", x = NULL, y = "Out-Degree") +
        theme(
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", color = "black")
        )
    }
    
    # Identify target genes and rank by incoming connections
    indegree_col <- if ("Indegree" %in% as.character(colnames(all_network_nodes))) "Indegree" else "InDegree"
    
    local_targets <- local_nodes %>%
      dplyr::filter(`Regulatory function` != "Regulator") %>%
      dplyr::arrange(desc(!!sym(indegree_col))) %>%
      head(12)
    
    # Placeholder if no target genes exist in this network
    if(nrow(local_targets) == 0) {

      target_bar_chart <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No target genes found", size = 3) +
        theme_void()

    } else {

      target_bar_chart <- ggplot(local_targets, aes(
        x = reorder(name, !!sym(indegree_col)),
        y = !!sym(indegree_col)
      )) +
        geom_col(fill = "darkorange3", width = 0.7) +
        coord_flip() +
        theme_minimal(base_size = 8) +
        labs(title = "Top Target Genes (In-Degree)", x = NULL, y = "In-Degree") +
        theme(
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", color = "black")
        )
    }

    # Combine network image and summary plots into a single report page
    network_img <- cowplot::ggdraw() +
      cowplot::draw_image(temp_png, scale = 1.0)

    page_header_text <- gsub("_", " ", net_name)

    title_theme <- cowplot::ggdraw() +
      cowplot::draw_label(page_header_text, fontface = 'bold', size = 20, x = 0.02, hjust = 0)
    
    # Bottom section shows TFs and target gene summaries side by side
    bottom_charts_row <- cowplot::plot_grid(
      target_bar_chart,
      tf_bar_chart,
      rel_widths = c(1, 1),
      nrow = 1
    )
    
    # Stack title, network image, and summary charts vertically
    full_page <- cowplot::plot_grid(
      title_theme,
      network_img,
      bottom_charts_row,
      rel_heights = c(0.5, 6.5, 3.0),
      ncol = 1
    )

    # Prefix filenames so sorting preserves intended page order
    page_prefix <- paste0(sprintf("%02d", rank_idx + 1), "_")

    page_path <- file.path(
      pdf_dir,
      paste0(page_prefix, "ranked_", rank_idx, ".pdf")
    )
    if(file.exists(page_path)) file.remove(page_path)
    
    ggsave(page_path, plot = full_page, width = 11, height = 8.5)
    
    pdf_compiled_list <- c(pdf_compiled_list, page_path)
    pages_generated <- pages_generated + 1
    cat("  -> Page", pages_generated, "saved:", basename(page_path), "\n")

    # Clean up temporary image after page is saved
    if (file.exists(temp_png)) file.remove(temp_png)
  }
  
  cat("--- Loop complete: generated", pages_generated, "network pages ---\n")
  cat("Total files in compiled list:", length(pdf_compiled_list), "\n")
  
  # Merge all pages into a single final report PDF
  sorted_pdf_files <- sort(unlist(pdf_compiled_list))
  
  if (length(sorted_pdf_files) > 0) {

    qpdf::pdf_combine(
      input = sorted_pdf_files,
      output = output_pdf_name
    )
    
    # Remove temp folder once everything is merged
    unlink(pdf_dir, recursive = TRUE)

    print(paste("made report: ", output_pdf_name))
  }
}