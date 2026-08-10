
library(rrvgo)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(beepr)
library(gridExtra)

#GEO Plot
#region

#Run gprofiler analysis in web

#Import the .csv from gprofiler into R
go_analysis <- read.csv("GO_results.csv")
source("utils/helpers.R")
all_BP <- read_GO(go_analysis)
#export as csv
all_bp_terms <- all_BP %>%
  dplyr::select(-intersections)
#write.csv(all_bp_terms, "all_BP_terms.csv")

gprofiler_selected <- go_analysis %>%
dplyr::filter(source == "GO:BP") %>%
arrange(desc(negative_log10_of_adjusted_p_value))%>%
mutate(
term_name = reorder(term_name, -negative_log10_of_adjusted_p_value)
) 

gprofiler_selected_t20 <- gprofiler_selected %>%
  dplyr::slice_head(n = 20)

# Plot top 20 enriched terms
print(ggplot(gprofiler_selected_t20, aes(x = term_name, y = negative_log10_of_adjusted_p_value)) +
  geom_point(aes(size = intersection_size), color = "#f8766d") +
  scale_size_continuous(range = c(3, 10), name = "Count") +
    labs(
    x = NULL,
    y = expression(-log[10] ~ "adjusted pvalue"),
    title = "Top 20 Enriched GO:BP Terms",
    subtitle = "EAE upregulated genes in motor neurons and RGC's",
    caption = paste0("produced on ", Sys.time())
  ) +
  coord_flip() + 
  theme_bw() +
  theme(
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ))



# ==============================================================================
# USER CHOICE: ReviGO vs simplifyEnrichment
# ==============================================================================

if (!exists("method") || !(method %in% c("1", "2"))) {
  beep(10)
  cat("\n=== Select GO term reduction method ===\n\n")
  cat("1. ReviGO - use with scripts 8-9 (parent-category approach)\n")
  cat("2. simplifyEnrichment - use with scripts 8-9 (parent-category approach)\n")

  method <- ""
  while (!(method %in% c("1", "2"))) {
    method <- readline(prompt = "Select method (1 or 2): ")
    method <- trimws(method)
    if (!(method %in% c("1", "2"))) {
      cat("Invalid input. Please enter 1 or 2.\n")
    }
  }
} else {
  cat(sprintf("Using pre-set method: %s\n", ifelse(method == "1", "ReviGO", "simplifyEnrichment")))
}

# ==============================================================================
# METHOD 1: ReviGO (Original)
# ==============================================================================

if (method == "1") {
  
  cat("\n✓ Running ReviGO...\n\n")
  
  simMatrix <- calculateSimMatrix(go_analysis$term_id,
                                  orgdb = "org.Mm.eg.db",
                                  ont = "BP",
                                  method = "Rel")
  
  scores <- setNames(go_analysis$negative_log10_of_adjusted_p_value, 
                      go_analysis$term_id)
  
  if (!exists("rrvgo_threshold")) rrvgo_threshold <- 0.7
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = rrvgo_threshold,
                                  orgdb = "org.Mm.eg.db")
  
  # Create parent_summary table
  parent_summary <- reducedTerms %>%
    group_by(parentTerm) %>%
    summarise(
      n_terms    = n(),
      mean_score = mean(score),
      max_score  = max(score),
      top_term   = term[which.max(score)],
      .groups = "drop"
    ) %>%
    arrange(desc(n_terms))
  
  print(parent_summary, n = nrow(parent_summary))
  
  # Generate ReviGO visualizations
  print(ggplot(parent_summary, aes(x = reorder(parentTerm, n_terms), y = n_terms)) +
    geom_col(fill = "#f8766d") +
    coord_flip() +
    labs(title = "Number of GO:BP terms per parent category", x = NULL, y = "Number of terms") +
    theme_bw(base_size = 10))
  
  print(scatterPlot(simMatrix, reducedTerms))
  print(treemapPlot(reducedTerms))
  print(heatmapPlot(simMatrix, reducedTerms, annotateParent = TRUE,
                    annotationLabel = "parentTerm", fontsize = 6))
  
  cat("\n✓ ReviGO analysis complete.\n")
  
}

# ==============================================================================
# METHOD 2: simplifyEnrichment (Interactive Threshold Selection)
# ==============================================================================

if (method == "2") {
  
  cat("\n✓ Running simplifyEnrichment with threshold exploration...\n\n")
  
  library(simplifyEnrichment)
  
  # Prepare data
  go_bp_data <- go_analysis %>%
    dplyr::filter(source == "GO:BP") %>%
    dplyr::select(term_id, term_name, negative_log10_of_adjusted_p_value)
  
  term_ids <- go_bp_data$term_id
  scores <- setNames(go_bp_data$negative_log10_of_adjusted_p_value, term_ids)
  
  cat("Calculating similarity matrix...\n")
  mat <- GO_similarity(term_ids, db = "org.Mm.eg.db")
  
  cat("Performing hierarchical clustering...\n")
  hc <- hclust(as.dist(1 - mat), method = "average")
  
  # Display dendrogram
  cat("Generating dendrogram...\n\n")
  
  thresholds_to_test <- seq(0.70, 0.95, by = 0.05)

  if (!exists("chosen_threshold")) {

    plot(hc, 
         main = "GO:BP Term Hierarchical Clustering Dendrogram",
         xlab = "GO Terms", 
         ylab = "1 - Semantic Similarity",
         cex.main = 1.3,
         cex.lab = 1.1)
    
    colors <- c("red", "orange", "gold", "green", "blue", "purple")
    for (i in seq_along(thresholds_to_test)) {
      abline(h = thresholds_to_test[i], col = colors[i], lty = 2, lwd = 1.5)
    }
    
    legend("topright", 
           legend = sprintf("h = %.2f", thresholds_to_test),
           col = colors,
           lty = 2,
           lwd = 1.5,
           cex = 0.9)
    
    # Test all thresholds
    cat("\nTesting thresholds from 0.70 to 0.95...\n")
    
    term_to_genes <- go_analysis %>%
      dplyr::filter(source == "GO:BP") %>%
      dplyr::select(term_id, intersections) %>%
      tidyr::separate_rows(intersections, sep = ",\\s*|\\s+") %>%
      dplyr::filter(!is.na(intersections), intersections != "") %>%
      dplyr::distinct()
    
    threshold_plots <- list()
    gene_threshold_plots <- list()
    threshold_summaries <- list()
    
    for (test_h in thresholds_to_test) {
      test_clusters <- cutree(hc, h = test_h)
      n_clust <- length(unique(test_clusters))
      cat(sprintf("  h = %.2f  ->  %d clusters\n", test_h, n_clust))
      
      test_mapping <- data.frame(
        go = term_ids,
        cluster_id = test_clusters,
        term_name = go_bp_data$term_name,
        score = scores,
        stringsAsFactors = FALSE
      )
      
      test_reps <- test_mapping %>%
        group_by(cluster_id) %>%
        slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
        dplyr::select(cluster_id, parentTerm = term_name) %>%
        ungroup()
      
      test_reduced <- test_mapping %>%
        left_join(test_reps, by = "cluster_id") %>%
        dplyr::select(go, parentTerm, score)
      
      test_summary <- test_reduced %>%
        group_by(parentTerm) %>%
        summarise(
          n_terms    = n(),
          mean_score = mean(score),
          max_score  = max(score),
          .groups = "drop"
        ) %>%
        mutate(top_term = parentTerm) %>%
        dplyr::select(parentTerm, n_terms, mean_score, max_score, top_term) %>%
        arrange(desc(n_terms))
      
      threshold_summaries[[as.character(test_h)]] <- test_summary
      
      gene_summary <- test_reduced %>%
        dplyr::left_join(term_to_genes, by = c("go" = "term_id")) %>%
        dplyr::select(parentTerm, gene = intersections) %>%
        dplyr::filter(!is.na(gene)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(parentTerm) %>%
        dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(desc(n_genes))
      
      gp <- ggplot(gene_summary, aes(x = reorder(parentTerm, n_genes), y = n_genes)) +
        geom_col(fill = "#E91E63") +
        coord_flip() +
        labs(
          title = sprintf("h = %.2f (%d genes)", test_h, sum(gene_summary$n_genes)),
          x = NULL,
          y = "Number of genes"
        ) +
        theme_bw(base_size = 9) +
        theme(
          axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.margin = margin(5, 5, 5, 5)
        )
      
      gene_threshold_plots[[as.character(test_h)]] <- gp
      
      p <- ggplot(test_summary, aes(x = reorder(parentTerm, n_terms), y = n_terms)) +
        geom_col(fill = "#2196F3") +
        coord_flip() +
        labs(
          title = sprintf("h = %.2f (%d clusters)", test_h, nrow(test_summary)),
          x = NULL, 
          y = "Number of terms"
        ) +
        theme_bw(base_size = 9) +
        theme(
          axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.margin = margin(5, 5, 5, 5)
        )
      
      threshold_plots[[as.character(test_h)]] <- p
    }
    
    # Display comparison grids
    cat("\nGenerating comparison grid (term counts)...\n\n")
    do.call(gridExtra::grid.arrange, c(threshold_plots, ncol = 3, nrow = 2))
    
    cat("\nGenerating comparison grid (gene counts)...\n\n")
    do.call(gridExtra::grid.arrange, c(gene_threshold_plots, ncol = 3, nrow = 2))
    
    # User selection
    cat("SELECT YOUR PREFERRED THRESHOLD\n")
    
    cat("Which threshold do you prefer?\n")
    cat("(Enter a value between 0.70 and 0.95, or press Enter for 0.85)\n\n")
    
    chosen_threshold <- ""
    default_threshold <- "0.80"
    while (TRUE) {
      chosen_threshold <- readline(prompt = "Enter threshold (0.70-0.95): ")
      chosen_threshold <- trimws(chosen_threshold)
      
      if (chosen_threshold == "") {
        chosen_threshold <- default_threshold
        cat(paste0("No input provided. Using default threshold: ", default_threshold, "\n"))
        break
      }
      
      chosen_h <- as.numeric(chosen_threshold)
      if (!is.na(chosen_h) && chosen_h > 0.00 && chosen_h <= 1.00) {
        break
      } else {
        cat("Invalid input.\n")
      }
    }

  } else {
    cat(sprintf("Using pre-set threshold: %s\n", chosen_threshold))
  }
  
  chosen_h <- as.numeric(chosen_threshold)
  
  # Use chosen threshold for final output
  cluster_ids <- cutree(hc, h = chosen_h)
  
  cluster_mapping <- data.frame(
    go = term_ids,
    cluster_id = cluster_ids,
    term_name = go_bp_data$term_name,
    score = scores,
    stringsAsFactors = FALSE
  )
  
  cluster_reps <- cluster_mapping %>%
    group_by(cluster_id) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    dplyr::select(cluster_id, parentTerm = term_name) %>%
    ungroup()
  
  reducedTerms <- cluster_mapping %>%
    left_join(cluster_reps, by = "cluster_id") %>%
    dplyr::select(go, parentTerm, score)
  
  parent_summary <- reducedTerms %>%
    group_by(parentTerm) %>%
    summarise(
      n_terms    = n(),
      mean_score = mean(score),
      max_score  = max(score),
      .groups = "drop"
    ) %>%
    mutate(top_term = parentTerm) %>%
    dplyr::select(parentTerm, n_terms, mean_score, max_score, top_term) %>%
    arrange(desc(n_terms))
  
  cat("\n✓ Using threshold h =", chosen_h, "which created", nrow(parent_summary), "clusters\n\n")
  print(parent_summary, n = nrow(parent_summary))
  
  # Final visualization
  cat("\nGenerating final visualization...\n\n")
  print(ggplot(parent_summary, aes(x = reorder(parentTerm, n_terms), y = n_terms)) +
    geom_col(fill = "#2196F3") +
    coord_flip() +
    labs(
      title = sprintf("GO:BP Clustering (h = %.2f, %d clusters)", chosen_h, nrow(parent_summary)),
      x = NULL, 
      y = "Number of terms"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 9)))
  
  final_gene_summary <- reducedTerms %>%
    dplyr::left_join(term_to_genes, by = c("go" = "term_id")) %>%
    dplyr::select(parentTerm, gene = intersections) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(parentTerm) %>%
    dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(desc(n_genes))
  
  print(ggplot(final_gene_summary, aes(x = reorder(parentTerm, n_genes), y = n_genes)) +
    geom_col(fill = "#E91E63") +
    coord_flip() +
    labs(
      title = sprintf("GO:BP Clustering (h = %.2f, %d genes)", chosen_h, sum(final_gene_summary$n_genes)),
      x = NULL,
      y = "Number of genes"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 9)))
  
  cat("\n✓ simplifyEnrichment analysis complete.\n")
  
}





