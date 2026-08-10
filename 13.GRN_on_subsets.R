source("iregulon_reformat.R")
library("readr")
library("ggplot2")
library("patchwork")
library("gprofiler2")
library("simplifyEnrichment")
library(png)

source("utils/helpers.R")

results_dir <- "filtered_genelists/iregulon_results"
gene_dir    <- "filtered_genelists"
out_dir     <- "filtered_genelists/iregulon_graphs"
dir.create(out_dir, showWarnings = FALSE)

for (tsv in list.files(results_dir, pattern = "\\.tsv$", full.names = TRUE)) {
  geneset   <- sub("^; Name[[:space:]]*", "", grep("^; Name", readLines(tsv, n = 1), value = TRUE))
  gene_file <- file.path(gene_dir, geneset)
  name      <- sub("\\.tsv$", "", basename(tsv))
  graph_dir <- file.path(out_dir, name)
  graph_dir <- sub("ireg_", "", graph_dir)  # Remove "ireg_" prefix from the directory name
  graph_dir <- sub("_genes", "", graph_dir)  
  dir.create(graph_dir, showWarnings = FALSE)

  res <- reformat_networks(tsv, gene_file)

  write_tsv(res$edges, file.path(graph_dir, "edges.tsv"))
  write_tsv(res$nodes, file.path(graph_dir, "nodes.tsv"))
  writeLines(sort(unique(as.character(res$tfs$id))), file.path(graph_dir, "regulators.txt"))
}

message("Done. Outputs written to ", out_dir)


run_analysis <- FALSE
while (!run_analysis) {
  user_input <- readline(prompt = "Run analysis on the iregulon graphs? (y/n): ")
  if (tolower(user_input) == "y") {
    run_analysis <- TRUE
  } else if (tolower(user_input) == "n") {
    message("Analysis skipped. You can run it later by executing the relevant script.")
    stop()
  } else {
    message("Invalid input. Please enter 'y' for yes or 'n' for no.")
  }
}

iregulon_results_dir <- "filtered_genelists/graph_output"

cat("SET GO:BP CLUSTERING THRESHOLD\n")
cat("(Enter a value between 0.00 and 1.00, or press Enter for 0.80)\n\n")
chosen_threshold <- ""
while (TRUE) {
  chosen_threshold <- readline(prompt = "Enter threshold (0.00-1.00): ")
  chosen_threshold <- trimws(chosen_threshold)

  if (chosen_threshold == "") {
    chosen_threshold <- "0.80"
    cat(paste0("No input provided. Using default threshold: ", chosen_threshold, "\n"))
    break
  }

  chosen_h <- as.numeric(chosen_threshold)
  if (!is.na(chosen_h) && chosen_h > 0.00 && chosen_h <= 1.00) {
    break
  } else {
    cat("Invalid input.\n")
  }
}
chosen_h <- as.numeric(chosen_threshold)
cat(sprintf("Using threshold h = %.2f for all GO:BP clustering.\n\n", chosen_h))

for(csv in list.files(iregulon_results_dir, pattern = "\\.csv$", full.names = TRUE)) {
  #get files
  name<- sub("\\_graph_output.csv$", "", basename(csv))
  regulators <- readLines(file.path("filtered_genelists/iregulon_graphs" , name, "regulators.txt"))
  cytoscape_result <- read_csv(csv)
  genelist <- readLines(file.path("filtered_genelists" , paste0(name,"_genes.txt")))
  print(colnames(cytoscape_result))
  cytoscape_result %>%
    # Fix the NA values before doing the final column selection
    dplyr::mutate(
      InDegree = tidyr::replace_na(InDegree, 0),
      OutDegree = tidyr::replace_na(OutDegree, 0)
    ) %>%
    dplyr::select(
      #columns to remove:
      -`shared name`,
      #regular ordering
      name,
      ensembl_id,
      `Regulatory function`,
      EAE_Combined_Score,
      Specificity,
      InDegree,
      OutDegree,
      
      starts_with("log2FoldChange"),
      starts_with("pvalue"),
      starts_with("padj"),
      starts_with("stat"),
      starts_with("lfcSE"),
      
      starts_with("baseMean"),
      
      starts_with("is_sig"),
      selected,
      
      everything()
    )
  #get avg in and out degree
  cytoscape_result <-cytoscape_result %>% 
    dplyr::mutate(
      avg_outdegree = OutDegree / sum(OutDegree),
      avg_indegree = InDegree / sum(InDegree)
    )
  
  assign(paste0(name, "_cytoscape_output"), cytoscape_result)

  tf_ranking <- cytoscape_result %>%
    dplyr::filter(`Regulatory function` == "Regulator") %>%
    dplyr::arrange(desc(avg_outdegree))

  p1 <- ggplot(tf_ranking, aes(x = reorder(name, -avg_outdegree), y = avg_outdegree))+
    geom_bar(stat = "identity", fill = "mediumpurple") +
    labs(title = paste0("Average OutDegree of Transcription Factors in ", name),
        subtitle = "EAE upregulated genes in motor neurons and RGC's",
        x = "Transcription Factor",
          y = "Average OutDegree",
        caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

  top_targets <- cytoscape_result %>%
    dplyr::slice_max(order_by = avg_indegree, n = 20)

  p2 <- ggplot(top_targets, aes(x = reorder(name, -avg_indegree), y = avg_indegree)) +
    geom_bar(stat = "identity", fill = "coral") +
    labs(title = paste0("Top 20 Regulated Proteins in ", name),
        x = "Regulated Protein",
        y = "Average In-Degree") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  



  p3 <- ggplot(tf_ranking, aes(x = log2FoldChange_MN, y = log2FoldChange_RGC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = OutDegree), color = "#6a3d9a", alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = name), size = 3) +
    labs(title = paste0("Transcription Factor Log2 Fold Change in ", name),
        x = "Log2 Fold Change (Motor Neurons)",
        y = "Log2 Fold Change (RGCs)",
        size = "Out-Degree (Targets)") +
    theme_bw()
  
  print(
    (p1 |p2)/p3
  )
  GO <- gost(
    query = genelist,
    organism = "mmusculus",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE,
    user_threshold = 0.05,
    correction_method = "g_SCS",
    domain_scope = "annotated",
    custom_bg = NULL,
    numeric_ns = "",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP")
  )$result 
  print(GO)
  print(colnames(GO))
  GO_selected <- GO %>%
    dplyr::rename(
      adjusted_p_value = p_value,
      intersections = intersection
    ) %>%
    dplyr::mutate(
      negative_log10_of_adjusted_p_value = -log10(adjusted_p_value)
    )%>%
    dplyr::filter(source == "GO:BP") %>%
      arrange(desc(negative_log10_of_adjusted_p_value))%>%
    mutate(
      term_name = reorder(term_name, -negative_log10_of_adjusted_p_value)
    )
    assign(paste0(name, "_GO"), GO)
  GO_t20 <- GO_selected %>%
    dplyr::slice_head(n = 20)
  print(colnames(GO_selected))
    
  print(ggplot(GO_t20, aes(x = term_name, y = negative_log10_of_adjusted_p_value)) +
  geom_point(aes(size = intersection_size), color = "#f8766d") +
  scale_size_continuous(range = c(3, 10), name = "Count") +
    labs(
    x = NULL,
    y = expression(-log[10] ~ "adjusted pvalue"),
    title = paste0("Top 20 Enriched GO:BP Terms", " in ", name),
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


  #simplify enrichment
  go_bp_data <- GO_selected %>%
    dplyr::select(term_id, term_name, negative_log10_of_adjusted_p_value)

  term_ids <- go_bp_data$term_id
  scores <- setNames(go_bp_data$negative_log10_of_adjusted_p_value, term_ids)

  cat("Calculating similarity matrix...\n")
  mat <- GO_similarity(term_ids, db = "org.Mm.eg.db")

  cat("Performing hierarchical clustering...\n")
  hc <- hclust(as.dist(1 - mat), method = "average")

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

  term_to_genes <- GO_selected %>%
    dplyr::select(term_id, intersections) %>%
    tidyr::separate_rows(intersections, sep = ",\\s*|\\s+") %>%
    dplyr::filter(!is.na(intersections), intersections != "") %>%
    dplyr::distinct()

  parent_summary <- reducedTerms %>%
    group_by(parentTerm) %>%
    summarise(n_terms = dplyr::n(), .groups = "drop") %>%
    arrange(desc(n_terms))

  cat(sprintf("\n✓ Using threshold h = %.2f which created %d clusters\n", chosen_h, nrow(parent_summary)))

  termsummaryplot <- (ggplot(parent_summary, aes(x = reorder(parentTerm, n_terms), y = n_terms)) +
    geom_col(fill = "#2196F3") +
    coord_flip() +
    labs(
      title = paste0(sprintf("GO:BP Clustering", chosen_h, nrow(parent_summary)),
            " in ", name),
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

  genesummaryplot <- (ggplot(final_gene_summary, aes(x = reorder(parentTerm, n_genes), y = n_genes)) +
    geom_col(fill = "#E91E63") +
    coord_flip() +
    labs(
      title = paste0(sprintf("GO:BP Clustering", chosen_h, sum(final_gene_summary$n_genes)), " in ", name),
      x = NULL,
      y = "Number of genes"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 9)))
    

  fullplot <- (
    (p1 / p3) | genesummaryplot
  ) + plot_layout(widths = c(2, 1))
  
  print(fullplot)
  assign(paste0(name,"_fullplot"), fullplot)
}

if(file.exists("filtered_genelists/analysis_plots.pdf")){
  file.remove("filtered_genelists/analysis_plots.pdf")
}

pdf("filtered_genelists/analysis_plots.pdf", width = 14, height = 10)

for (plot_name in ls(pattern = "_fullplot$")) {
  name = sub("_fullplot$", "", plot_name) %>% tolower()
  print(get(plot_name))

  png_path <- paste0("filtered_genelists/", name, ".png")
  
  if (file.exists(png_path)) {
    grid.newpage() 
    
    # 1. Read and draw the image
    img <- readPNG(png_path)
    grid.raster(img)
    
    # 2. Add a customized title text overlay
    grid.text(
      label = name,                  # The text to display (e.g., your gene name)
      x = 0.5,                       # Horizontal center of the page (0 = left, 1 = right)
      y = 0.95,                      # Vertical position near the top (0 = bottom, 1 = top)
      gp = gpar(
        fontsize = 24,               # Make text large and clear
        fontface = "bold",           # Bold font
        col = "black"                # Text colour
      )
    )
    
  }
}

dev.off()
