library(nichenetr)
library(tidyverse)

# Load NicheNet data (ligand_target_matrix, lr_network, weighted_networks)
source("nichenet-analysis/get_data.R")

# Drop all-zero ligand columns (e.g. Igfl3) which produce degenerate AUPR scores
all_zero <- colnames(ligand_target_matrix)[apply(ligand_target_matrix, 2, function(x) all(x == 0))]
if (length(all_zero) > 0) {
  message("Excluding ligands with no target scores (degenerate): ", paste(all_zero, collapse = ", "))
  ligand_target_matrix <- ligand_target_matrix[, !colnames(ligand_target_matrix) %in% all_zero]
}

# Read 4 gene lists, dedupe, filter to matrix genes
files <- c("GO_cell_death", "GO_metabolism", "LLM_cell_death", "LLM_metabolism")
genelists <- list()
for (f in files) {
  genes <- unique(readLines(file.path("filtered_genelists", paste0(f, "_genes.txt"))))
  genelists[[f]] <- genes[genes %in% rownames(ligand_target_matrix)]
}

# Load expression data for background gene definitions
exp_map <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Get target genes that exist in the ligand-target matrix
target_genes <- rownames(ligand_target_matrix)

# Prompt user for background gene set
cat("\n=== Background Gene Set Options ===\n")
cat("1. All genes\n")
cat("2. All detected genes (baseMean > 0 in both MN and RGC)\n")
cat("3. All upregulated genes (significant & positive FC in both MN and RGC)\n\n")
background_choice <- as.integer(readline(prompt = "Enter choice (1-3): "))
if (is.na(background_choice)) {
  background_choice <- as.integer(readLines(file("stdin"), n = 1))
}

# Define background genes based on user choice
if (background_choice == 1) {
  background_genes <- unique(exp_map$symbol)
  bg_label <- "all_genes"

} else if (background_choice == 2) {
  mn_detected <- exp_map$symbol[exp_map$baseMean_MN > 0]
  rgc_detected <- exp_map$symbol[exp_map$baseMean_RGC > 0]
  background_genes <- unique(intersect(mn_detected, rgc_detected))
  bg_label <- "detected"

} else if (background_choice == 3) {
  mn_up <- exp_map$symbol[exp_map$is_sig_MN == TRUE & exp_map$log2FoldChange_MN > 0]
  rgc_up <- exp_map$symbol[exp_map$is_sig_RGC == TRUE & exp_map$log2FoldChange_RGC > 0]
  background_genes <- unique(intersect(mn_up, rgc_up))
  bg_label <- "upregulated"

} else {
  stop("Invalid choice. Please enter 1, 2, or 3.")
}

# Filter background to genes in the ligand-target matrix and remove NAs
background_genes <- unique(background_genes[background_genes %in% target_genes])
background_genes <- background_genes[!is.na(background_genes)]

if (length(background_genes) == 0) {
  stop("Background gene set is empty. Try a different option.")
}

message("Background genes: ", length(background_genes))
for (name in names(genelists)) {
  message("  ", name, " geneset size: ", length(genelists[[name]]))
}

# Get all ligands from the ligand-target matrix
all_ligands <- colnames(ligand_target_matrix)

# Create output directory
output_dir <- "nichenet-analysis/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Run ligand activity inference for each gene list
all_results <- list()

for (name in names(genelists)) {
  message("\nAnalyzing: ", name)

  geneset <- genelists[[name]]

  activities <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = all_ligands
  )

  all_results[[name]] <- activities

  write.csv(activities,
            file.path(output_dir, paste0(name, "_", bg_label, "_ligand_activity.csv")),
            row.names = FALSE)
}

# Combine all results into one dataframe
combined <- bind_rows(all_results, .id = "gene_list")

write.csv(combined,
          file.path(output_dir, paste0("all_ligand_activities_", bg_label, ".csv")),
          row.names = FALSE)

# ============================================================
# BAR PLOT: top 50 ligands per gene list, ranked by AUPR
# ============================================================

for (name in names(genelists)) {
  ranked <- combined %>%
    filter(gene_list == name) %>%
    arrange(desc(aupr)) %>%
    head(50)

  p <- ggplot(ranked, aes(x = reorder(test_ligand, aupr), y = aupr, fill = aupr)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Top 50 Ligands —", name, "(", bg_label, "background )"),
         x = "Ligand", y = "AUPR") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal()

  ggsave(file.path(output_dir, paste0("barplot_top50_", name, "_", bg_label, ".pdf")),
         p, width = 10, height = 12)
}

message("\nDone! Results saved to: ", output_dir)
