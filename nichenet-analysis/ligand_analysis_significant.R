library(nichenetr)
library(tidyverse)

# Load local NicheNet mouse data
ligand_target_matrix <- readRDS("nichenet_data/ligand_target_matrix_nsga2r_final_mouse.rds")

# Drop all-zero ligand columns (e.g. Igfl3) which produce degenerate AUPR scores
all_zero <- colnames(ligand_target_matrix)[apply(ligand_target_matrix, 2, function(x) all(x == 0))]
if (length(all_zero) > 0) {
  message("Excluding ligands with no target scores (degenerate): ", paste(all_zero, collapse = ", "))
  ligand_target_matrix <- ligand_target_matrix[, !colnames(ligand_target_matrix) %in% all_zero]
}

# Read 4 gene lists, dedupe, filter to matrix
files <- c("GO_cell_death", "GO_metabolism", "LLM_cell_death", "LLM_metabolism")
genelists <- list()
for (f in files) {
  genes <- unique(readLines(file.path("filtered_genelists", paste0(f, "_genes.txt"))))
  genelists[[f]] <- genes[genes %in% rownames(ligand_target_matrix)]
}

# Load expression data for background
exp_map <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Prompt user for background gene set
cat("\n=== Background Gene Set Options ===\n")
cat("1. All genes\n")
cat("2. All detected genes (baseMean > 0 in both MN and RGC)\n")
cat("3. All upregulated genes (significant & positive FC in both MN and RGC)\n\n")
background_choice <- as.integer(readline(prompt = "Enter choice (1-3): "))
if (is.na(background_choice)) {
  background_choice <- as.integer(readLines(file("stdin"), n = 1))
}

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

# Filter background to matrix genes
background_genes <- unique(background_genes[background_genes %in% rownames(ligand_target_matrix)])
background_genes <- background_genes[!is.na(background_genes)]

if (length(background_genes) == 0) {
  stop("Background gene set is empty. Try a different option.")
}

message("Background genes: ", length(background_genes))

all_ligands <- colnames(ligand_target_matrix)
output_dir <- "nichenet-analysis/results_significant"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_iter <- 100
n_cores <- 16

for (name in names(genelists)) {
  message("\nAnalyzing: ", name, " (geneset size ", length(genelists[[name]]), ")")

  geneset <- genelists[[name]]

  # Observed AUPR
  observed <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = all_ligands
  )

  # Permutation null distribution
  perms <- bootstrap_ligand_activity_analysis(
    expressed_genes_receiver = background_genes,
    geneset_oi = geneset,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = all_ligands,
    n_iter = n_iter,
    n_cores = n_cores
  )

  # Build null AUPR matrix per ligand
  zero <- rep(0, length(all_ligands))
  names(zero) <- all_ligands
  perm_aupr <- lapply(perms, function(res) {
    m <- zero
    m[res$test_ligand] <- res$aupr
    m
  })
  perm_matrix <- do.call(rbind, perm_aupr)

  # Empirical p-value and FDR
  observed$p_value <- observed$test_ligand %>%
    sapply(function(lig) mean(perm_matrix[, lig] >= observed$aupr[observed$test_ligand == lig]))
  observed$p_value[observed$p_value == 0] <- 1 / n_iter
  observed$fdr <- p.adjust(observed$p_value, method = "BH")

  write.csv(observed, file.path(output_dir, paste0(name, "_", bg_label, "_significant.csv")),
            row.names = FALSE)

  # Bar plot of significant ligands
  significant <- observed %>% filter(fdr < 0.05) %>% arrange(desc(aupr))
  message("  Significant ligands (FDR < 0.05): ", nrow(significant))

  p <- ggplot(significant, aes(x = reorder(test_ligand, aupr), y = aupr, fill = aupr)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Significant Ligands -", name, "(", bg_label, "background)"),
         x = "Ligand", y = "AUPR") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal()

  ggsave(file.path(output_dir, paste0("barplot_significant_", name, "_", bg_label, ".pdf")),
         p, width = 10, height = max(6, nrow(significant) * 0.35))
}

message("\nDone! Results saved to: ", output_dir)
