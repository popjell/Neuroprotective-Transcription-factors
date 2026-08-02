library(dplyr)
library(tidyr)
library(ggplot2)
library(fs)
library(stringr)
library(beepr)
library(org.Mm.eg.db)
library(rrvgo)
library(gridExtra)
library(RCy3)
library(igraph)
library(fgsea)

go_analysis <- read.csv("GO_results.csv", stringsAsFactors = FALSE)
intersected_genes <- readLines("intersected_upregulated_ensembl_ids.txt")
intersected_genes <- intersected_genes[intersected_genes != ""]
source("iregulon_reformat.R")

skip_prompts <- TRUE
layout_applied <- TRUE
run_cytoscape <- FALSE
skip_cytoscape_pdf <- TRUE

methods <- c("1", "2")
thresholds <- seq(0.70, 0.95, by = 0.05)
gene_modes <- c(FALSE, TRUE)
base_dirs <- c("subsetted options", "subsetted options_unique")
method_dirs <- c("rrvgo", "simplifyEnrichment")
method_names <- c("1" = "rrvgo", "2" = "simplifyEnrichment")

for (bd in base_dirs) {
  for (md in method_dirs) {
    for (th in thresholds) {
      dir.create(file.path(bd, md, sprintf("%.2f", th)), recursive = TRUE, showWarnings = FALSE)
    }
  }
}

total <- length(gene_modes) * length(methods) * length(thresholds)
n <- 0

cleanup_cytoscape <- function() {
  tryCatch({
    nets <- RCy3::getNetworkList()
    for (nm in nets) tryCatch(RCy3::deleteNetwork(nm), error = function(e) NULL)
  }, error = function(e) NULL)
}

save_outputs <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (exists("reducedTerms") && exists("parent_summary")) {
    save(reducedTerms, parent_summary, file = file.path(out_dir, "go_results.RData"))
  }
  if (dir.exists("parent_subnetworks")) {
    tgt <- file.path(out_dir, "parent_subnetworks")
    if (dir.exists(tgt)) unlink(tgt, recursive = TRUE)
    file.copy("parent_subnetworks", out_dir, recursive = TRUE)
  }
  if (dir.exists("parent_category_genes")) {
    tgt <- file.path(out_dir, "parent_category_genes")
    if (dir.exists(tgt)) unlink(tgt, recursive = TRUE)
    file.copy("parent_category_genes", out_dir, recursive = TRUE)
  }
  for (f in c("tf_target_enrichment.csv", "subnetwork_report.pdf")) {
    if (file.exists(f)) file.copy(f, out_dir, overwrite = TRUE)
  }
  if (exists("tf_pathway_plot")) {
    ggsave(file.path(out_dir, "dotplot.pdf"), plot = tf_pathway_plot, width = 14, height = 8, dpi = 300)
    ggsave(file.path(out_dir, "dotplot.png"), plot = tf_pathway_plot, width = 14, height = 8, dpi = 300)
  }
}

for (gene_mode in gene_modes) {
  base_dir <- if (gene_mode) "subsetted options_unique" else "subsetted options"
  use_unique <- gene_mode

  for (method in methods) {
    m_name <- method_names[method]

    for (threshold in thresholds) {
      n <- n + 1
      th_label <- sprintf("%.2f", threshold)
      cat(sprintf("\n=== %d/%d: %s/%s/%s ===\n", n, total, base_dir, m_name, th_label))

      if (method == "1") {
        rrvgo_threshold <- threshold
        if (exists("chosen_threshold")) rm(chosen_threshold)
      } else {
        chosen_threshold <- th_label
        if (exists("rrvgo_threshold")) rm(rrvgo_threshold)
      }

      cat("  script 6...\n")
      tryCatch(
        source("6.GO_analysis.R", local = FALSE, print.eval = FALSE),
        error = function(e) { cat("  ERROR in script 6:", conditionMessage(e), "\n") }
      )

      out_dir <- file.path(base_dir, m_name, th_label)

      cleanup_cytoscape()

      cat("  script 8...\n")
      tryCatch(
        source("8.generate_subregulatory_networks.R", local = FALSE, print.eval = FALSE),
        error = function(e) { cat("  ERROR in script 8:", conditionMessage(e), "\n") }
      )

      save_outputs(out_dir)

      cat("  script 9...\n")
      tryCatch(
        source("9.visualize_subnetworks.R", local = FALSE, print.eval = FALSE),
        error = function(e) { cat("  ERROR in script 9:", conditionMessage(e), "\n") }
      )

      save_outputs(out_dir)

      cat(sprintf("  saved to %s/\n", out_dir))
    }
  }
}

cat("\nDone! All", total, "combinations complete.\n")
beep(1)
