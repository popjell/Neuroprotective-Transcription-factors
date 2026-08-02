################################################################################
# Installation Script for Neuroprotective Transcription Factors Repository
# Date: 2026-06-12
# Purpose: Install all required R packages for the analysis pipeline
# 
# Usage: 
#   source("install_packages.R")
#
# This script will:
# 1. Install all CRAN packages needed for the analysis
# 2. Install all Bioconductor packages with automatic dependency resolution
# 3. Report installation success/failure for each package
################################################################################

message("\n", strrep("=", 80))
message("Starting package installation for Neuroprotective Transcription Factors")
message("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message(strrep("=", 80), "\n")

# Helper function to safely install packages
install_if_needed <- function(package, source = "CRAN") {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("  Installing ", package, "...")
    tryCatch({
      if (source == "Bioconductor") {
        # Ensure BiocManager is available
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          message("    -> BiocManager not found. Installing...")
          install.packages("BiocManager", quiet = TRUE)
        }
        BiocManager::install(package, update = FALSE, ask = FALSE, quiet = TRUE)
      } else {
        install.packages(package, quiet = TRUE)
      }
      
      if (requireNamespace(package, quietly = TRUE)) {
        message("  ✓ ", package, " installed successfully")
        return(TRUE)
      } else {
        message("  ✗ ", package, " installation failed (package not available)")
        return(FALSE)
      }
    }, error = function(e) {
      message("  ✗ ", package, " installation failed: ", conditionMessage(e))
      return(FALSE)
    })
  } else {
    message("  ✓ ", package, " already installed")
    return(TRUE)
  }
}

# Track installation results
results <- list(
  cran_success = character(),
  cran_failed = character(),
  bioc_success = character(),
  bioc_failed = character()
)

################################################################################
# CRAN PACKAGES
################################################################################

message("\n", strrep("-", 80))
message("Installing CRAN packages (data manipulation, visualization, utilities)")
message(strrep("-", 80), "\n")

cran_packages <- c(
  # Data manipulation & core (used across multiple scripts)
  "dplyr",                    # Data manipulation (7 script uses)
  "tidyverse",                # Meta-package: ggplot2, dplyr, tidyr, readr (5 uses)
  "tidyr",                    # Data tidying
  "tibble",                   # Modern data frames
  "readr",                    # CSV file reading
  "fs",                       # File system operations (2 uses)
  "stringr",                   # String manipulation (2 uses)
  
  # Visualization packages
  "magick",                   # Image processing (used in 9.visualize_subnetworks.R)
  "pheatmap",                 # Heatmap visualization (4 uses)
  "ggplot2",                  # Plotting (implicit in tidyverse, but explicit for Motor Neuron/RGC analyses)
  "ggrepel",
  "VennDiagram",              # Venn diagram (5.gene_intersection.R)
  "grid",                     # Graphics utilities (VennDiagram dependency)
  "cowplot",                   # Plot arrangement 
  "gridExtra",                 # Additional grid utilities
  
  # Utilities
  "readxl",                   # Excel file reading (2 uses)
  "R.utils",                  # General utilities (2 uses)
  "futile.logger",             # Logging control (5.gene_intersection.R)
  "purrr",                     # Functional programming (used in 9.visualize_subnetworks.R) 
  "qpdf",                     # PDF manipulation (used in 9.visualize_subnetworks.R)
  "beepr",                    # Audio notifications (used in 0.Pipeline_execution.R)
  #shiny app stuff
  "shiny", "DT", "ggplot2", "treemap", "wordcloud"
)

for (pkg in cran_packages) {
  success <- install_if_needed(pkg, source = "CRAN")
  if (success) {
    results$cran_success <- c(results$cran_success, pkg)
  } else {
    results$cran_failed <- c(results$cran_failed, pkg)
  }
}
################################################################################
# BIOCONDUCTOR PACKAGES
################################################################################

message("\n", strrep("-", 80))
message("Installing Bioconductor packages (genomics & enrichment analysis)")
message(strrep("-", 80), "\n")

# Ensure BiocManager is installed first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager", quiet = TRUE)
}

bioc_packages <- c(
  "DESeq2", "GEOquery", "EnsDb.Mmusculus.v79", "biomaRt",
  "rrvgo", "org.Mm.eg.db", "GSEABase", "ExperimentHub", 
  "msigdbr", "GSVA", "RCy3", "simplifyEnrichment", "fgsea",
  "clusterProfiler"
)

# Identify which ones are actually missing
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  
  # Bulk install forces BiocManager to keep the repos active for all dependencies
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE, quiet = TRUE)
}

# Verify and populate your results tracking
for (pkg in bioc_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    results$bioc_success <- c(results$bioc_success, pkg)
  } else {
    results$bioc_failed <- c(results$bioc_failed, pkg)
  }
}
################################################################################
# SUMMARY REPORT
################################################################################

message("\n", strrep("=", 80))
message("INSTALLATION SUMMARY")
message(strrep("=", 80), "\n")

total_success <- length(results$cran_success) + length(results$bioc_success)
total_failed <- length(results$cran_failed) + length(results$bioc_failed)
total_packages <- length(cran_packages) + length(bioc_packages)

message("CRAN Packages:      ", length(results$cran_success), "/", length(cran_packages), " installed")
if (length(results$cran_failed) > 0) {
  message("  Failed: ", paste(results$cran_failed, collapse = ", "))
}

message("Bioconductor:       ", length(results$bioc_success), "/", length(bioc_packages), " installed")
if (length(results$bioc_failed) > 0) {
  message("  Failed: ", paste(results$bioc_failed, collapse = ", "))
}

message("\nTotal: ", total_success, "/", total_packages, " packages successfully installed")

if (total_failed == 0) {
  message("\n✓ All packages installed successfully!")
} else {
  message("\n⚠ Warning: ", total_failed, " package(s) failed to install.")
  message("  Please check internet connection and package availability.")
}

message(strrep("=", 80), "\n")
