# Complete script dependency set, including packages omitted by install_packages.R.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
cran <- c("beepr", "cluster", "cowplot", "dplyr", "fs", "futile.logger",
          "ggplot2", "ggplotify", "ggpubr", "ggrepel", "ggsignif", "gprofiler2",
          "gridExtra", "gtools", "httr", "igraph", "Matrix", "msigdbr",
          "patchwork", "pheatmap", "png", "purrr", "qpdf", "R.utils",
          "RColorBrewer", "readr", "readxl", "reticulate", "scales", "Seurat",
          "stringr", "tibble", "tidyr", "tidyverse", "VennDiagram", "withr",
          "XML", "BiocManager", "remotes", "magick", "shiny", "DT",
          "treemap", "wordcloud")
bioc <- c("AnnotationDbi", "AUCell", "BiocParallel", "biomaRt",
          "clusterProfiler", "DESeq2", "edgeR", "EnsDb.Mmusculus.v79",
          "ExperimentHub", "fgsea", "GEOquery", "GO.db", "GSEABase", "GSVA",
          "JASPAR2024", "limma", "org.Mm.eg.db", "RCy3", "rrvgo",
          "simplifyEnrichment", "singscore", "SummarizedExperiment", "TFBSTools", "UCell")
missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing(cran))) install.packages(missing(cran))
if (!requireNamespace("BiocManager", quietly = TRUE)) stop("BiocManager installation failed")
# Ubuntu's prebuilt BioC packages include a BiocVersion from R 4.4.
# Select the release supported by this R before resolving missing packages.
bioc_release <- switch(paste(R.version$major, strsplit(R.version$minor, "[.]")[[1]][1], sep = "."),
                       "4.5" = "3.22", "4.6" = "3.23",
                       stop("Review the Bioconductor release for this R version"))
BiocManager::install(version = bioc_release, update = FALSE, ask = FALSE)
# Align Ubuntu's older BioC dependency stack with the chosen release.
bioc_repositories <- BiocManager::repositories(version = bioc_release)
bioc_available <- available.packages(repos = bioc_repositories[grepl("bioconductor", bioc_repositories)])
installed <- installed.packages()
common <- intersect(rownames(installed), rownames(bioc_available))
align <- common[vapply(common, function(p) {
  package_version(installed[p, "Version"]) != package_version(bioc_available[p, "Version"])
}, logical(1))]
if (length(align)) BiocManager::install(align, force = TRUE, update = FALSE, ask = FALSE, Ncpus = 2)
if (length(missing(bioc))) BiocManager::install(missing(bioc), update = FALSE, ask = FALSE, Ncpus = 2)
if (!requireNamespace("nichenetr", quietly = TRUE)) {
  remotes::install_github("saeyslab/nichenetr", upgrade = "never", dependencies = NA)
}
failed <- missing(c(cran, bioc, "nichenetr"))
if (length(failed)) stop("Dependencies still missing: ", paste(failed, collapse = ", "))
dir.create("setup/logs", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")]),
          "setup/logs/r-installed-packages.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "setup/logs/r-session-info.txt")
message("All declared R dependencies load successfully.")
