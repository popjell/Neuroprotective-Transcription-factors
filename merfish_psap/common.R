# Shared setup for the merfish_psap scripts. Run them from merfish_psap/ with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
# Environment variables: RUN (names the results/ and figures/ subfolders, default gse284005)
library(ggplot2)
library(Matrix)
library(data.table)

set.seed(1010)
run <- Sys.getenv("RUN", "gse284005")
res_dir <- file.path("results", run)
fig_dir <- file.path("figures", run)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res <- function(file) file.path(res_dir, file)

min_gm_neurons <- 50      # neurons per donor in grey matter for a donor to count
neuron_label <- "Neuron"  # author majorCluster_final label
gm_label <- "GM"          # author Region_banksy label for grey matter

# Inputs from other analyses (read only)
rank_genes_csv <- "../rank_genes.csv"  # at the project root, not scRNAseq/_keep_input/
trobisch_donors <- c("CO39", "CO14", "CO36", "CO25", "MS60", "MS229", "MS371")

# Marker sets used to describe the MERFISH panel
markers <- list(
  Neuron_pan = c("RBFOX3", "SNAP25", "SYT1", "STMN2", "NEFL"),
  Excitatory = c("SLC17A7", "SLC17A6", "SATB2"),
  Inhibitory = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP", "RELN", "LAMP5"),
  Layer = c("CUX2", "RORB", "THEMIS", "FEZF2", "TLE4", "NR4A2", "BCL11B"),
  Microglia = c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "C1QA", "AIF1", "CD68", "TREM2"),
  Lipid_microglia = c("GPNMB", "PLIN2", "ABCA1", "SPP1", "APOE", "LPL"),
  Astrocyte = c("GFAP", "AQP4", "ALDH1L1", "GJA1", "SLC1A3"),
  Oligodendrocyte = c("MBP", "PLP1", "MOBP", "MOG", "MAG"))

theme_set(theme_bw(base_size = 12))
orange <- "#D9702E"
blue <- "#2F78A8"
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
}
no_overwrite <- function(file) if (file.exists(file)) stop(file, " exists; move it aside first")
