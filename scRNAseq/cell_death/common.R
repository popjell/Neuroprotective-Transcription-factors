# Shared setup for the cell_death scripts. Run them from scRNAseq/cell_death/.
# GENE (default PSAP) is tested against every set in GENESETS (a .gmt file, default data/yu_tier2_primary.gmt).
# Outputs go to results/<GENE>_<set>/ and figures/<GENE>_<set>/. SET_LABEL names the sets in the figure title.
# data/exclude_genes.txt: the 31 Rank 2 genes and every gene of the earlier gene sets, kept out of the matched-gene pool.
library(ggplot2)
library(Matrix)

set.seed(1010)
gene <- Sys.getenv("GENE", "PSAP")
genesets <- Sys.getenv("GENESETS", "data/yu_tier2_primary.gmt")
set_label <- Sys.getenv("SET_LABEL", "cell-death program")
stopifnot(file.exists(genesets))
run <- paste(gene, tools::file_path_sans_ext(basename(genesets)), sep = "_")
res_dir <- file.path("results", run)
fig_dir <- file.path("figures", run)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res <- function(file) file.path(res_dir, file)

min_cells <- 20      # cells per donor x subtype
min_donors <- 2      # control donors and MS donors per subtype
detect_frac <- 0.05  # a set is tested in a subtype if at least detect_genes of its genes are detected in this fraction of cells
detect_genes <- 10
n_null <- 2000       # matched genes shared by all subtypes
n_boot <- 1000       # donor bootstrap replicates
max_rank <- 1500     # UCell
norm_scale <- 1e4

theme_set(theme_bw(base_size = 12))
orange <- "#D9702E"
caption_theme <- theme(plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot", plot.title.position = "plot")
wrap <- function(...) paste(strwrap(paste(...), width = 110), collapse = "\n")
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
}
no_overwrite <- function(file) if (file.exists(file)) stop(file, " exists; move it aside first")

read_gmt <- function(file) {
  lines <- strsplit(readLines(file), "\t")
  setNames(lapply(lines, function(x) x[-(1:2)]), vapply(lines, `[`, "", 1))
}
clipz <- function(r) atanh(pmin(pmax(r, -0.9999), 0.9999))

# DerSimonian-Laird random-effects combination across donors. Z: donors x genes (Fisher z), V: variance per donor. Returns one r per gene.
dl_r <- function(Z, V) {
  k <- nrow(Z)
  w <- 1 / V
  sw <- sum(w)
  z_fixed <- colSums(w * Z) / sw
  Q <- colSums(w * sweep(Z, 2, z_fixed)^2)
  tau2 <- pmax(0, (Q - (k - 1)) / (sw - sum(w^2) / sw))
  w_re <- 1 / outer(V, tau2, "+")
  tanh(colSums(w_re * Z) / colSums(w_re))
}
