# Shared setup for the pseudotime scripts. Run them from scRNAseq/pseudotime/.
# GENESET picks the gene set: a rank label in ../../rank_genes.csv (default "Rank 2") or a file of symbols (txt, one per line, or csv with a gene / human_symbol column).
# MIN_DETECT is the smallest fraction of nuclei with a count for a gene to enter the module (default 0.05).
# 00 and 01 write to results/; 02 and 03 write to results/<geneset>/ and figures/<geneset>/.
library(ggplot2)

set.seed(1010)
data_dir <- "../data/Schirmer_2019"
cell_types <- c("EN-L2-3-A", "EN-L2-3-B")   # the authors' pseudotime exists for these nuclei only
min_cells <- 20                              # nuclei per donor for a donor to be counted
geneset_arg <- Sys.getenv("GENESET", "Rank 2")
min_detect <- as.numeric(Sys.getenv("MIN_DETECT", "0.05"))
geneset_name <- gsub("[^A-Za-z0-9]+", "_", sub("\\.[A-Za-z]+$", "", basename(geneset_arg)))
res_dir <- file.path("results", geneset_name)
fig_dir <- file.path("figures", geneset_name)
res <- function(file) {
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(res_dir, file)
}
shared <- function(file) file.path("results", file)

theme_set(theme_bw(base_size = 12))
orange <- "#D9702E"
blue <- "#2F78A8"
wrap <- function(...) paste(strwrap(paste(...), width = 150), collapse = "\n")
save_fig <- function(p, name, width, height) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height, bg = "white")
}

# Gene symbols of the chosen set, upper-cased and de-duplicated
load_geneset <- function() {
  if (file.exists(geneset_arg)) {
    if (grepl("\\.(csv|tsv)$", geneset_arg)) {
      x <- read.csv(geneset_arg, sep = if (grepl("\\.tsv$", geneset_arg)) "\t" else ",")
      col <- intersect(c("gene", "human_symbol"), names(x))
      genes <- x[[if (length(col)) col[1] else 1]]
    } else {
      genes <- readLines(geneset_arg)
    }
  } else {
    ranks <- read.csv("../../rank_genes.csv")
    stopifnot("GENESET is neither a file nor a rank label in rank_genes.csv" = geneset_arg %in% ranks$rank)
    genes <- ranks$human_symbol[ranks$rank == geneset_arg & ranks$use_primary_analysis == "Y"]
  }
  genes <- unique(toupper(trimws(genes)))
  genes[!is.na(genes) & genes != ""]
}
