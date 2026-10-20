# Step 1 inspection: format, dimensions, integer check, metadata columns, donor table. Writes results/inspect.txt. Needs the figshare files in _input/.
source("common.R")

out <- file.path(res_dir, "inspect.txt")
no_overwrite(out)
sink(out, split = TRUE)

cat("== files ==\n")
print(file.info(list.files(in_dir, full.names = TRUE))[, "size", drop = FALSE])
cat("\nmtx header:\n")
print(readLines(gzfile(file.path(in_dir, "RGC_RNA_matrix.mtx.gz")), n = 3))

x <- load_counts()
genes <- attr(x, "genes")
cells <- attr(x, "cells")
cat("\n== matrix ==\n")
cat("dims (genes x cells):", dim(x), "\n")
cat("nonzero entries:", length(x@x), "\n")
cat("non-integer values:", sum(x@x != round(x@x)), "\n")
cat("negative values:", sum(x@x < 0), "\n")
cat("count range:", range(x@x), "\n")
cat("feature file columns:", ncol(genes), "\n")
print(head(genes, 5))
sym <- genes[[ncol(genes) - 1]]
cat("first gene ids:", head(genes[[1]], 5), "\n")
cat("duplicated symbols (col 2):", sum(duplicated(genes[[2]])), "\n")

meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), check.names = FALSE, stringsAsFactors = FALSE)
cat("\n== metadata ==\n")
cat("rows:", nrow(meta), " columns:", ncol(meta), "\n")
for (nm in names(meta)) {
  v <- meta[[nm]]
  u <- unique(v)
  cat("\n[", nm, "] class=", class(v)[1], " n_unique=", length(u), " n_NA=", sum(is.na(v)), "\n", sep = "")
  if (length(u) <= 40) print(table(v, useNA = "ifany")) else if (is.numeric(v)) print(summary(v)) else cat("e.g.:", head(u, 8), "\n")
}

cat("\n== barcode / metadata matching ==\n")
bc_col <- names(meta)[sapply(meta, function(v) is.character(v) && all(cells[[1]] %in% v))]
cat("metadata columns containing all barcodes:", bc_col, "\n")

cat("\n== gene list presence ==\n")
rk <- read.csv(file.path(in_dir, "rank_genes.csv"), stringsAsFactors = FALSE)
g31 <- unique(rk$human_symbol[rk$rank == "Rank 2" & rk$use_primary_analysis == "Y"])
cat("Rank 2 primary genes:", length(g31), "\n")
syms <- genes[[2]]
print(data.frame(gene = g31, present = g31 %in% syms))
cat("TBC1D13 present:", "TBC1D13" %in% syms, "\n")
cat("microglia/astro/oligo markers present:\n")
mk <- c("P2RY12", "CSF1R", "C1QB", "CD74", "GFAP", "GLUL", "RLBP1", "MBP", "PLP1")
print(setNames(mk %in% syms, mk))

cat("\n== trobisch/schirmer donor id check ==\n")
tro <- c("CO39", "CO14", "CO36", "CO25", "MS60", "MS229", "MS371")
hits <- sapply(meta, function(v) any(as.character(v) %in% tro))
cat("columns containing any Trobisch/Schirmer id:", names(hits)[hits], "\n")
sink()
cat("wrote", out, "\n")
