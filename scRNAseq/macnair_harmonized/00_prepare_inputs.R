# Cache Macnair grey matter neurons as one sparse counts file per sample, and the Trobisch cortex neurons (reference); write the sample table.
# Needs data/ (Zenodo 8338963 files) and ../ms-cross-regional/cor. Streams the 4.6 GB matrix once.
source("common.R")
library(Seurat)

out_samples <- res("samples.csv")
out_t <- res("trobisch_neurons.rds")
no_overwrite(out_samples)
no_overwrite(out_t)
dir.create(res("cache"), showWarnings = FALSE)

# Metadata and genes; row order of col_data is the column order of the matrix
mtx <- file.path(macnair_dir, "ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz")
meta <- read.csv(gzfile(file.path(macnair_dir, "ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz")))
genes <- read.csv(gzfile(file.path(macnair_dir, "ms_lesions_snRNAseq_row_data_2023-09-12.txt.gz")))
header <- readLines(gzfile(mtx), n = 2)
dims <- as.integer(strsplit(header[2], " ")[[1]])
integer_counts <- grepl("integer", header[1])
stopifnot(integer_counts, dims[1] == nrow(genes), dims[2] == nrow(meta), !anyDuplicated(meta$cell_id))

# Grey matter: matter == GM and lesion_type is control GM, NAGM or GML (the two must agree); neurons by Macnair type_broad
stopifnot(identical(meta$matter == "GM", meta$lesion_type %in% c("GM", "NAGM", "GML")))
meta$condition <- ifelse(meta$diagnosis == "CTR", "Control", "MS")
meta$tissue <- c(GM = "control GM", NAGM = "NAGM", GML = "GML")[meta$lesion_type]
neuron_types <- c("Excitatory neurons", "Inhibitory neurons")
stopifnot(neuron_types %in% meta$type_broad)
keep <- which(meta$matter == "GM" & meta$type_broad %in% neuron_types)
cells <- meta[keep, ]
cells$col <- keep
stopifnot(!any(tapply(cells$condition, cells$individual_id_anon, function(x) length(unique(x)) > 1)),
          all(cells$tissue[cells$condition == "Control"] == "control GM"), all(cells$tissue[cells$condition == "MS"] != "control GM"))
samples <- unique(cells$sample_id_anon)  # order of appearance in the matrix columns
cells$sample_index <- match(cells$sample_id_anon, samples)
write.table(cells[, c("col", "sample_index")], res("kept_columns.txt"), row.names = FALSE, col.names = FALSE)

# Gene symbols shared by several rows are summed
symbol <- factor(genes$symbol)
collapse <- sparseMatrix(i = as.integer(symbol), j = seq_along(symbol), x = 1, dimnames = list(levels(symbol), NULL))

# Write one sample: genes x cells (symbols summed); library size must equal the author's log10_counts
flush <- function(k, parts) {
  cs <- cells[cells$sample_index == k, ]
  p <- do.call(rbind, parts)
  x <- sparseMatrix(i = p$gene, j = match(p$col, cs$col), x = p$count, dims = c(nrow(genes), nrow(cs)))
  x <- collapse %*% x
  colnames(x) <- cs$cell_id
  stopifnot(all(abs(log10(colSums(x)) - cs$log10_counts) < 1e-6))
  out <- file.path(res("cache"), paste0(samples[k], ".rds"))
  no_overwrite(out)
  saveRDS(list(counts = x, cells = cs), out)
  cat("  ", samples[k], ncol(x), "neurons\n")
}

# Stream the matrix through awk (kept columns only); samples must arrive one after another
cmd <- paste0("zcat ", mtx, " | awk 'NR==FNR{k[$1]=$2; next} /^%/{next} !h{h=1; next} ($2 in k){print $1, $2, $3, k[$2]}' ",
              res("kept_columns.txt"), " -")
con <- pipe(cmd, "r")
pending <- list()
done <- integer()
repeat {
  chunk <- scan(con, what = list(gene = 0L, col = 0L, count = 0, sample = 0L), nmax = chunk_records, quiet = TRUE)
  if (!length(chunk$gene)) break
  stopifnot(!is.unsorted(chunk$sample), !any(unique(chunk$sample) %in% done))
  for (k in unique(chunk$sample)) {
    i <- chunk$sample == k
    pending[[as.character(k)]] <- c(pending[[as.character(k)]], list(data.frame(gene = chunk$gene[i], col = chunk$col[i], count = chunk$count[i])))
  }
  for (k in setdiff(names(pending), as.character(max(chunk$sample)))) {
    flush(as.integer(k), pending[[k]])
    done <- c(done, as.integer(k))
    pending[[k]] <- NULL
  }
}
close(con)
for (k in names(pending)) flush(as.integer(k), pending[[k]])
stopifnot(length(list.files(res("cache"))) == length(samples))

# Sample table (donor-level covariates are constant within donor)
d <- unique(cells[, c("individual_id_anon", "diagnosis", "sample_source", "sex", "age_at_death", "pmi_minutes")])
stopifnot(!anyDuplicated(d$individual_id_anon))
s <- unique(cells[, c("sample_id_anon", "individual_id_anon", "condition", "diagnosis", "tissue", "sample_source", "sex", "age_at_death", "pmi_minutes")])
names(s) <- c("sample", "donor", "condition", "diagnosis", "tissue", "bank", "sex", "age", "pmi_minutes")
s$n_neurons <- as.integer(table(cells$sample_id_anon)[s$sample])
s$n_flagged_exclude_pseudobulk <- as.integer(tapply(cells$exclude_pseudobulk, cells$sample_id_anon, sum)[s$sample])
write.csv(s, out_samples, row.names = FALSE)
writeLines(c(paste("grey matter cells:", sum(meta$matter == "GM")), paste("raw integer counts:", integer_counts),
             paste("matrix dims:", paste(dims, collapse = " x ")), paste("neuron types:", paste(neuron_types, collapse = ", "))), res("macnair_dims.txt"))
cat("Macnair GM neurons:", nrow(cells), "cells,", nrow(s), "samples,", nrow(d), "donors\n")
rm(cells, meta)

# Trobisch cortex neurons: reference for SingleR
tmeta <- read.delim(file.path(trobisch_dir, "meta.tsv"), check.names = FALSE)
tcounts <- GetAssayData(readRDS(file.path(trobisch_dir, "seurat.rds")), assay = "RNA", layer = "counts")
stopifnot(all(tmeta$Cell %in% colnames(tcounts)))
tcounts <- tcounts[, tmeta$Cell]
stopifnot(all(tcounts@x == round(tcounts@x)), !anyDuplicated(rownames(tcounts)))
trobisch_labels <- c("L2-EN", "L3-EN", "L2-L3-EN", "L4-LRRK1-EN", "L4-GABRG1-EN", "L4-MME-EN", "L4-NPFFR2-EN",
                     "L5-L6-ADAMTSL1-EN", "L5-L6-FRMD6-EN", "L5-L6-HTR2C-EN", "L5-L6-NTNG2-EN", "motoEN",
                     "VIP-IN", "SST-IN", "PVALB-IN", "PVALB-CA8-IN", "RELN-IN", "KIT-IN")
stopifnot(all(trobisch_labels %in% tmeta$cell_subtypes))
tkeep <- tmeta$cell_subtypes %in% trobisch_labels
saveRDS(list(counts = tcounts[, tmeta$Cell[tkeep]], labels = tmeta$cell_subtypes[tkeep]), out_t)
cat("Trobisch reference:", sum(tkeep), "neurons in", length(trobisch_labels), "labels\n")
