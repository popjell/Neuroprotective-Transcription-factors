# Cache Schirmer neurons (all author EN-/IN- labels) and Trobisch cortex neurons (18 author subtypes) as donor-labeled sparse counts; write the input inventory.
# Needs ../data/Schirmer_2019 and ../ms-cross-regional/cor. All 19 Schirmer donors are kept; the 7 shared with Trobisch are flagged.
source("common.R")
library(Seurat)

out_s <- res("schirmer_neurons.rds")
out_t <- res("trobisch_neurons.rds")
no_overwrite(out_s)
no_overwrite(out_t)

# Schirmer: counts, metadata, sample -> donor (samples from one donor are pooled)
meta <- read.delim(file.path(schirmer_dir, "meta.tsv"))
barcodes <- readLines(file.path(schirmer_dir, "raw/barcodes.tsv"))
genes <- read.delim(file.path(schirmer_dir, "raw/genes.tsv"), header = FALSE, col.names = c("ensembl", "symbol"))
stopifnot(all(barcodes %in% meta$cell))
meta <- meta[match(barcodes, meta$cell), ]
counts <- as(readMM(file.path(schirmer_dir, "raw/matrix.mtx")), "CsparseMatrix")
stopifnot(ncol(counts) == nrow(meta), nrow(counts) == nrow(genes))
integer_counts <- all(counts@x == round(counts@x))
stopifnot(integer_counts)

donors <- read.csv(file.path(schirmer_dir, "sample_donor.csv"))
stopifnot(all(meta$sample %in% donors$sample))
meta$donor <- donors$donor[match(meta$sample, donors$sample)]
meta$condition <- ifelse(grepl("^MS", meta$sample), "MS", "Control")
stopifnot(identical(meta$condition == "MS", meta$diagnosis == "MS"))
stopifnot(!any(tapply(meta$condition, meta$donor, function(x) length(unique(x)) > 1)))

# Neurons by the author labels; genes sharing a symbol are summed
neuron_labels <- sort(grep("^(EN|IN)-", unique(meta$cell_type), value = TRUE))
keep <- meta$cell_type %in% neuron_labels
cells <- data.frame(cell = meta$cell[keep], author_label = meta$cell_type[keep], sample = meta$sample[keep],
                    donor = meta$donor[keep], condition = meta$condition[keep])
cells$shared_with_trobisch <- cells$donor %in% shared_donors
symbol <- factor(genes$symbol)
collapse <- sparseMatrix(i = as.integer(symbol), j = seq_along(symbol), x = 1, dimnames = list(levels(symbol), NULL))
x <- collapse %*% counts[, keep]
colnames(x) <- cells$cell
saveRDS(list(counts = x, cells = cells), out_s)
cat("Schirmer neurons:", ncol(x), "cells,", nrow(x), "genes,", length(unique(cells$donor)), "donors\n")
rm(counts, x)

# Trobisch cortex: neuron subtypes, donor from the sample mapping
tmeta <- read.delim(file.path(trobisch_dir, "meta.tsv"), check.names = FALSE)
tcounts <- GetAssayData(readRDS(file.path(trobisch_dir, "seurat.rds")), assay = "RNA", layer = "counts")
stopifnot(all(tmeta$Cell %in% colnames(tcounts)), all(tmeta$sample_id %in% names(trobisch_donor)))
tcounts <- tcounts[, tmeta$Cell]
stopifnot(all(tcounts@x == round(tcounts@x)), !anyDuplicated(rownames(tcounts)))
trobisch_labels <- c("L2-EN", "L3-EN", "L2-L3-EN", "L4-LRRK1-EN", "L4-GABRG1-EN", "L4-MME-EN", "L4-NPFFR2-EN",
                     "L5-L6-ADAMTSL1-EN", "L5-L6-FRMD6-EN", "L5-L6-HTR2C-EN", "L5-L6-NTNG2-EN", "motoEN",
                     "VIP-IN", "SST-IN", "PVALB-IN", "PVALB-CA8-IN", "RELN-IN", "KIT-IN")
stopifnot(all(trobisch_labels %in% tmeta$cell_subtypes))
tkeep <- tmeta$cell_subtypes %in% trobisch_labels
tcells <- data.frame(cell = tmeta$Cell[tkeep], label = tmeta$cell_subtypes[tkeep], sample = tmeta$sample_id[tkeep],
                     donor = unname(trobisch_donor[tmeta$sample_id[tkeep]]), condition = tmeta$condition[tkeep])
saveRDS(list(counts = tcounts[, tcells$cell], cells = tcells), out_t)
cat("Trobisch neurons:", nrow(tcells), "cells,", length(unique(tcells$donor)), "donors\n")

# Shared donors and inventory
shared <- data.frame(donor = shared_donors,
                     trobisch_samples = sapply(shared_donors, function(d) paste(names(trobisch_donor)[trobisch_donor == d], collapse = "+")),
                     schirmer_samples = sapply(shared_donors, function(d) paste(donors$sample[donors$donor == d], collapse = "+")))
write.csv(shared, res("shared_donors.csv"), row.names = FALSE)
spd <- aggregate(sample ~ donor + condition, unique(meta[, c("sample", "donor", "condition")]), length)
names(spd)[3] <- "n_samples"
spd$neurons <- as.integer(table(cells$donor)[spd$donor])
spd$shared_with_trobisch <- spd$donor %in% shared_donors
inv <- c("# Input inventory", "",
         "Schirmer 2019 (all 19 donors used; the 7 shared with Trobisch are flagged, not excluded)",
         paste0("- path: ", normalizePath(file.path(schirmer_dir, "raw/matrix.mtx")), " (", round(file.size(file.path(schirmer_dir, "raw/matrix.mtx")) / 1e9, 2),
                " GB), with meta.tsv and sample_donor.csv in ", normalizePath(schirmer_dir)),
         paste0("- cells: ", nrow(meta), "; genes: ", nrow(genes), "; neurons: ", nrow(cells)),
         paste0("- raw integer counts: ", integer_counts),
         paste0("- metadata columns: ", paste(names(meta)[names(meta) %in% names(read.delim(file.path(schirmer_dir, "meta.tsv"), nrows = 1))], collapse = ", ")),
         paste0("- neuronal author labels counted: ", paste(neuron_labels, collapse = ", ")),
         "- donor comes from sample_donor.csv; samples from one donor are summed", "",
         "Samples per donor", "", "| donor | condition | samples | neurons | shared with Trobisch |", "|---|---|---|---|---|",
         sprintf("| %s | %s | %d | %d | %s |", spd$donor, spd$condition, spd$n_samples, spd$neurons, spd$shared_with_trobisch), "",
         "Trobisch cortex (reference and comparator)",
         paste0("- path: ", normalizePath(trobisch_dir), "/seurat.rds and meta.tsv (label column cell_subtypes)"),
         paste0("- neuron cells: ", nrow(tcells), " in ", length(trobisch_labels), " labels; donors: ", length(unique(tcells$donor)), "; integer counts: TRUE"))
writeLines(inv, res("input_inventory.md"))
