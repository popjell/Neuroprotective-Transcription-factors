# Transfer the 18 Trobisch neuronal subtypes to Macnair grey matter neurons with SingleR (31 Rank 2 genes removed), then assign coarse groups by label name.
# Needs results/cache/*.rds, results/samples.csv and results/trobisch_neurons.rds from 00_prepare_inputs.R.
source("common.R")
library(SingleR)
library(scuttle)
library(SingleCellExperiment)

out <- res("label_transfer_counts.csv")
no_overwrite(out)
no_overwrite(res("labels.csv"))
no_overwrite(res("pseudobulk_sample_group.rds"))

t <- readRDS(res("trobisch_neurons.rds"))
samples <- read.csv(res("samples.csv"))
drop <- rank2_genes()

# Coarse group map for the 18 reference labels
labels <- sort(unique(t$labels))
stopifnot(length(labels) == 18)
group_map <- data.frame(label = labels, coarse_group = coarse_group(labels))
group_map$coarse_group[is.na(group_map$coarse_group)] <- "excluded"
write.csv(group_map, res("coarse_group_map.csv"), row.names = FALSE)

# Reference without the Rank 2 genes: log-normalized counts aggregated to pseudobulk per label (defaults)
cat("Rank 2 genes found in the reference:", sum(drop %in% rownames(t$counts)), "of 31 (all removed)\n")
ref <- logNormCounts(SingleCellExperiment(list(counts = t$counts[!rownames(t$counts) %in% drop, ]), colData = DataFrame(label = t$labels)))
agg <- aggregateReference(ref, ref$label)

# One SingleR run per sample (each cell is scored on its own ranks, so this equals one joint run); raw counts keep all genes for pseudobulk
labelled <- list()
pb <- list()
for (sm in samples$sample) {
  x <- readRDS(res("cache", paste0(sm, ".rds")))
  stopifnot(nrow(x$cells) == samples$n_neurons[samples$sample == sm])
  query <- logNormCounts(SingleCellExperiment(list(counts = x$counts[!rownames(x$counts) %in% drop, ])))
  pred <- SingleR(test = logcounts(query), ref = agg, labels = agg$label, de.method = "wilcox")
  stopifnot(identical(rownames(pred), x$cells$cell_id))
  cells <- x$cells[, c("cell_id", "sample_id_anon", "individual_id_anon", "condition", "tissue", "sample_source", "type_broad", "type_fine")]
  cells$label <- ifelse(is.na(pred$pruned.labels), "unassigned", pred$pruned.labels)
  cells$coarse_group <- group_map$coarse_group[match(cells$label, group_map$label)]
  cells$coarse_group[cells$label == "unassigned"] <- "unassigned"
  labelled[[sm]] <- cells
  for (g in c("upper_EN", "deep_EN", "IN")) {
    i <- cells$coarse_group == g
    if (any(i)) pb[[paste(sm, g, sep = "|")]] <- list(counts = Matrix::rowSums(x$counts[, i, drop = FALSE]), n = sum(i))
  }
  cat(sm, nrow(cells), "neurons,", sum(cells$label == "unassigned"), "unassigned\n")
}
cells <- do.call(rbind, labelled)
rownames(cells) <- NULL
names(cells)[2:3] <- c("sample", "donor")
write.csv(cells, res("labels.csv"), row.names = FALSE)
saveRDS(list(counts = sapply(pb, `[[`, "counts"), n_cells = sapply(pb, `[[`, "n")), res("pseudobulk_sample_group.rds"))

# Sample x label x coarse group counts
counts <- as.data.frame(table(sample = cells$sample, label = cells$label), stringsAsFactors = FALSE)
counts$coarse_group <- cells$coarse_group[match(counts$label, cells$label)]
counts <- cbind(counts, samples[match(counts$sample, samples$sample), c("donor", "condition", "tissue", "bank")])
write.csv(counts, out, row.names = FALSE)
cat("Macnair GM neurons:", nrow(cells), " unassigned:", sum(cells$label == "unassigned"), sprintf("(%.1f%%)\n", 100 * mean(cells$label == "unassigned")))
print(table(cells$coarse_group))
