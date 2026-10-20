# Transfer the 18 Trobisch neuronal subtypes to Schirmer neurons with SingleR, then assign coarse groups by label name.
# Needs results/schirmer_neurons.rds and results/trobisch_neurons.rds from 00_prepare_inputs.R.
source("common.R")
library(SingleR)
library(scuttle)
library(SingleCellExperiment)

out <- res("label_transfer_counts.csv")
no_overwrite(out)

s <- readRDS(res("schirmer_neurons.rds"))
t <- readRDS(res("trobisch_neurons.rds"))

# Coarse group map for all 18 reference labels
labels <- sort(unique(t$cells$label))
stopifnot(length(labels) == 18)
group_map <- data.frame(label = labels, coarse_group = coarse_group(labels))
group_map$coarse_group[is.na(group_map$coarse_group)] <- "excluded"
write.csv(group_map, res("coarse_group_map.csv"), row.names = FALSE)
print(group_map)

# Reference: log-normalized Trobisch counts aggregated to pseudobulk per label (defaults); query: log-normalized Schirmer neurons
ref <- logNormCounts(SingleCellExperiment(list(counts = t$counts), colData = DataFrame(label = t$cells$label)))
agg <- aggregateReference(ref, ref$label)
query <- logNormCounts(SingleCellExperiment(list(counts = s$counts)))
pred <- SingleR(test = logcounts(query), ref = agg, labels = agg$label, de.method = "wilcox")
stopifnot(identical(rownames(pred), s$cells$cell))

cells <- s$cells
cells$label <- ifelse(is.na(pred$pruned.labels), "unassigned", pred$pruned.labels)
cells$coarse_group <- group_map$coarse_group[match(cells$label, group_map$label)]
cells$coarse_group[cells$label == "unassigned"] <- "unassigned"
write.csv(cells, res("labels.csv"), row.names = FALSE)

# Donor x label x coarse group counts
counts <- as.data.frame(table(donor = cells$donor, label = cells$label), stringsAsFactors = FALSE)
counts$coarse_group <- cells$coarse_group[match(counts$label, cells$label)]
counts$condition <- cells$condition[match(counts$donor, cells$donor)]
counts$shared_with_trobisch <- counts$donor %in% shared_donors
write.csv(counts, out, row.names = FALSE)

cat("Schirmer neurons:", nrow(cells), " unassigned:", sum(cells$label == "unassigned"), sprintf("(%.1f%%)\n", 100 * mean(cells$label == "unassigned")))
print(table(cells$coarse_group))
print(round(prop.table(table(cells$donor, cells$coarse_group), 1), 2))
print(table(cells$donor, cells$coarse_group))
print(table(author = cells$author_label, transferred = cells$coarse_group))
