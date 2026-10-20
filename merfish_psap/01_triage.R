# Triage of GSE284005 for the PSAP check: panel, neurons, donors, data type, annotation, donor overlap. Needs data/merfish_input.rds.
# Reads (read only): ../rank_genes.csv for the Rank 2 module (Rank 2, use_primary_analysis == "Y", unique, TBC13 -> TBC1D13).
source("common.R")

d <- readRDS("data/merfish_input.rds")
meta <- d$meta
tri <- file.path(res_dir, "triage")
dir.create(tri, showWarnings = FALSE)

rg <- read.csv(rank_genes_csv)
rank2 <- unique(rg$human_symbol[rg$rank == "Rank 2" & rg$use_primary_analysis == "Y"])
rank2[rank2 == "TBC13"] <- "TBC1D13"
stopifnot(length(rank2) == 31)

# 1: panel, Rank 2 genes and markers on it
mk <- stack(markers)
panel <- data.frame(gene = sort(rownames(d$counts)))
panel$in_rank2 <- panel$gene %in% rank2
panel$marker_class <- vapply(panel$gene, function(g) paste(mk$ind[mk$values == g], collapse = ";"), "")
write.csv(panel, file.path(tri, "panel_genes.csv"), row.names = FALSE)
rank2_panel <- data.frame(gene = rank2, on_panel = rank2 %in% rownames(d$counts))
write.csv(rank2_panel, file.path(tri, "rank2_on_panel.csv"), row.names = FALSE)
markers_panel <- data.frame(marker_class = mk$ind, gene = mk$values, on_panel = mk$values %in% rownames(d$counts))
write.csv(markers_panel, file.path(tri, "markers_on_panel.csv"), row.names = FALSE)
cat("panel:", nrow(panel), "genes; PSAP on panel:", "PSAP" %in% panel$gene, "; ITM2C:", "ITM2C" %in% panel$gene,
    "; Rank 2 on panel:", sum(rank2_panel$on_panel), "of", length(rank2), "\n")

# 2: neurons per sample and region
neu <- meta[meta$major == neuron_label, ]
by_region <- as.data.frame.matrix(table(neu$sample, neu$region))
by_region <- cbind(sample = rownames(by_region), by_region)
write.csv(by_region, file.path(tri, "neurons_by_donor_region.csv"), row.names = FALSE)
cat("neurons:", nrow(neu), "; in GM:", sum(neu$region == gm_label), "\n")

# 3: donors
donors <- d$samples[, c("sample", "donor", "condition", "sex", "age", "ms_duration_years", "tissue")]
donors$gm_neurons <- as.integer(table(factor(neu$sample[neu$region == gm_label], donors$sample)))
donors$gm_ge_min <- donors$gm_neurons >= min_gm_neurons
write.csv(donors, file.path(tri, "donors.csv"), row.names = FALSE)
gm_n <- table(donors$condition[donors$gm_ge_min])
cat("donors with >=", min_gm_neurons, "GM neurons:", paste(names(gm_n), gm_n, collapse = ", "), "\n")

# 4: data type
x <- d$counts
data_type <- data.frame(item = c("counts_all_integer", "counts_max", "median_transcripts_per_cell", "coordinate_columns", "transcript_level_coordinates"),
                        value = c(all(x@x == round(x@x)), max(x@x), median(Matrix::colSums(x)), "x,y (cell centroids)", FALSE))
write.csv(data_type, file.path(tri, "data_type.csv"), row.names = FALSE)

# 5: annotation
write.csv(as.data.frame(table(region = meta$region, major = meta$major)), file.path(tri, "region_by_celltype.csv"), row.names = FALSE)
write.csv(as.data.frame(table(region = meta$region, sample = meta$sample)), file.path(tri, "cells_by_region_sample.csv"), row.names = FALSE)
write.csv(as.data.frame(table(region = neu$region, sub = neu$sub)), file.path(tri, "neuron_subtypes_by_region.csv"), row.names = FALSE)

# 6: donor overlap by ID (the study gives no donor IDs, so this can only exclude matching labels)
overlap <- data.frame(sample = donors$sample, matches_trobisch_id = toupper(donors$sample) %in% trobisch_donors)
write.csv(overlap, file.path(tri, "donor_overlap.csv"), row.names = FALSE)

# go/no-go
n_gm_donors <- sum(donors$gm_neurons > 0)
n_gm_ctrl <- sum(donors$gm_ge_min & donors$condition == "Control")
go <- data.frame(criterion = c("PSAP on panel", "raw integer counts", "donors with GM neurons >= 3", "GM controls with >= min neurons",
                               "transcript coordinates"),
                 value = c("PSAP" %in% panel$gene, data_type$value[1], n_gm_donors, n_gm_ctrl, FALSE),
                 pass = c("PSAP" %in% panel$gene, data_type$value[1] == "TRUE", n_gm_donors >= 3, n_gm_ctrl >= 3, FALSE))
write.csv(go, file.path(tri, "go_no_go.csv"), row.names = FALSE)
print(go)
