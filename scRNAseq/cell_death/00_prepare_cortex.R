# Trobisch cortex neurons: eligible donor x subtype cells, per-cell covariates, detection by subtype -> data/cortex_neurons.rds
source("common.R")
library(Seurat)

out <- "data/cortex_neurons.rds"
no_overwrite(out)

counts <- GetAssayData(readRDS("../ms-cross-regional/cor/seurat.rds"), assay = "RNA", layer = "counts")
meta <- read.delim("../ms-cross-regional/cor/meta.tsv", check.names = FALSE)
donors <- c(CR_ctrl1 = "CO39", CR_ctrl2 = "CO14", CR_ctrl3 = "CO36", CR_ctrl4 = "CO25",
            CR_ms1 = "MS60", CR_ms2 = "MS229", CR_ms3 = "MS371", CR_ms4 = "MS371")
neuron_types <- c("L2-EN", "L3-EN", "L2-L3-EN", "L4-LRRK1-EN", "L4-GABRG1-EN", "L4-MME-EN", "L4-NPFFR2-EN",
                  "L5-L6-ADAMTSL1-EN", "L5-L6-FRMD6-EN", "L5-L6-HTR2C-EN", "L5-L6-NTNG2-EN", "motoEN",
                  "VIP-IN", "SST-IN", "PVALB-IN", "PVALB-CA8-IN", "RELN-IN", "KIT-IN")

cells <- data.frame(cell = meta$Cell, subtype = meta$cell_subtypes, sample_id = meta$sample_id, condition = meta$condition)
cells <- cells[cells$subtype %in% neuron_types, ]
cells$donor <- unname(donors[cells$sample_id])

# Subtypes with at least min_donors control and MS donors holding at least min_cells cells
per_donor <- as.data.frame(table(donor = cells$donor, subtype = cells$subtype), stringsAsFactors = FALSE)
per_donor$condition <- cells$condition[match(per_donor$donor, cells$donor)]
per_donor$eligible <- per_donor$Freq >= min_cells
enough <- sapply(unique(per_donor$subtype), function(s) {
  e <- per_donor[per_donor$subtype == s & per_donor$eligible, ]
  sum(e$condition == "Control") >= min_donors && sum(e$condition == "MS") >= min_donors
})
subtypes <- sort(unique(per_donor$subtype)[enough])
keep <- cells$subtype %in% subtypes & paste(cells$donor, cells$subtype) %in% paste(per_donor$donor, per_donor$subtype)[per_donor$eligible]
analysis <- cells[keep, ]

x <- counts[, analysis$cell]
analysis$nUMI <- colSums(x)
analysis$ngen <- colSums(x > 0)
# Fraction of cells detecting each gene, from all cells of the subtype
frac <- sapply(subtypes, function(s) rowMeans(counts[, cells$cell[cells$subtype == s], drop = FALSE] > 0))

cat("subtypes:", length(subtypes), " cells:", nrow(analysis), " donors:", length(unique(analysis$donor)), "\n")
saveRDS(list(counts = x, cells = analysis, frac = frac, subtypes = subtypes), out)
