# Part A: glial markers as negative controls in the existing neuron DESeq2 results, plus donor-level marker and PSAP CPM.
# Needs ../results/Cortex/<subtype>/deseq2.csv, ../results/fit_status.csv, ../cell_death/data/cortex_neurons.rds, ../../rank_genes.csv (all read only).
source("common.R")

# The 15 markers must not be in the 31-gene Rank 2 list
rank <- read.csv(rank_genes_csv)
rank2 <- unique(rank$human_symbol[rank$rank == "Rank 2" & rank$use_primary_analysis == "Y"])
stopifnot(length(rank2) == 31, !any(unlist(markers) %in% rank2))
cat("Rank 2 genes:", length(rank2), " markers in Rank 2:", sum(unlist(markers) %in% rank2), "\n")

# Markers in each fitted subtype's results, then the contamination rule
fitted <- read.csv("../results/fit_status.csv")
fitted <- fitted$subtype[fitted$region == "Cortex" & fitted$status == "fitted"]
de <- do.call(rbind, lapply(fitted, function(s) {
  d <- read.csv(file.path(baseline_dir, s, "deseq2.csv"))
  d$subtype <- s
  d
}))
marker_neg <- marker_table(de)
flags <- marker_flags(marker_neg)
write.csv(marker_neg, res("marker_negative_controls.csv"), row.names = FALSE)
write.csv(flags, res("marker_flags.csv"), row.names = FALSE)
print(flags, row.names = FALSE)
cat("subtype x glial type combinations flagged:", sum(flags$flag), "of", nrow(flags), "\n")

# Donor-level marker and PSAP CPM: counts pooled per donor within subtype (donors with at least min_cells cells)
neurons <- readRDS(neurons_rds)
cells <- neurons$cells
cpm <- function(genes, idx) sum(neurons$counts[rownames(neurons$counts) %in% genes, idx]) / sum(cells$nUMI[idx]) * cpm_scale
donor_cpm <- do.call(rbind, lapply(fitted, function(s) {
  k <- cells[cells$subtype == s, ]
  n <- table(k$donor)
  do.call(rbind, lapply(names(n)[n >= min_cells], function(dn) {
    idx <- which(cells$subtype == s & cells$donor == dn)
    data.frame(donor = dn, condition = cells$condition[idx[1]], subtype = s, n_cells = length(idx),
               microglia_cpm = cpm(markers$Microglia, idx), astro_cpm = cpm(markers$Astrocyte, idx),
               oligo_cpm = cpm(markers$Oligodendrocyte, idx), psap_cpm = cpm("PSAP", idx))
  }))
}))
write.csv(donor_cpm, res("donor_marker_psap.csv"), row.names = FALSE)

# Descriptive only: Spearman rho across donors, no p-values
rho <- do.call(rbind, lapply(split(donor_cpm, donor_cpm$subtype), function(d)
  data.frame(subtype = d$subtype[1], n_donors = nrow(d), spearman_rho_microglia_psap = cor(d$microglia_cpm, d$psap_cpm, method = "spearman"))))
write.csv(rho, res("spearman_microglia_psap.csv"), row.names = FALSE)
print(rho, row.names = FALSE)
