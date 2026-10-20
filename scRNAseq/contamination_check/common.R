# Shared setup for the contamination_check scripts. Run them from scRNAseq/contamination_check/ (see run_all.sh).
# RUN (default trobisch_cortex) names the results/ and figures/ subfolders.
.libPaths(c("Rlib", .libPaths()))  # celda is installed here
library(ggplot2)
library(Matrix)

set.seed(1010)
run <- Sys.getenv("RUN", "trobisch_cortex")
res_dir <- file.path("results", run)
fig_dir <- file.path("figures", run)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
res <- function(file) file.path(res_dir, file)

min_cells <- 20        # cells per donor x subtype
min_donors <- 3        # control donors and MS donors per subtype
min_count <- 10        # pseudobulk count threshold for the gene filter
min_pseudobulks <- 2   # pseudobulks that must reach min_count
cpm_scale <- 1e6       # counts per million
norm_scale <- 1e4      # log-normalised expression = log1p(count / nUMI * norm_scale)

# Inputs from other analyses (read only)
baseline_dir <- "../results/Cortex"                          # existing per-subtype DESeq2 results
neurons_rds <- "../cell_death/data/cortex_neurons.rds"       # neurons of eligible subtypes (column `cell`)
rank_genes_csv <- "../../rank_genes.csv"
rank_xlsx <- "../../Gene_Lists_Rank1_and_Rank2_Pathways.xlsx"
cortex_dir <- "../ms-cross-regional/cor"

donors <- c(CR_ctrl1 = "CO39", CR_ctrl2 = "CO14", CR_ctrl3 = "CO36", CR_ctrl4 = "CO25",
            CR_ms1 = "MS60", CR_ms2 = "MS229", CR_ms3 = "MS371", CR_ms4 = "MS371")  # MS371 is pooled

# Author cell_subtypes -> broad type (decontX z, and the neuron/glia split)
broad_of <- c(setNames(rep("Neuron", 18),
                       c("L2-EN", "L3-EN", "L2-L3-EN", "L4-LRRK1-EN", "L4-GABRG1-EN", "L4-MME-EN", "L4-NPFFR2-EN", "L5-L6-ADAMTSL1-EN",
                         "L5-L6-FRMD6-EN", "L5-L6-HTR2C-EN", "L5-L6-NTNG2-EN", "motoEN", "VIP-IN", "SST-IN", "PVALB-IN", "PVALB-CA8-IN",
                         "RELN-IN", "KIT-IN")),
              MG = "Microglia", "GM-AS" = "Astrocyte", "WM-AS" = "Astrocyte",
              "OL-A" = "Oligodendrocyte", "OL-Stress" = "Oligodendrocyte", "WM-OL" = "Oligodendrocyte", OPC = "OPC",
              END = "Vascular", PER = "Vascular", "B cells" = "Immune", "T cells" = "Immune")

# Fixed glial marker panel (negative controls)
markers <- list(Microglia = c("C1QA", "C1QB", "C1QC", "CSF1R", "P2RY12"),
                Astrocyte = c("GFAP", "AQP4", "ALDH1L1", "GJA1", "SLC1A3"),
                Oligodendrocyte = c("MBP", "PLP1", "MOBP", "MOG", "MAG"))

theme_set(theme_bw(base_size = 12))
orange <- "#D9702E"
blue <- "#2F78A8"
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
}
no_overwrite <- function(file) if (file.exists(file)) stop(file, " exists; move it aside first")

# Donor-level pseudobulk DESeq2, MS vs Control, one fit per subtype. cells needs subtype, donor, condition; counts is genes x cells.
# Returns the per-gene results and one status row per subtype.
run_de <- function(counts, cells) {
  suppressPackageStartupMessages(library(DESeq2))
  stopifnot(identical(colnames(counts), cells$cell), all(counts@x == round(counts@x)))
  out <- list()
  status <- list()
  for (s in sort(unique(cells$subtype))) {
    k <- cells[cells$subtype == s, ]
    n <- table(k$donor)
    use <- names(n)[n >= min_cells]
    k <- k[k$donor %in% use, ]
    cond <- k$condition[match(use, k$donor)]
    nc <- sum(cond == "Control")
    nm <- sum(cond == "MS")
    ok <- nc >= min_donors && nm >= min_donors
    status[[s]] <- data.frame(subtype = s, n_cells = nrow(k), control_donors = nc, ms_donors = nm, fitted = ok)
    if (!ok) next
    # Sum integer counts per donor
    ind <- sparseMatrix(i = seq_len(nrow(k)), j = match(k$donor, use), x = 1, dims = c(nrow(k), length(use)))
    y <- as.matrix(counts[, k$cell] %*% ind)
    colnames(y) <- use
    y <- y[rowSums(y >= min_count) >= min_pseudobulks, , drop = FALSE]
    d <- data.frame(row.names = use, condition = factor(cond, levels = c("Control", "MS")))
    dds <- DESeq(DESeqDataSetFromMatrix(y, d, ~condition), quiet = TRUE)
    r <- as.data.frame(results(dds, contrast = c("condition", "MS", "Control")))
    r$gene <- rownames(r)
    r$subtype <- s
    out[[s]] <- r
    cat(s, ": ", nc, " control vs ", nm, " MS donors, ", nrow(y), " genes tested\n", sep = "")
  }
  list(de = do.call(rbind, out), status = do.call(rbind, status))
}

# Glial marker table and the contamination rule for one set of DESeq2 results (columns subtype, gene, log2FoldChange, padj)
marker_table <- function(de) {
  m <- stack(markers)
  names(m) <- c("gene", "glial_type")
  m$glial_type <- as.character(m$glial_type)
  g <- merge(data.frame(subtype = unique(de$subtype)), m)
  g <- merge(g, de[, c("subtype", "gene", "log2FoldChange", "padj")], all.x = TRUE)
  g$status <- ifelse(is.na(g$log2FoldChange), "removed by count filter", "tested")
  names(g)[names(g) == "log2FoldChange"] <- "log2FC"
  g[order(g$subtype, g$glial_type, g$gene), c("subtype", "glial_type", "gene", "log2FC", "padj", "status")]
}
# Signal if, among markers that passed the filter, median log2FC > 0 and at least 2 have padj < 0.05
marker_flags <- function(mt) {
  do.call(rbind, lapply(split(mt, list(mt$subtype, mt$glial_type), drop = TRUE), function(d) {
    t <- d[d$status == "tested", ]
    n_sig <- sum(t$padj < 0.05, na.rm = TRUE)
    med <- if (nrow(t)) median(t$log2FC) else NA
    data.frame(subtype = d$subtype[1], glial_type = d$glial_type[1], n_markers_tested = nrow(t), median_log2FC = med,
               `n_padj_lt_0.05` = n_sig, flag = !is.na(med) && med > 0 && n_sig >= 2, check.names = FALSE)
  }))
}
