# Prioritize glial ligands linked to MS-up neuronal genes in the full atlas.
library(Matrix)
library(nichenetr)
library(ggplot2)

dir.create("results/ligand", showWarnings = FALSE, recursive = TRUE)
lr <- readRDS("data/lr_network_human_21122021.rds")
ligands <- unique(lr$from)
sender_types <- c("MG", "GM-AS", "WM-AS", "WM-AS-GFPT2", "BG")
neurons <- unique(read.delim("data/meta.tsv")$cell_subtypes_intg)
neurons <- neurons[grepl("-EN$|-IN$|^PC$|^motoEN$", neurons)]

# Record ligand detection in each glial group, separately by region and disease.
sender_rows <- list()
for (source in c("data/cortex_senders.rds", "data/new_regions.rds")) {
  x <- readRDS(source)
  counts <- x$counts
  meta <- x$meta
  for (region in unique(meta$cns_region)) {
    for (type in sender_types) {
      for (condition in c("Control", "MS")) {
        cells <- which(meta$cns_region == region &
                         meta$cell_subtypes_intg == type & meta$condition == condition)
        if (length(cells) < 20) next
        genes <- intersect(ligands, rownames(counts))
        z <- counts[genes, cells, drop = FALSE]
        donor_ids <- unique(meta$sample_id[cells])
        donor_n <- rowSums(vapply(donor_ids, function(id) {
          rowSums(z[, meta$sample_id[cells] == id, drop = FALSE]) > 0
        }, logical(length(genes))))
        sender_rows[[length(sender_rows) + 1]] <- data.frame(
          region = region, sender = type, condition = condition,
          ligand = genes, cells = length(cells), donors = length(donor_ids),
          fraction = rowSums(z > 0) / length(cells), donors_detected = donor_n)
      }
    }
  }
  rm(x, counts, meta)
  gc()
}
sender <- do.call(rbind, sender_rows)
write.csv(sender, "results/ligand/sender_expression.csv", row.names = FALSE)

# Keep ligand candidates more detectable in glia than in pooled neurons.
neuron_rows <- list()
for (source in c("data/cortex_neurons.rds", "data/new_regions.rds")) {
  x <- readRDS(source)
  for (region in unique(x$meta$cns_region)) {
    cells <- which(x$meta$cns_region == region &
                     x$meta$cell_subtypes_intg %in% neurons &
                     x$meta$condition == "MS")
    if (length(cells) < 100) next
    genes <- intersect(ligands, rownames(x$counts))
    neuron_rows[[length(neuron_rows) + 1]] <- data.frame(
      region = region, ligand = genes,
      neuron_fraction = rowSums(x$counts[genes, cells, drop = FALSE] > 0) /
        length(cells))
  }
  rm(x)
  gc()
}
neuron_fraction <- do.call(rbind, neuron_rows)
write.csv(neuron_fraction, "results/ligand/neuron_expression.csv", row.names = FALSE)
sender <- merge(sender, neuron_fraction, by = c("region", "ligand"))

# Use the saved donor-level neuronal counts and MS-vs-control results.
pb_cortex <- readRDS("results/cortex_integrated_pseudobulk.rds")
pb_new <- readRDS("results/integrated_neuron_pseudobulk.rds")
de_cortex <- read.csv("results/cortex_integrated_de.csv")
de_cortex$region <- "Cortex"
de_new <- read.csv("results/integrated_neuron_de.csv")
de <- rbind(de_cortex[, names(de_new)], de_new)
rank <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank <- rank[rank$use_primary_analysis == "Y", c("rank", "human_symbol")]
rank_label <- tapply(rank$rank, rank$human_symbol, paste, collapse = ", ")

# Read the reference once; each neuronal group is scored against the same network.
lt <- readRDS("data/ligand_target_matrix_nsga2r_final.rds")
raw_cortex <- readRDS("data/cortex_neurons.rds")
raw_new <- readRDS("data/new_regions.rds")
status_rows <- list()
activity_rows <- list()
link_rows <- list()
support_rows <- list()
receptor_rows <- list()
for (region in c("Cortex", "Cerebellum", "Spinal_cord")) {
  pb <- if (region == "Cortex") pb_cortex else pb_new
  metadata <- pb$metadata
  if (region != "Cortex") metadata <- metadata[metadata$cns_region == region, ]
  for (type in neurons) {
    group <- metadata[metadata$subtype == type & metadata$cells >= 10, ]
    n_control <- sum(group$condition == "Control")
    n_ms <- sum(group$condition == "MS")
    result <- data.frame(region = region, subtype = type,
                         control_donors = n_control, ms_donors = n_ms,
                         targets = 0, candidate_ligands = 0,
                         status = "too few donor samples")
    if (n_control < 2 || n_ms < 2) {
      status_rows[[length(status_rows) + 1]] <- result
      next
    }
    ms <- rownames(group)[group$condition == "MS"]
    z <- pb$counts[, ms, drop = FALSE]
    cpm <- t(t(z) / colSums(z) * 1e6)
    expressed <- rownames(z)[rowSums(cpm > 1) >= 2]
    result$status <- "too few MS-up genes"
    d <- de[de$region == region & de$subtype == type, ]
    targets <- intersect(d$gene[d$FDR < 0.05 & d$logFC > 0], expressed)
    targets <- intersect(targets, rownames(lt))
    result$targets <- length(targets)
    if (length(targets) < 20) {
      status_rows[[length(status_rows) + 1]] <- result
      next
    }
    source_ligands <- unique(sender$ligand[
      sender$region == region & sender$condition == "MS" &
        sender$donors_detected >= 2 & sender$fraction >= 0.10 &
        sender$fraction >= 2 * sender$neuron_fraction])
    raw_neurons <- if (region == "Cortex") raw_cortex else raw_new
    raw_cells <- which(raw_neurons$meta$cns_region == region &
                         raw_neurons$meta$cell_subtypes_intg == type &
                         raw_neurons$meta$condition == "MS")
    receptors <- intersect(expressed, lr$to)
    receptor_fraction <- rowSums(raw_neurons$counts[
      receptors, raw_cells, drop = FALSE] > 0) / length(raw_cells)
    receptor_rows[[length(receptor_rows) + 1]] <- data.frame(
      region = region, subtype = type, receptor = receptors,
      fraction = receptor_fraction, cells = length(raw_cells))
    receptors <- receptors[receptor_fraction >= 0.05]
    possible <- intersect(source_ligands, lr$from[lr$to %in% receptors])
    possible <- intersect(possible, colnames(lt))
    result$candidate_ligands <- length(possible)
    result$status <- if (length(possible) >= 5) "scored" else "too few ligand-receptor pairs"
    status_rows[[length(status_rows) + 1]] <- result
    if (result$status != "scored") next
    background <- intersect(expressed, rownames(lt))
    activity <- predict_ligand_activities(targets, background, lt, possible)
    activity$region <- region
    activity$subtype <- type
    activity$targets <- length(targets)
    activity_rows[[length(activity_rows) + 1]] <- activity
    best <- head(activity$test_ligand[order(-activity$aupr_corrected)], 20)
    support <- merge(lr[lr$from %in% best & lr$to %in% receptors,
                        c("from", "to")],
                     sender[sender$region == region & sender$condition == "MS" &
                              sender$donors_detected >= 2 & sender$fraction >= 0.10 &
                              sender$fraction >= 2 * sender$neuron_fraction,
                            c("sender", "ligand", "fraction", "neuron_fraction",
                              "donors_detected")],
                     by.x = "from", by.y = "ligand")
    support$region <- region
    support$subtype <- type
    support_rows[[length(support_rows) + 1]] <- unique(support)
    rank_up <- intersect(rank$human_symbol, d$gene[d$logFC > 0])
    rank_up <- intersect(rank_up, expressed)
    for (ligand in best) {
      links <- get_weighted_ligand_target_links(ligand, rank_up, lt, n = 250)
      links <- links[!is.na(links$target), ]
      if (nrow(links) == 0) next
      links$region <- region
      links$subtype <- type
      links$rank <- rank_label[links$target]
      links$target_logFC <- d$logFC[match(links$target, d$gene)]
      links$target_FDR <- d$FDR[match(links$target, d$gene)]
      links$ligand_aupr_corrected <- activity$aupr_corrected[
        match(ligand, activity$test_ligand)]
      link_rows[[length(link_rows) + 1]] <- links
    }
    message(region, " / ", type, ": ", length(targets), " targets; ",
            length(possible), " candidate ligands")
  }
}
status <- do.call(rbind, status_rows)
write.csv(status, "results/ligand/subtype_status.csv", row.names = FALSE)
if (length(receptor_rows))
  write.csv(do.call(rbind, receptor_rows),
            "results/ligand/receptor_expression.csv", row.names = FALSE)
if (length(activity_rows)) {
  activity <- do.call(rbind, activity_rows)
  activity <- activity[order(activity$region, activity$subtype,
                             -activity$aupr_corrected), ]
  write.csv(activity, "results/ligand/ligand_activity.csv", row.names = FALSE)
  if (length(link_rows)) {
    links <- do.call(rbind, link_rows)
    write.csv(links, "results/ligand/rank_gene_links.csv", row.names = FALSE)
  }
  if (length(support_rows))
    write.csv(do.call(rbind, support_rows),
              "results/ligand/ligand_receptor_support.csv", row.names = FALSE)
  top <- do.call(rbind, lapply(split(activity, interaction(activity$region,
                                     activity$subtype, drop = TRUE)), head, 5))
  top$group <- paste(top$region, top$subtype, sep = " / ")
  p <- ggplot(top, aes(x = test_ligand, y = aupr_corrected)) +
    geom_col(fill = "#345f72") + coord_flip() +
    facet_wrap(~group, scales = "free_y", ncol = 4) +
    labs(x = NULL, y = "NicheNet corrected AUPR",
         title = "Candidate glial ligands for MS-up neuronal genes") +
    theme_bw(base_size = 10)
  ggsave("results/ligand/top_ligands.png", p, width = 13, height = 12, dpi = 180)
}
