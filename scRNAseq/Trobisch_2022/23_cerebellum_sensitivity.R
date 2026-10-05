# Exploratory cerebellum run using nominal MS-up genes; not FDR-confirmed.
library(Matrix)
library(nichenetr)
library(ggplot2)

pb <- readRDS("results/integrated_neuron_pseudobulk.rds")
raw <- readRDS("data/new_regions.rds")
de <- read.csv("results/integrated_neuron_de.csv")
sender <- read.csv("results/ligand/sender_expression.csv")
neuron <- read.csv("results/ligand/neuron_expression.csv")
sender <- merge(sender, neuron, by = c("region", "ligand"))
lr <- readRDS("data/lr_network_human_21122021.rds")
lt <- readRDS("data/ligand_target_matrix_nsga2r_final.rds")
rows <- list()
for (type in c("GRAN-EN", "VIP-IN")) {
  group <- subset(pb$metadata, cns_region == "Cerebellum" &
                    subtype == type & cells >= 10)
  if (sum(group$condition == "Control") < 2 ||
      sum(group$condition == "MS") < 2) next
  ms <- rownames(group)[group$condition == "MS"]
  z <- pb$counts[, ms, drop = FALSE]
  cpm <- t(t(z) / colSums(z) * 1e6)
  expressed <- rownames(z)[rowSums(cpm > 1) >= 2]
  d <- subset(de, region == "Cerebellum" & subtype == type)
  targets <- intersect(d$gene[d$PValue < 0.05 & d$logFC > 0], expressed)
  targets <- intersect(targets, rownames(lt))
  if (length(targets) < 20) next
  ligands <- unique(sender$ligand[
    sender$region == "Cerebellum" & sender$condition == "MS" &
      sender$donors_detected >= 2 & sender$fraction >= 0.10 &
      sender$fraction >= 2 * sender$neuron_fraction])
  cells <- which(raw$meta$cns_region == "Cerebellum" &
                   raw$meta$cell_subtypes_intg == type & raw$meta$condition == "MS")
  receptors <- intersect(expressed, lr$to)
  fraction <- rowSums(raw$counts[receptors, cells, drop = FALSE] > 0) /
    length(cells)
  receptors <- receptors[fraction >= 0.05]
  ligands <- intersect(ligands, lr$from[lr$to %in% receptors])
  ligands <- intersect(ligands, colnames(lt))
  if (length(ligands) < 5) next
  a <- predict_ligand_activities(targets, intersect(expressed, rownames(lt)),
                                 lt, ligands)
  a$region <- "Cerebellum"
  a$subtype <- type
  a$nominal_targets <- length(targets)
  rows[[length(rows) + 1]] <- a
  message(type, ": ", length(targets), " nominal targets; ",
          length(ligands), " candidate ligands")
}
if (length(rows)) {
  result <- do.call(rbind, rows)
  result <- result[order(result$subtype, -result$aupr_corrected), ]
  write.csv(result, "results/ligand/cerebellum_nominal_activity.csv",
            row.names = FALSE)
  top <- do.call(rbind, lapply(split(result, result$subtype), head, 8))
  p <- ggplot(top, aes(test_ligand, aupr_corrected)) +
    geom_col(fill = "#7b8975") + coord_flip() +
    facet_wrap(~subtype, scales = "free_y") +
    labs(x = NULL, y = "NicheNet corrected AUPR",
         title = "Cerebellar candidates from nominal MS-up genes") +
    theme_bw(base_size = 11)
  ggsave("results/ligand/cerebellum_nominal_ligands.png", p,
         width = 10, height = 6, dpi = 180)
}
