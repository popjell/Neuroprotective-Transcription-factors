# Check every eligible ligand against MS-up Rank 1/2 genes directly.
library(Matrix)

activity <- read.csv("results/ligand/ligand_activity.csv")
de_cortex <- read.csv("results/cortex_integrated_de.csv")
de_cortex$region <- "Cortex"
de_new <- read.csv("results/integrated_neuron_de.csv")
de <- rbind(de_cortex[, names(de_new)], de_new)
rank <- subset(read.csv("../Schirmer_2019/data/rank_genes.csv"),
               use_primary_analysis == "Y")
rank_label <- tapply(rank$rank, rank$human_symbol, paste, collapse = ", ")
lt <- readRDS("data/ligand_target_matrix_nsga2r_final.rds")
rows <- list()
groups <- unique(activity[, c("region", "subtype")])
for (i in seq_len(nrow(groups))) {
  group <- groups[i, ]
  a <- subset(activity, region == group$region & subtype == group$subtype)
  d <- subset(de, region == group$region & subtype == group$subtype &
                FDR < 0.05 & logFC > 0 & gene %in% rank$human_symbol)
  if (nrow(d) == 0) next
  for (target in intersect(d$gene, rownames(lt))) {
    ligands <- intersect(a$test_ligand, colnames(lt))
    weight <- as.numeric(lt[target, ligands])
    target_rank <- colSums(sweep(lt[, ligands, drop = FALSE],
                                  2, weight, ">")) + 1
    rows[[length(rows) + 1]] <- data.frame(
      region = group$region, subtype = group$subtype,
      target = target, rank = unname(rank_label[target]),
      target_logFC = d$logFC[match(target, d$gene)],
      target_FDR = d$FDR[match(target, d$gene)],
      ligand = ligands, ligand_activity = a$aupr_corrected[
        match(ligands, a$test_ligand)],
      target_weight = weight, target_rank_in_ligand = target_rank)
  }
}
if (length(rows))
  write.csv(do.call(rbind, rows),
            "results/ligand/significant_rank_target_links.csv", row.names = FALSE)
