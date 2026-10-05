# Show donor-level ligand detection for the leading L2-L3 neuron candidates.
library(Matrix)
library(ggplot2)

activity <- read.csv("results/ligand/ligand_activity.csv")
support <- read.csv("results/ligand/ligand_receptor_support.csv")
top <- head(subset(activity, region == "Cortex" & subtype == "L2-L3-EN")$test_ligand, 5)
pairs <- unique(subset(support, region == "Cortex" & subtype == "L2-L3-EN" &
                         from %in% top)[, c("from", "sender")])
names(pairs) <- c("ligand", "sender")
raw <- readRDS("data/cortex_senders.rds")
rows <- list()
for (i in seq_len(nrow(pairs))) {
  ligand <- pairs$ligand[i]
  sender <- pairs$sender[i]
  cells <- which(raw$meta$cell_subtypes_intg == sender)
  for (donor in unique(raw$meta$sample_id[cells])) {
    use <- cells[raw$meta$sample_id[cells] == donor]
    if (length(use) < 20) next
    rows[[length(rows) + 1]] <- data.frame(
      ligand = ligand, sender = sender, donor = donor,
      condition = raw$meta$condition[use[1]], cells = length(use),
      fraction = mean(raw$counts[ligand, use] > 0))
  }
}
donors <- do.call(rbind, rows)
donors$condition <- factor(donors$condition, levels = c("Control", "MS"))
write.csv(donors, "results/ligand/top_ligand_donors.csv", row.names = FALSE)

p <- ggplot(donors, aes(condition, fraction, color = condition)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, color = "grey40", fill = "white") +
  geom_jitter(width = 0.08, size = 2) +
  facet_wrap(~ligand + sender, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Control" = "#4b7190", "MS" = "#b75542")) +
  labs(x = NULL, y = "Fraction of nuclei with ligand detected",
       title = "Candidate ligand detection by donor") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave("results/ligand/top_ligand_donors.png", p,
       width = 10, height = 7, dpi = 180)
