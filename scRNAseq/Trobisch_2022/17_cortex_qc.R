# Plot donor-level nucleus quality measures for the reused cortex.
library(ggplot2)

meta <- readRDS("data/cortex_neurons.rds")$meta
qc <- aggregate(cbind(nCount_RNA, nFeature_RNA, mito_percent) ~
                  sample_id + condition, meta, median)
write.csv(qc, "results/cortex_donor_qc.csv", row.names = FALSE)
long <- reshape(qc, varying = c("nCount_RNA", "nFeature_RNA", "mito_percent"),
                v.names = "value", timevar = "measure",
                times = c("UMIs", "Genes", "Mito percent"), direction = "long")
p <- ggplot(long, aes(condition, value, color = condition)) +
  geom_point(position = position_jitter(width = 0.08), size = 3) +
  facet_wrap(~measure, scales = "free_y") + theme_minimal() +
  labs(title = "Cortical donor-level nucleus quality",
       subtitle = "Median per nucleus; four control and four MS donors",
       x = NULL, y = "Median value")
ggsave("plots/cortex_donor_qc.png", p, width = 9, height = 4, dpi = 180)
