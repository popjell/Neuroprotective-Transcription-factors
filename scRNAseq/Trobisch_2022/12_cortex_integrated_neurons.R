# Test the integrated cortical neuronal labels using the matched raw counts.
# These eight cortical donors overlap Schirmer and are not independent replication.
library(Matrix)
library(edgeR)
library(ggplot2)

x <- readRDS("data/cortex_neurons.rds")
counts <- x$counts
meta <- x$meta
stopifnot(all(colSums(counts) == meta$nCount_RNA))
neurons <- sort(unique(read.delim("data/meta.tsv")$cell_subtypes_intg))
neurons <- neurons[grepl("-EN$|-IN$|^PC$|^motoEN$", neurons)]

key <- paste(meta$cell_subtypes_intg, meta$sample_id, sep = "__")
keys <- sort(unique(key))
membership <- sparseMatrix(i = seq_along(key), j = match(key, keys), x = 1,
                           dims = c(length(key), length(keys)))
pb <- counts %*% membership
colnames(pb) <- keys
info <- meta[match(keys, key), c("sample_id", "condition", "cell_subtypes_intg")]
names(info)[3] <- "subtype"
info$cells <- tabulate(match(key, keys), length(keys))
rownames(info) <- keys
saveRDS(list(counts = pb, metadata = info), "results/cortex_integrated_pseudobulk.rds")
rm(counts, meta, membership, x)
gc()

rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
status <- list()
gene_results <- list()
set_results <- list()
donor_scores <- list()

for (subtype in neurons) {
  use <- info$subtype == subtype
  d <- info[use, , drop = FALSE]
  n_control <- sum(d$condition == "Control")
  n_ms <- sum(d$condition == "MS")
  state <- if (!nrow(d)) "No cells" else if (n_control < 2 || n_ms < 2)
    "Fewer than two donors in one condition" else "Tested"
  status[[subtype]] <- data.frame(
    subtype = subtype, cells = sum(d$cells),
    control_donors = n_control, ms_donors = n_ms,
    smallest_donor_cells = if (nrow(d)) min(d$cells) else 0,
    status = state)
  if (state != "Tested") next

  group <- factor(d$condition, levels = c("Control", "MS"))
  y <- pb[, use, drop = FALSE]
  keep <- filterByExpr(y, group = group)
  if (sum(keep) < 2) {
    status[[subtype]]$status <- "Too few expressed genes"
    next
  }
  result <- tryCatch({
    z <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
    design <- model.matrix(~group)
    z <- estimateDisp(z, design)
    fit <- glmQLFit(z, design)
    de <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
    de$gene <- rownames(de)
    de$subtype <- subtype

    logcpm <- cpm(z, log = TRUE, prior.count = 1)
    for (rank_name in names(sets)) {
      genes <- intersect(sets[[rank_name]], rownames(logcpm))
      if (length(genes) < 3) next
      score <- colMeans(logcpm[genes, , drop = FALSE])
      donor_scores[[paste(subtype, rank_name)]] <- data.frame(
        subtype = subtype, rank = rank_name, sample_id = d$sample_id,
        condition = d$condition, cells = d$cells, score = score)
      model <- lm(score ~ group)
      set_results[[paste(subtype, rank_name)]] <- data.frame(
        subtype = subtype, rank = rank_name, genes_tested = length(genes),
        control_donors = n_control, ms_donors = n_ms,
        ms_minus_control = unname(coef(model)[2]),
        p_value = summary(model)$coefficients[2, 4])
    }
    de
  }, error = function(e) {
    status[[subtype]]$status <<- paste("Model failed:", conditionMessage(e))
    NULL
  })
  if (!is.null(result)) gene_results[[subtype]] <- result
}

status <- do.call(rbind, status)
write.csv(status, "results/cortex_integrated_test_status.csv", row.names = FALSE)
if (length(gene_results)) {
  de <- do.call(rbind, gene_results)
  de$q_global <- p.adjust(de$PValue, method = "BH")
  write.csv(de, "results/cortex_integrated_de.csv", row.names = FALSE)
  rank_de <- merge(rank_genes, de, by.x = "human_symbol", by.y = "gene")
  write.csv(rank_de, "results/cortex_integrated_rank_genes.csv", row.names = FALSE)
  p <- ggplot(de, aes(logFC, -log10(PValue), color = q_global < 0.05)) +
    geom_point(size = 0.35, alpha = 0.6) + facet_wrap(~subtype) +
    scale_color_manual(values = c("grey75", "#C84B31"), na.value = "grey75") +
    theme_minimal() + labs(title = "Trobisch cortical neuronal pseudobulk",
                           subtitle = "Overlapping Schirmer donors; global FDR across all tests",
                           x = "MS-control log2 fold change", y = "-log10(p)",
                           color = "FDR < 0.05")
  ggsave("plots/cortex_integrated_volcano.png", p, width = 13, height = 9, dpi = 180)
}
if (length(set_results)) {
  tests <- do.call(rbind, set_results)
  tests$q_global <- p.adjust(tests$p_value, method = "BH")
  write.csv(tests, "results/cortex_integrated_rank_sets.csv", row.names = FALSE)
  write.csv(do.call(rbind, donor_scores),
            "results/cortex_integrated_donor_scores.csv", row.names = FALSE)
}

bars <- status
bars$significant_genes <- 0L
if (length(gene_results)) {
  hits <- aggregate(as.integer(de$FDR < 0.05) ~ subtype, de, sum)
  bars$significant_genes <- hits[[2]][match(bars$subtype, hits$subtype)]
  bars$significant_genes[is.na(bars$significant_genes)] <- 0L
}
p <- ggplot(bars, aes(subtype, significant_genes, fill = status == "Tested")) +
  geom_col() + coord_flip() + theme_minimal() +
  scale_fill_manual(values = c("grey75", "#3B82F6"), guide = "none") +
  labs(title = "Trobisch cortical integrated neuronal groups",
       subtitle = "Overlapping Schirmer donors; FDR < 0.05 within group",
       x = NULL, y = "MS-control genes")
ggsave("plots/cortex_integrated_de_counts.png", p, width = 8, height = 7, dpi = 180)
