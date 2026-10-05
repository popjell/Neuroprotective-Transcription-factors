# Test every neuronal label from the authors' full integrated atlas.
library(Matrix)
library(edgeR)
library(ggplot2)

atlas <- read.delim("data/meta.tsv")
neurons <- sort(unique(atlas$cell_subtypes_intg[
  grepl("-EN$|-IN$|^PC$|^motoEN$", atlas$cell_subtypes_intg)]))
stopifnot(length(neurons) == 20)

# Show every label in every region, including groups with no cells.
samples <- unique(atlas[, c("sample_id", "condition", "cns_region")])
coverage <- merge(data.frame(subtype = neurons), samples, by = NULL)
observed <- aggregate(cellId ~ cell_subtypes_intg + sample_id, atlas, length)
coverage$cells <- observed$cellId[match(
  paste(coverage$subtype, coverage$sample_id),
  paste(observed$cell_subtypes_intg, observed$sample_id))]
coverage$cells[is.na(coverage$cells)] <- 0
coverage$raw_counts_available <- coverage$cns_region != "Cortex"
write.csv(coverage, "results/integrated_neuron_coverage.csv", row.names = FALSE)

# The cached raw-count objects cover the new cerebellum and spinal-cord donors.
x <- readRDS("data/new_regions.rds")
stopifnot(all(colSums(x$counts) == x$meta$nCount_RNA))
stopifnot(all(x$counts@x == round(x$counts@x)))
meta <- x$meta[x$meta$cell_subtypes_intg %in% neurons, ]
counts <- x$counts[, rownames(meta), drop = FALSE]
rm(x)
gc()

# Sum raw counts once per donor and integrated neuronal subtype.
key <- paste(meta$cns_region, meta$cell_subtypes_intg, meta$sample_id, sep = "__")
keys <- sort(unique(key))
membership <- sparseMatrix(i = seq_along(key), j = match(key, keys), x = 1,
                           dims = c(length(key), length(keys)))
pb <- counts %*% membership
colnames(pb) <- keys
info <- meta[match(keys, key), c("sample_id", "condition", "cns_region",
                                  "cell_subtypes_intg")]
names(info)[4] <- "subtype"
info$cells <- tabulate(match(key, keys), length(keys))
rownames(info) <- keys
saveRDS(list(counts = pb, metadata = info), "results/integrated_neuron_pseudobulk.rds")
rm(counts, meta, membership)
gc()

rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
gene_results <- list()
set_results <- list()
eligibility <- list()
donor_scores <- list()

for (region in c("Cerebellum", "Spinal_cord")) {
  for (subtype in neurons) {
    use <- info$cns_region == region & info$subtype == subtype
    d <- info[use, , drop = FALSE]
    n_control <- sum(d$condition == "Control")
    n_ms <- sum(d$condition == "MS")
    status <- if (!nrow(d)) "No cells in this region" else if (
      n_control < 2 || n_ms < 2) "Fewer than two donors in one condition" else "Tested"
    eligibility[[paste(region, subtype)]] <- data.frame(
      region = region, subtype = subtype, cells = sum(d$cells),
      control_donors = n_control, ms_donors = n_ms,
      smallest_donor_cells = if (nrow(d)) min(d$cells) else 0,
      status = status)
    if (status != "Tested") next

    group <- factor(d$condition, levels = c("Control", "MS"))
    y <- pb[, use, drop = FALSE]
    keep <- filterByExpr(y, group = group)
    if (sum(keep) < 2) {
      eligibility[[paste(region, subtype)]]$status <- "Too few expressed genes"
      next
    }

    # edgeR tests MS versus control using donors, not individual nuclei.
    result <- tryCatch({
      z <- calcNormFactors(DGEList(y[keep, , drop = FALSE]))
      design <- model.matrix(~group)
      z <- estimateDisp(z, design)
      fit <- glmQLFit(z, design)
      de <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
      de$gene <- rownames(de)
      de$region <- region
      de$subtype <- subtype

      # A donor-level score is the mean normalized expression of detected list genes.
      logcpm <- cpm(z, log = TRUE, prior.count = 1)
      for (rank_name in names(sets)) {
        genes <- intersect(sets[[rank_name]], rownames(logcpm))
        if (length(genes) < 3) next
        score <- colMeans(logcpm[genes, , drop = FALSE])
        donor_scores[[paste(region, subtype, rank_name)]] <- data.frame(
          region = region, subtype = subtype, rank = rank_name,
          sample_id = d$sample_id, condition = d$condition,
          cells = d$cells, score = score)
        model <- lm(score ~ group)
        set_results[[paste(region, subtype, rank_name)]] <- data.frame(
          region = region, subtype = subtype, rank = rank_name,
          genes_tested = length(genes), control_donors = n_control,
          ms_donors = n_ms, ms_minus_control = unname(coef(model)[2]),
          p_value = summary(model)$coefficients[2, 4])
      }
      de
    }, error = function(e) {
      eligibility[[paste(region, subtype)]]$status <<- paste("Model failed:", conditionMessage(e))
      NULL
    })
    if (!is.null(result)) gene_results[[paste(region, subtype)]] <- result
  }
}

status <- do.call(rbind, eligibility)
status <- rbind(status, data.frame(
  region = "Cortex", subtype = neurons,
  cells = vapply(neurons, function(s) sum(coverage$cells[
    coverage$cns_region == "Cortex" & coverage$subtype == s]), numeric(1)),
  control_donors = vapply(neurons, function(s) sum(coverage$cells[
    coverage$cns_region == "Cortex" & coverage$subtype == s &
      coverage$condition == "Control"] > 0), integer(1)),
  ms_donors = vapply(neurons, function(s) sum(coverage$cells[
    coverage$cns_region == "Cortex" & coverage$subtype == s &
      coverage$condition == "MS"] > 0), integer(1)),
  smallest_donor_cells = NA_real_, status = "Raw counts not yet matched"))
write.csv(status, "results/integrated_neuron_test_status.csv", row.names = FALSE)

if (length(gene_results)) {
  de <- do.call(rbind, gene_results)
  de$q_global <- p.adjust(de$PValue, method = "BH")
  write.csv(de, "results/integrated_neuron_de.csv", row.names = FALSE)
  rank_de <- merge(rank_genes, de, by.x = "human_symbol", by.y = "gene")
  write.csv(rank_de, "results/integrated_neuron_rank_genes.csv", row.names = FALSE)
  p <- ggplot(de, aes(logFC, -log10(PValue), color = q_global < 0.05)) +
    geom_point(size = 0.35, alpha = 0.6) + facet_wrap(~region + subtype) +
    scale_color_manual(values = c("grey75", "#C84B31"), na.value = "grey75") +
    theme_minimal() + labs(title = "Trobisch new-region neuronal pseudobulk",
                           subtitle = "MS versus control; global FDR across all tests",
                           x = "MS-control log2 fold change", y = "-log10(p)",
                           color = "FDR < 0.05")
  ggsave("plots/integrated_neuron_volcano.png", p, width = 13, height = 8, dpi = 180)
}
if (length(set_results)) {
  sets_out <- do.call(rbind, set_results)
  sets_out$q_global <- p.adjust(sets_out$p_value, method = "BH")
  write.csv(sets_out, "results/integrated_neuron_rank_sets.csv", row.names = FALSE)
  scores <- do.call(rbind, donor_scores)
  write.csv(scores, "results/integrated_neuron_donor_scores.csv", row.names = FALSE)
  p <- ggplot(scores, aes(condition, score, color = condition)) +
    geom_point(position = position_jitter(width = 0.1), size = 1.8) +
    facet_grid(rank ~ subtype, scales = "free_y") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Trobisch integrated neuronal groups", subtitle = "One point per donor", x = NULL,
         y = "Mean normalized rank-gene expression")
  ggsave("plots/integrated_neuron_rank_scores.png", p, width = 15, height = 6, dpi = 180)
}

bars <- status[status$region != "Cortex", c("region", "subtype", "status")]
bars$significant_genes <- 0L
if (length(gene_results)) {
  hits <- aggregate(as.integer(de$FDR < 0.05) ~ region + subtype, de, sum)
  bars$significant_genes <- hits[[3]][match(paste(bars$region, bars$subtype),
                                          paste(hits$region, hits$subtype))]
  bars$significant_genes[is.na(bars$significant_genes)] <- 0L
}
p <- ggplot(bars, aes(subtype, significant_genes, fill = status == "Tested")) +
  geom_col() + coord_flip() + facet_wrap(~region, scales = "free_y") +
  scale_fill_manual(values = c("grey75", "#3B82F6"), guide = "none") +
  theme_minimal() + labs(title = "MS-control genes by integrated neuronal group",
                         subtitle = "FDR < 0.05 within group; grey means not estimable",
                         x = NULL, y = "Genes")
ggsave("plots/integrated_neuron_de_counts.png", p, width = 10, height = 9, dpi = 180)
