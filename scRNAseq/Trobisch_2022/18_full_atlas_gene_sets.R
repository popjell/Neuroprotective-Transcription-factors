# Test Rank 1 and Rank 2 across all Trobisch neuronal donors.
# The combined test adjusts for region; regional tests show where signals arise.
library(Matrix)
library(edgeR)
library(limma)
library(fgsea)
library(ggplot2)

cortex <- readRDS("results/cortex_integrated_pseudobulk.rds")
other <- readRDS("results/integrated_neuron_pseudobulk.rds")
stopifnot(identical(rownames(cortex$counts), rownames(other$counts)))
cortex$metadata$cns_region <- "Cortex"
counts <- cbind(cortex$counts, other$counts)
info <- rbind(cortex$metadata, other$metadata)
stopifnot(ncol(counts) == nrow(info),
          !anyDuplicated(c(unique(cortex$metadata$sample_id),
                           unique(other$metadata$sample_id))))
rm(cortex, other)

rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
labels <- sort(unique(read.delim("data/meta.tsv")$cell_subtypes_intg))
labels <- labels[grepl("-EN$|-IN$|^PC$|^motoEN$", labels)]

results <- list()
set.seed(1)
for (scope in c("Combined", "Cortex", "Cerebellum", "Spinal_cord")) {
for (subtype in labels) {
  # Tiny donor-subtype groups do not provide a reliable expression profile.
  all_use <- info$subtype == subtype &
    (scope == "Combined" | info$cns_region == scope)
  use <- all_use & info$cells >= 10
  d <- info[use, ]
  group <- factor(d$condition, levels = c("Control", "MS"))
  regions <- paste(sort(unique(d$cns_region)), collapse = ", ")
  design <- if (!nrow(d)) NULL else if (length(unique(d$cns_region)) > 1)
    model.matrix(~factor(d$cns_region) + group) else model.matrix(~group)
  status <- if (sum(group == "Control") < 2 || sum(group == "MS") < 2)
    "Too few donors" else if (qr(design)$rank < ncol(design) ||
                              nrow(d) - ncol(design) < 2)
    "Disease effect not estimable" else "Tested"

  rows <- lapply(names(sets), function(rank_name) {
    data.frame(scope = scope, subtype = subtype, rank = rank_name,
               method = c("GSEA", "ROAST"),
               total_cells = sum(info$cells[all_use]),
               small_groups = sum(all_use & info$cells < 10),
               control_donors = sum(group == "Control"),
               ms_donors = sum(group == "MS"), regions = regions,
               genes_tested = NA_integer_, direction = NA_character_,
               effect = NA_real_, p_value = NA_real_, status = status)
  })
  names(rows) <- names(sets)
  if (status != "Tested") {
    results[[paste(scope, subtype)]] <- do.call(rbind, rows)
    next
  }

  y <- counts[, use, drop = FALSE]
  # Keep genes expressed in any region with at least four usable donors.
  # A single cross-region filter can discard cortex-expressed genes.
  region_samples <- split(seq_len(nrow(d)), d$cns_region)
  keep <- Reduce(`|`, lapply(region_samples, function(i) {
    if (length(i) < 4) return(rep(FALSE, nrow(y)))
    filterByExpr(y[, i, drop = FALSE], group = group[i])
  }))
  if (sum(keep) < 2) {
    z <- do.call(rbind, rows)
    z$status <- "Too few expressed genes"
    results[[paste(scope, subtype)]] <- z
    next
  }
  y <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  stats <- fit$t[, "groupMS"]
  names(stats) <- rownames(v$E)
  stats <- stats[is.finite(stats)]
  testable <- lapply(sets, intersect, y = names(stats))
  for (rank_name in names(sets)) {
    rows[[rank_name]]$genes_tested <- length(testable[[rank_name]])
    if (length(testable[[rank_name]]) < 3)
      rows[[rank_name]]$status <- "Fewer than three expressed list genes"
  }
  valid <- testable[lengths(testable) >= 3]
  if (length(valid)) {
    fg <- fgseaMultilevel(valid, stats, minSize = 3, maxSize = 500,
                         eps = 0, BPPARAM = BiocParallel::SerialParam())
    roast <- mroast(v, lapply(valid, function(g) which(rownames(v$E) %in% g)),
                    design, contrast = "groupMS", nrot = 9999)
    for (rank_name in names(valid)) {
      f <- fg[fg$pathway == rank_name, ]
      r <- roast[rank_name, ]
      rows[[rank_name]]$direction[1] <- if (f$NES > 0) "MS-up" else "MS-down"
      rows[[rank_name]]$effect[1] <- f$NES
      rows[[rank_name]]$p_value[1] <- f$pval
      rows[[rank_name]]$direction[2] <- paste0("MS-", tolower(r$Direction))
      rows[[rank_name]]$p_value[2] <- r$PValue
    }
  }
  results[[paste(scope, subtype)]] <- do.call(rbind, rows)
  rm(y, v, fit)
  gc(verbose = FALSE)
  message("Finished ", scope, " ", subtype)
}
}

results <- do.call(rbind, results)
rownames(results) <- NULL
combined <- results[results$scope == "Combined", ]
regional <- results[results$scope != "Combined", ]
combined$q_method <- ave(combined$p_value, combined$method,
                         FUN = function(p) p.adjust(p, method = "BH"))
regional$q_method <- ave(regional$p_value, regional$method,
                         FUN = function(p) p.adjust(p, method = "BH"))
regional$q_region <- ave(regional$p_value,
                         interaction(regional$method, regional$scope),
                         FUN = function(p) p.adjust(p, method = "BH"))
write.csv(combined, "results/full_atlas_gene_set_tests.csv", row.names = FALSE)
write.csv(regional, "results/full_atlas_region_gene_set_tests.csv", row.names = FALSE)

# Show significant, nonsignificant, and untestable comparisons separately.
combined$result <- ifelse(is.na(combined$q_method), "Untestable",
                          ifelse(combined$q_method < 0.05, "q < 0.05", "q >= 0.05"))
combined$subtype <- factor(combined$subtype, levels = rev(labels))
p <- ggplot(combined, aes(rank, subtype, fill = result)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = c("q < 0.05" = "#28689C",
                               "q >= 0.05" = "#C7D3DD",
                               "Untestable" = "#F0F0F0")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank()) +
  labs(title = "Trobisch full-atlas neuronal gene-set tests",
       subtitle = "Usable donor-subtype samples; adjusted for region where needed",
       fill = NULL)
ggsave("plots/full_atlas_gene_sets.png", p, width = 9, height = 7, dpi = 180)
