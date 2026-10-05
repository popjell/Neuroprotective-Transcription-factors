# Check whether each rank list rises relative to other expressed genes.
library(Matrix)
library(edgeR)
library(limma)

atlas <- read.delim("data/meta.tsv")
neurons <- sort(unique(atlas$cell_subtypes_intg[
  grepl("-EN$|-IN$|^PC$|^motoEN$", atlas$cell_subtypes_intg)]))
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)

files <- c("results/integrated_neuron_pseudobulk.rds",
           "results/cortex_integrated_pseudobulk.rds")
out <- list()
for (file in files[file.exists(files)]) {
  x <- readRDS(file)
  pb <- x$counts
  meta <- x$metadata
  if (!"cns_region" %in% names(meta)) meta$cns_region <- "Cortex"
  for (region in unique(meta$cns_region)) {
    for (subtype in neurons) {
      use <- meta$cns_region == region & meta$subtype == subtype
      d <- meta[use, , drop = FALSE]
      state <- if (sum(d$condition == "Control") < 2 ||
                   sum(d$condition == "MS") < 2) "Too few donors" else "Tested"
      test <- NULL
      if (state == "Tested") {
        group <- factor(d$condition, levels = c("Control", "MS"))
        y <- pb[, use, drop = FALSE]
        keep <- filterByExpr(y, group = group)
        if (sum(keep) < 2) state <- "Too few expressed genes" else {
          test <- tryCatch({
            z <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
            design <- model.matrix(~group)
            v <- voom(z, design, plot = FALSE)
            index <- lapply(sets, function(s) which(rownames(v$E) %in% s))
            index <- index[lengths(index) >= 3]
            if (!length(index)) NULL else {
              result <- camera(v, index, design, contrast = 2)
              data.frame(rank = rownames(result), genes_tested = result$NGenes,
                         direction = result$Direction, p_value = result$PValue)
            }
          }, error = function(e) NULL)
        }
      }
      for (rank_name in names(sets)) {
        row <- if (!is.null(test)) test[test$rank == rank_name, ] else data.frame()
        out[[paste(region, subtype, rank_name)]] <- data.frame(
          region = region, subtype = subtype, rank = rank_name,
          genes_tested = if (nrow(row)) row$genes_tested else NA_integer_,
          direction = if (nrow(row)) row$direction else NA_character_,
          p_value = if (nrow(row)) row$p_value else NA_real_,
          status = if (nrow(row)) "Tested" else if (state == "Tested")
            "Too few list genes" else state)
      }
    }
  }
  rm(x, pb, meta)
  gc()
}

result <- do.call(rbind, out)
result$q_global <- p.adjust(result$p_value, method = "BH")
result$q_new_regions <- NA_real_
new <- result$region != "Cortex"
result$q_new_regions[new] <- p.adjust(result$p_value[new], method = "BH")
write.csv(result, "results/integrated_neuron_camera.csv", row.names = FALSE)
