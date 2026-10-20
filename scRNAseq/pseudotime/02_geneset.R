# Module for one gene set: drop genes that are not in the data or are detected in too few nuclei, then score the rest in every nucleus
library(Matrix)
library(UCell)
source("common.R")

nuclei <- read.csv(shared("nuclei.csv"))
counts <- readRDS(shared("counts.rds"))
stopifnot(identical(nuclei$nucleus, colnames(counts)))

geneset <- load_geneset()
idx <- match(geneset, toupper(rownames(counts)))
detected <- rep(NA_real_, length(geneset))
detected[!is.na(idx)] <- Matrix::rowMeans(counts[idx[!is.na(idx)], , drop = FALSE] > 0)
genes <- data.frame(gene = geneset, data_symbol = rownames(counts)[idx], in_dataset = !is.na(idx), fraction_nuclei_detected = detected)
genes$kept <- genes$in_dataset & genes$fraction_nuclei_detected >= min_detect
genes$reason <- ifelse(!genes$in_dataset, "not in dataset", ifelse(!genes$kept, paste0("detected in <", 100 * min_detect, "% of nuclei"), ""))
write.csv(genes, res("geneset_genes.csv"), row.names = FALSE)

module <- genes$data_symbol[genes$kept]
cat(geneset_name, ":", length(geneset), "genes,", sum(!genes$in_dataset), "not in the dataset,", sum(genes$in_dataset & !genes$kept), "detected in <", 100 * min_detect, "% of nuclei,", length(module), "kept\n")
stopifnot("fewer than 3 genes left in the module" = length(module) >= 3)

# Module score (UCell on raw counts) and log-normalised expression (counts per 10,000) of the module genes
score <- ScoreSignatures_UCell(counts, features = list(module = module), ncores = 1)[, 1]
expr <- as.matrix(counts[module, , drop = FALSE] %*% Diagonal(x = 1e4 / Matrix::colSums(counts)))
expr <- log1p(expr)

# Ten pseudotime deciles of about the same number of nuclei
pt <- nuclei$author_pseudotime
decile <- cut(pt, quantile(pt, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
write.csv(data.frame(nucleus = nuclei$nucleus, donor = nuclei$donor, condition = nuclei$condition, author_pseudotime = pt, decile = decile, module_ucell = unname(score)),
          res("module_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(1:10, function(d) {
  i <- decile == d
  data.frame(decile = d, n_nuclei = sum(i), mean_pseudotime = mean(pt[i]), mean_module = mean(score[i]), se_module = sd(score[i]) / sqrt(sum(i)))
})), res("decile_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(1:10, function(d) data.frame(gene = module, decile = d, mean_expression = rowMeans(expr[, decile == d, drop = FALSE])))),
          res("decile_gene_means.csv"), row.names = FALSE)
