# Check the cortical rank signals against donor age, RNA quality, and PMI.
# Each factor is adjusted separately because there are only eight donors.
library(Matrix)
library(edgeR)
library(limma)

x <- readRDS("data/cortex_neurons.rds")
donors <- unique(x$meta[, c("sample_id", "age", "rin", "pmi")])
stopifnot(!anyDuplicated(donors$sample_id))
qc <- aggregate(cbind(nCount_RNA, nFeature_RNA, mito_percent) ~ sample_id,
                x$meta, median)
donors <- merge(donors, qc, by = "sample_id")
rm(x)
pb <- readRDS("results/cortex_integrated_pseudobulk.rds")
counts <- pb$counts
meta <- pb$metadata
rank_genes <- read.csv("../Schirmer_2019/data/rank_genes.csv")
rank_genes <- unique(rank_genes[rank_genes$use_primary_analysis == "Y",
                                c("rank", "human_symbol")])
sets <- split(rank_genes$human_symbol, rank_genes$rank)
gene_out <- list()
set_out <- list()
camera_out <- list()

for (factor_name in c("age", "rin", "pmi")) {
  for (subtype in sort(unique(meta$subtype))) {
    use <- meta$subtype == subtype
    d <- meta[use, ]
    if (sum(d$condition == "Control") < 2 || sum(d$condition == "MS") < 2) next
    d$condition <- factor(d$condition, levels = c("Control", "MS"))
    d$covariate <- donors[[factor_name]][match(d$sample_id, donors$sample_id)]
    design <- model.matrix(~covariate + condition, d)
    if (qr(design)$rank < ncol(design)) next
    y <- counts[, use, drop = FALSE]
    keep <- filterByExpr(y, group = d$condition)
    if (sum(keep) < 2) next
    z <- normLibSizes(DGEList(y[keep, , drop = FALSE]))
    z <- estimateDisp(z, design)
    fit <- glmQLFit(z, design)
    de <- topTags(glmQLFTest(fit, coef = 3), n = Inf, sort.by = "none")$table
    de$gene <- rownames(de)
    de$subtype <- subtype
    de$adjusted_for <- factor_name
    gene_out[[paste(factor_name, subtype)]] <- de

    logcpm <- cpm(z, log = TRUE, prior.count = 1)
    v <- voom(z, design, plot = FALSE)
    index <- lapply(sets, function(s) which(rownames(v$E) %in% s))
    index <- index[lengths(index) >= 3]
    if (length(index)) {
      cam <- camera(v, index, design, contrast = 3)
      camera_out[[paste(factor_name, subtype)]] <- data.frame(
        adjusted_for = factor_name, subtype = subtype, rank = rownames(cam),
        genes_tested = cam$NGenes, direction = cam$Direction,
        p_value = cam$PValue)
    }
    for (rank_name in names(sets)) {
      genes <- intersect(sets[[rank_name]], rownames(logcpm))
      if (length(genes) < 3) next
      score <- colMeans(logcpm[genes, , drop = FALSE])
      model <- lm(score ~ covariate + condition, d)
      set_out[[paste(factor_name, subtype, rank_name)]] <- data.frame(
        adjusted_for = factor_name, subtype = subtype, rank = rank_name,
        genes_tested = length(genes),
        ms_minus_control = unname(coef(model)[3]),
        p_value = summary(model)$coefficients[3, 4])
    }
  }
}

genes <- do.call(rbind, gene_out)
genes$q_global <- ave(genes$PValue, genes$adjusted_for,
                      FUN = function(p) p.adjust(p, method = "BH"))
rank_hits <- merge(rank_genes, genes, by.x = "human_symbol", by.y = "gene")
write.csv(rank_hits, "results/cortex_rank_gene_sensitivity.csv", row.names = FALSE)

sets_out <- do.call(rbind, set_out)
sets_out$q_global <- ave(sets_out$p_value, sets_out$adjusted_for,
                         FUN = function(p) p.adjust(p, method = "BH"))
write.csv(sets_out, "results/cortex_rank_set_sensitivity.csv", row.names = FALSE)

camera_out <- do.call(rbind, camera_out)
camera_out$q_global <- ave(camera_out$p_value, camera_out$adjusted_for,
                           FUN = function(p) p.adjust(p, method = "BH"))
write.csv(camera_out, "results/cortex_camera_sensitivity.csv", row.names = FALSE)

# Leave one donor out for every primary rank-set result passing FDR.
primary <- read.csv("results/cortex_integrated_rank_sets.csv")
scores <- read.csv("results/cortex_integrated_donor_scores.csv")
stable <- list()
for (i in which(primary$q_global < 0.05)) {
  d <- scores[scores$subtype == primary$subtype[i] &
                scores$rank == primary$rank[i], ]
  effects <- sapply(seq_len(nrow(d)), function(j) {
    z <- d[-j, ]
    unname(coef(lm(score ~ factor(condition, levels = c("Control", "MS")), z))[2])
  })
  stable[[i]] <- data.frame(
    subtype = primary$subtype[i], rank = primary$rank[i],
    donors = nrow(d), same_direction = sum(sign(effects) ==
      sign(primary$ms_minus_control[i])),
    smallest_effect = min(effects), largest_effect = max(effects))
}
write.csv(do.call(rbind, stable), "results/cortex_rank_set_leave_one_out.csv",
          row.names = FALSE)

# Check the broad, subtype-balanced UCell result with the same donor factors.
ucell <- read.csv("results/cortex_integrated_ucell_donors.csv")
coverage <- aggregate(sample_id ~ subtype, ucell,
                      function(z) length(unique(z)))
complete <- coverage$subtype[coverage$sample_id == length(unique(ucell$sample_id))]
balanced <- ucell[ucell$subtype %in% complete, ]
broad <- aggregate(score ~ sample_id + condition + rank, balanced, mean)
broad <- merge(broad, donors, by = "sample_id")
broad_sensitivity <- list()
for (factor_name in c("age", "rin", "pmi", "nCount_RNA",
                      "nFeature_RNA", "mito_percent")) {
  for (rank_name in names(sets)) {
    d <- broad[broad$rank == rank_name, ]
    d$condition <- factor(d$condition, levels = c("Control", "MS"))
    d$covariate <- d[[factor_name]]
    model <- lm(score ~ covariate + condition, d)
    broad_sensitivity[[paste(factor_name, rank_name)]] <- data.frame(
      adjusted_for = factor_name, rank = rank_name,
      ms_minus_control = unname(coef(model)[3]),
      p_value = summary(model)$coefficients[3, 4])
  }
}
broad_sensitivity <- do.call(rbind, broad_sensitivity)
broad_sensitivity$q_two_lists <- ave(broad_sensitivity$p_value,
  broad_sensitivity$adjusted_for, FUN = function(p) p.adjust(p, method = "BH"))
write.csv(broad_sensitivity, "results/cortex_broad_ucell_sensitivity.csv",
          row.names = FALSE)
