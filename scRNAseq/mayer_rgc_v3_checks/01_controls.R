# Positive controls (frozen expectations from Mayer et al.): CFH pooled; directional gene-set tests in the resilient and susceptible groups. Reuses the v2 Analysis D model. Read-only on v2 and baseline.
source("../mayer_rgc_v2/common_v2.R")
source("../mayer_rgc_v2/fit_lib.R")
v2_res <- "../mayer_rgc_v2/results"
n_roast <- 9999

# 1a. CFH, pooled MS vs Control (baseline DESeq2 and Analysis B), read from existing results
base <- read.csv(file.path(base_res, "all_RGC", "deseq2.csv"), stringsAsFactors = FALSE)
bt <- read.csv(file.path(v2_res, "B_all_genes.csv"), stringsAsFactors = FALSE)
i1 <- match("CFH", base$gene); i2 <- match("CFH", bt$gene)  # NA if the gene was removed by that model's filter
cfh <- rbind(
  data.frame(model = "baseline pooled DESeq2", log2FC = base$log2FC[i1], se = base$lfcSE[i1], p = base$pvalue[i1], padj = base$padj[i1]),
  data.frame(model = "Analysis B subtype-adjusted limma-voom", log2FC = bt$log2FC[i2], se = bt$se[i2], p = bt$p[i2], padj = bt$padj[i2]))
cfh$present <- c(!is.na(i1), !is.na(i2))
write.csv(cfh, file.path(out_root, "controls_CFH_pooled.csv"), row.names = FALSE)
print(cfh)
for (i in 1:2) register_test(paste0("PC_CFH_", c("pooledDESeq2", "B")[i]), 3, "controls", "CFH", cfh$log2FC[i], cfh$se[i], cfh$p[i], 10, 12, paste0(cfh$model[i], "; expected up; padj ", signif(cfh$padj[i], 3), if (!cfh$present[i]) "; CFH removed by the gene filter" else ""))
log_line("1a CFH: ", paste(paste0(cfh$model, " log2FC ", round(cfh$log2FC, 3), " SE ", round(cfh$se, 3), " p ", signif(cfh$p, 3), " padj ", signif(cfh$padj, 3)), collapse = " | "))

# group model: limma-voom ~ group*condition, block = retina (same code as v2 Analysis D)
pb <- readRDS(file.path(v2_res, "pseudobulk.rds"))
info <- eligible_units(pb$pbg$info, min_cells, min_per_cond)
design <- model.matrix(~ label * condition, info)
colnames(design) <- c("intercept", "label_susceptible", "conditionMS", "interaction")
m <- fit_voom_block(pb$pbg$pb, info, design)
res_tab <- coef_table(m, "conditionMS")
fs <- eBayes(contrasts.fit(m$fit, cbind(sus = c(0, 0, 1, 1))))
tts <- topTable(fs, coef = 1, number = Inf, sort.by = "none")
sus_tab <- data.frame(gene = rownames(tts), log2FC = tts$logFC, se = sqrt(fs$s2.post) * fs$stdev.unscaled[, 1], t = tts$t, p = tts$P.Value)
# the fit must reproduce the v2 Analysis D PSAP estimates
v2r <- read.csv(file.path(v2_res, "D_resilient_MS_effect_all_genes.csv")); v2s <- read.csv(file.path(v2_res, "D_susceptible_MS_effect_all_genes.csv"))
stopifnot(abs(res_tab$log2FC[res_tab$gene == "PSAP"] - v2r$log2FC[v2r$gene == "PSAP"]) < 1e-8, abs(sus_tab$log2FC[sus_tab$gene == "PSAP"] - v2s$log2FC[v2s$gene == "PSAP"]) < 1e-8)
log_line("controls: group model reproduces v2 Analysis D (", m$n_genes, " genes after filterByExpr; duplicateCorrelation ", round(m$cor, 4), ")")

sets <- list(
  resilient_induced_up = list(group = "resilient", expect = "up", genes = c("JUN", "JUND", "FOSB", "ACSL4", "SQSTM1", "NFKBIA", "NFKBIZ", "TGFB2")),
  resilient_identity_down = list(group = "resilient", expect = "down", genes = c("TBR1", "POU4F1", "EBF2", "DLGAP1", "SORCS1", "GRIN2C", "GPT2", "PSAT1")),
  susceptible_up = list(group = "susceptible", expect = "up", genes = c("DCN", "SCUBE1", "BMP5", "MT-ATP8", "MT-ND3", "MT-ND4", "MT-ND4L")))
gene_rows <- list(); set_rows <- list()
for (nm in names(sets)) {
  s <- sets[[nm]]
  tab <- if (s$group == "resilient") res_tab else sus_tab
  cvec <- if (s$group == "resilient") c(0, 0, 1, 0) else c(0, 0, 1, 1)
  in_mat <- s$genes %in% pb$genes
  test <- s$genes %in% rownames(m$v)
  gr <- data.frame(set = nm, group = s$group, expected = s$expect, gene = s$genes, status = ifelse(!in_mat, "missing: not in matrix", ifelse(!test, "missing: removed by filterByExpr", "testable")))
  gr <- cbind(gr, tab[match(s$genes, tab$gene), c("log2FC", "se", "p")])
  gr$matches_expected <- ifelse(gr$status == "testable", ifelse(s$expect == "up", gr$log2FC > 0, gr$log2FC < 0), NA)
  gene_rows[[nm]] <- gr
  idx <- s$genes[test]
  set.seed(1010)
  ro <- roast(m$v, index = idx, design = design, contrast = cvec, block = m$block, correlation = m$cor, nrot = n_roast)$p.value
  fr <- fry(m$v, index = list(set = idx), design = design, contrast = cvec, block = m$block, correlation = m$cor)
  est_mean <- mean(gr$log2FC[gr$status == "testable"]); est_med <- median(gr$log2FC[gr$status == "testable"])
  p_dir <- if (s$expect == "up") ro["Up", "P.Value"] else ro["Down", "P.Value"]
  fry_dir <- if ((fr$Direction == "Up") == (s$expect == "up")) fr$PValue / 2 else 1 - fr$PValue / 2
  set_rows[[nm]] <- data.frame(set = nm, group = s$group, expected = s$expect, n_in_set = length(s$genes), n_testable = length(idx), n_matching_expected = sum(gr$matches_expected, na.rm = TRUE),
                               missing = paste(s$genes[!test], collapse = ","), mean_log2FC = est_mean, median_log2FC = est_med, roast_p_directional = p_dir, roast_p_up = ro["Up", "P.Value"], roast_p_down = ro["Down", "P.Value"],
                               roast_p_two_sided = ro["UpOrDown", "P.Value"], roast_p_mixed = ro["Mixed", "P.Value"], fry_direction = fr$Direction, fry_p_two_sided = fr$PValue, fry_p_directional = fry_dir,
                               direction_reproduced = (est_mean > 0) == (s$expect == "up"))
  register_test(paste0("PC_set_", nm), 3, "controls", nm, est_mean, NA, p_dir, 10, 12, paste0("ROAST directional p (", s$expect, "), ", n_roast, " rotations, ", length(idx), " of ", length(s$genes), " genes testable; estimate = mean log2FC; ", s$group, " group MS effect"))
  log_line("controls ", nm, ": testable ", length(idx), "/", length(s$genes), "; matching expected ", sum(gr$matches_expected, na.rm = TRUE), "; mean log2FC ", round(est_mean, 3), "; ROAST directional p ", signif(p_dir, 3), "; missing: ", ifelse(any(!test), paste(s$genes[!test], collapse = ","), "none"))
}
gt <- do.call(rbind, gene_rows); st <- do.call(rbind, set_rows)
rownames(gt) <- NULL; rownames(st) <- NULL
write.csv(gt, file.path(out_root, "controls_gene_results.csv"), row.names = FALSE)
write.csv(st, file.path(out_root, "controls_set_tests.csv"), row.names = FALSE)
print(st[, c("set", "n_testable", "n_matching_expected", "mean_log2FC", "roast_p_directional", "fry_p_directional", "direction_reproduced")], digits = 3)
