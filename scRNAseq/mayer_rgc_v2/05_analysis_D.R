# Analysis D: resilient vs susceptible groups (author-defined, 27-label). limma-voom ~group*condition with retina as block (interaction is Tier 2, one-sided), per-group DESeq2, and sensitivities.
source("common_v2.R")
source("fit_lib.R")
suppressPackageStartupMessages(library(DESeq2))

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
genes <- pb$genes
info_g <- pb$pbg$info
n_by_cond <- function(info) c(control = length(unique(info$retina[info$condition == "Control"])), ms = length(unique(info$retina[info$condition == "MS"])))

# eligible units: >= 20 cells; a group is fitted only if >= 3 control and >= 3 MS retinas remain
info <- eligible_units(info_g, min_cells, min_per_cond)
log_line("D: groups fitted ", paste(sort(unique(info$label)), collapse = ","), "; units ", nrow(info), "; retinas ", paste(n_by_cond(info), collapse = " Control / "), " MS")
stopifnot(setequal(unique(info$label), c("resilient", "susceptible")))
nb <- n_by_cond(info)
for (g in c("resilient", "susceptible")) log_line("D: ", g, " units ", sum(info$label == g), " (", sum(info$label == g & info$condition == "Control"), " Control, ", sum(info$label == g & info$condition == "MS"), " MS)")

# limma-voom, reference group = resilient (alphabetical), so conditionMS = MS effect in resilient; interaction = susceptible effect minus resilient effect
design <- model.matrix(~ label * condition, info)
colnames(design) <- c("intercept", "label_susceptible", "conditionMS", "interaction")
m <- fit_voom_block(pb$pbg$pb, info, design)
log_line("D: genes after filterByExpr ", m$n_genes, "; duplicateCorrelation ", round(m$cor, 4))
res_tab <- coef_table(m, "conditionMS")
int_tab <- coef_table(m, "interaction")
cm <- cbind(susceptible_MS_effect = c(0, 0, 1, 1))
fs <- eBayes(contrasts.fit(m$fit, cm))
tts <- topTable(fs, coef = 1, number = Inf, sort.by = "none")
sus_tab <- data.frame(gene = rownames(tts), log2FC = tts$logFC, se = sqrt(fs$s2.post) * fs$stdev.unscaled[, 1], t = tts$t, p = tts$P.Value, padj = tts$adj.P.Val, row.names = NULL)
write.csv(res_tab, file.path(out_root, "D_resilient_MS_effect_all_genes.csv"), row.names = FALSE)
write.csv(sus_tab, file.path(out_root, "D_susceptible_MS_effect_all_genes.csv"), row.names = FALSE)
write.csv(int_tab, file.path(out_root, "D_interaction_all_genes.csv"), row.names = FALSE)

g <- function(tab, gene) tab[tab$gene == gene, ]
# interaction hypothesis: MS effect larger in resilient than susceptible, i.e. interaction (susceptible - resilient) < 0
one_sided_neg <- function(r) ifelse(r$t < 0, r$p / 2, 1 - r$p / 2)
ri <- g(int_tab, "PSAP")
register_test("D_PSAP_interaction", 2, "D", "PSAP", ri$log2FC, ri$se, one_sided_neg(ri), nb["control"], nb["ms"], paste0("interaction (susceptible minus resilient MS effect); ONE-SIDED p for resilient > susceptible (interaction < 0); two-sided p = ", signif(ri$p, 4)))
rr <- g(res_tab, "PSAP")
register_test("D_PSAP_resilient", 2, "D", "PSAP", rr$log2FC, rr$se, rr$p, nb["control"], nb["ms"], "MS effect within resilient group, limma-voom ~group*condition, two-sided")
rs <- g(sus_tab, "PSAP")
register_test("D_PSAP_susceptible", 3, "D", "PSAP", rs$log2FC, rs$se, rs$p, nb["control"], nb["ms"], "MS effect within susceptible group, same model, two-sided")
for (gn in "ITM2C") {
  a <- g(int_tab, gn); b <- g(res_tab, gn); cc <- g(sus_tab, gn)
  register_test(paste0("D_", gn, "_interaction"), 3, "D", gn, a$log2FC, a$se, a$p, nb["control"], nb["ms"], "interaction, two-sided")
  register_test(paste0("D_", gn, "_resilient"), 3, "D", gn, b$log2FC, b$se, b$p, nb["control"], nb["ms"], "MS effect within resilient")
  register_test(paste0("D_", gn, "_susceptible"), 3, "D", gn, cc$log2FC, cc$se, cc$p, nb["control"], nb["ms"], "MS effect within susceptible")
}
log_line("D: PSAP resilient MS effect ", round(rr$log2FC, 4), " SE ", round(rr$se, 4), " p ", signif(rr$p, 4), "; susceptible ", round(rs$log2FC, 4), " p ", signif(rs$p, 4),
         "; interaction ", round(ri$log2FC, 4), " SE ", round(ri$se, 4), " one-sided p ", signif(one_sided_neg(ri), 4))

# DESeq2 per unit set, as the authors did (design ~condition, counts >= 10 in >= 2 pseudobulks)
fit_deseq <- function(counts, inf) {
  cd <- data.frame(row.names = inf$unit, condition = factor(inf$condition, levels = c("Control", "MS")))
  keep <- rowSums(counts[, inf$unit] >= 10) >= 2
  dds <- DESeq(DESeqDataSetFromMatrix(counts[keep, inf$unit], cd, ~condition), quiet = TRUE)
  r <- results(dds, contrast = c("condition", "MS", "Control"))
  data.frame(gene = rownames(r), log2FC = r$log2FoldChange, se = r$lfcSE, p = r$pvalue, padj = r$padj, row.names = NULL)
}
ds <- list()
run_deseq <- function(id, label, counts, inf) {
  inf <- inf[inf$n_cells >= min_cells, ]
  nn <- n_by_cond(inf)
  if (nn["control"] < min_per_cond || nn["ms"] < min_per_cond) { log_line("D: skipped ", id, " (", nn["control"], " Control, ", nn["ms"], " MS retinas with >= ", min_cells, " cells)"); return(invisible(NULL)) }
  tb <- fit_deseq(counts, inf)
  r <- tb[tb$gene == "PSAP", ]
  register_test(id, 3, "D", "PSAP", r$log2FC, r$se, r$p, nn["control"], nn["ms"], paste0("DESeq2 ~condition, ", label))
  ds[[id]] <<- data.frame(analysis = id, description = label, n_control = nn["control"], n_ms = nn["ms"], psap_log2FC = r$log2FC, psap_se = r$se, psap_p = r$p, psap_padj = r$padj)
  write.csv(tb, file.path(out_root, paste0(id, "_all_genes.csv")), row.names = FALSE)
  log_line("D: ", id, " PSAP log2FC ", round(r$log2FC, 4), " SE ", round(r$se, 4), " p ", signif(r$p, 4), " (", nn["control"], " vs ", nn["ms"], ")")
}
for (gr in c("resilient", "susceptible")) run_deseq(paste0("D_deseq2_", gr), paste(gr, "group, all member subtypes"), pb$pbg$pb, info_g[info_g$label == gr, ])
# sensitivities (Tier 3): resilient without RGC13 and RGC26 = PG_ON only; each susceptible member subtype alone
i27 <- pb$pb27$info
for (s in c("PG_ON", "MG_ON4", "MG_OFF2", "ipRGC")) run_deseq(paste0("D_sens_deseq2_", s), paste(s, "alone (27-label units)"), pb$pb27$pb, i27[i27$label == s, ])
write.csv(do.call(rbind, ds), file.path(out_root, "D_deseq2_summary.csv"), row.names = FALSE)
