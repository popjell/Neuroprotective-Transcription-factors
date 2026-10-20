# Analysis G (Tier 3): method concordance for PSAP, ITM2C, CNDP2 on the pooled retina pseudobulk (baseline DESeq2, edgeR QL, limma-trend) plus Analyses B and C.
source("common_v2.R")
suppressPackageStartupMessages({ library(edgeR); library(limma) })

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
cts <- pb$pb_retina
cond <- factor(pb$retina$condition[match(colnames(cts), pb$retina$retina)], levels = c("Control", "MS"))
keep <- rowSums(cts >= 10) >= 2                      # same filter as the baseline DESeq2
cts <- cts[keep, ]
gs <- c("PSAP", "ITM2C", "CNDP2")
design <- model.matrix(~ cond)
nc <- sum(cond == "Control"); nm <- sum(cond == "MS")

# baseline DESeq2, read-only
base <- read.csv(file.path(base_res, "all_RGC", "genes.csv"), stringsAsFactors = FALSE)
rows <- list(data.frame(method = "baseline DESeq2 (pooled)", gene = gs, log2FC = base$log2FC[match(gs, base$gene)], se = base$lfcSE[match(gs, base$gene)], p = base$pvalue[match(gs, base$gene)]))

# edgeR quasi-likelihood; SE derived from |logFC| / sqrt(F) for the one-df contrast
y <- calcNormFactors(DGEList(cts))
y <- estimateDisp(y, design)
qf <- glmQLFit(y, design)
q <- glmQLFTest(qf, coef = 2)$table
rows[[2]] <- data.frame(method = "edgeR quasi-likelihood (pooled)", gene = gs, log2FC = q[gs, "logFC"], se = abs(q[gs, "logFC"]) / sqrt(q[gs, "F"]), p = q[gs, "PValue"])

# limma-trend on log2 CPM (prior 1)
lc <- cpm(y, log = TRUE, prior.count = 1)
ft <- eBayes(lmFit(lc, design), trend = TRUE)
tt <- topTable(ft, coef = 2, number = Inf, sort.by = "none")
rows[[3]] <- data.frame(method = "limma-trend (pooled)", gene = gs, log2FC = tt[gs, "logFC"], se = sqrt(ft$s2.post[gs]) * ft$stdev.unscaled[gs, 2], p = tt[gs, "P.Value"])

for (i in 2:3) for (g in gs) {
  r <- rows[[i]][rows[[i]]$gene == g, ]
  register_test(paste0("G_", c("edgeRQL", "limmatrend")[i - 1], "_", g), 3, "G", g, r$log2FC, r$se, r$p, nc, nm, paste(r$method, "on pooled retina pseudobulk"))
}
# Analyses B and C
for (a in c("B", "C")) {
  tb <- read.csv(file.path(out_root, paste0(a, "_all_genes.csv")), stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <- data.frame(method = c(B = "Analysis B subtype-adjusted limma-voom", C = "Analysis C composition-standardized limma")[a], gene = gs, log2FC = tb$log2FC[match(gs, tb$gene)], se = tb$se[match(gs, tb$gene)], p = tb$p[match(gs, tb$gene)])
}
res <- do.call(rbind, rows)
res$ci_low <- res$log2FC - 1.96 * res$se; res$ci_high <- res$log2FC + 1.96 * res$se
write.csv(res, file.path(out_root, "G_concordance.csv"), row.names = FALSE)
print(res, digits = 3)
log_line("G: concordance written for ", paste(gs, collapse = ","), " across ", length(unique(res$method)), " methods")
