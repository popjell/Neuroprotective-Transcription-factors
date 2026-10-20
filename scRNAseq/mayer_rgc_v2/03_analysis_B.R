# Analysis B: subtype-adjusted limma-voom model (retina x 22-label units, block = retina). Tier 1 = PSAP; Tier 2 = ITM2C and ROAST on the 31 genes. Sensitivities are Tier 3.
source("common_v2.R")
source("fit_lib.R")

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
genes <- pb$genes
info_all <- pb$pb22$info
counts <- pb$pb22$pb
g31 <- read_gene_list(genes)
n_roast <- 9999  # ROAST rotations

run_B <- function(info, quality_weights = TRUE) {
  info <- eligible_units(info, min_cells, min_per_cond)
  design <- model.matrix(~ label + condition, info)
  m <- fit_voom_block(counts, info, design, quality_weights)
  m$info <- info
  m$tab <- coef_table(m, "conditionMS")
  m
}
n_by_cond <- function(info) c(control = length(unique(info$retina[info$condition == "Control"])), ms = length(unique(info$retina[info$condition == "MS"])))

# main fit
m <- run_B(info_all)
info <- m$info
nb <- n_by_cond(info)
log_line("B: units ", nrow(info), "; subtypes kept ", length(unique(info$label)), " (", paste(sort(unique(info$label)), collapse = ","), "); dropped ",
         paste(setdiff(unique(info_all$label), unique(info$label)), collapse = ","))
log_line("B: retinas with >= 1 unit ", nb["control"], " Control, ", nb["ms"], " MS; genes after filterByExpr ", m$n_genes, "; duplicateCorrelation ", round(m$cor, 4))
write.csv(m$tab, file.path(out_root, "B_all_genes.csv"), row.names = FALSE)
write.csv(info, file.path(out_root, "B_units_used.csv"), row.names = FALSE)
saveRDS(list(genes = m$tab$gene, info = info), file.path(out_root, "B_objects.rds"))

p_of <- function(tab, g) tab[tab$gene == g, ]
one_sided <- function(r) ifelse(r$log2FC > 0, r$p / 2, 1 - r$p / 2)

# Tier 1: PSAP; Tier 2: ITM2C
r <- p_of(m$tab, "PSAP")
register_test("B_PSAP", 1, "B", "PSAP", r$log2FC, r$se, r$p, nb["control"], nb["ms"], paste0("subtype-adjusted limma-voom; two-sided p; one-sided (MS>Control) p = ", signif(one_sided(r), 4)))
r2 <- p_of(m$tab, "ITM2C")
register_test("B_ITM2C", 2, "B", "ITM2C", r2$log2FC, r2$se, r2$p, nb["control"], nb["ms"], "subtype-adjusted limma-voom; two-sided")
log_line("B: PSAP log2FC ", round(r$log2FC, 4), " SE ", round(r$se, 4), " p ", signif(r$p, 4), " padj ", signif(r$padj, 4), "; ITM2C log2FC ", round(r2$log2FC, 4), " p ", signif(r2$p, 4))

# ROAST on the 31-gene set (and without PSAP), restricted to genes passing the filter
roast_set <- function(idx, label, tier) {
  idx <- intersect(idx, rownames(m$v))
  set.seed(1010)
  ci <- which(colnames(m$design) == "conditionMS")
  res <- tryCatch(roast(m$v, index = idx, design = m$design, contrast = ci, block = m$block, correlation = m$cor, nrot = n_roast)$p.value, error = function(e) { log_line("ROAST failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    out <- data.frame(set = label, n_genes = length(idx), method = "roast", p_up = res["Up", "P.Value"], p_down = res["Down", "P.Value"], p_two_sided = res["UpOrDown", "P.Value"], p_mixed = res["Mixed", "P.Value"])
    register_test(paste0("B_roast_", label), tier, "B", label, NA, NA, out$p_mixed, nb["control"], nb["ms"], paste0("ROAST mixed-direction p, ", n_roast, " rotations, ", length(idx), " genes; up p ", signif(out$p_up, 3)))
  } else {
    fr <- fry(m$v, index = idx, design = m$design, contrast = ci)
    out <- data.frame(set = label, n_genes = length(idx), method = "fry", p_up = NA, p_down = NA, p_two_sided = fr$PValue, p_mixed = fr$PValue.Mixed)
    register_test(paste0("B_roast_", label), tier, "B", label, NA, NA, out$p_mixed, nb["control"], nb["ms"], "ROAST failed; fry mixed p used")
  }
  out
}
ro <- rbind(roast_set(g31, "set31", 2), roast_set(setdiff(g31, "PSAP"), "set31_no_PSAP", 3))
write.csv(ro, file.path(out_root, "B_roast.csv"), row.names = FALSE)
print(ro)

# heterogeneity: subtype x condition interaction F-test for PSAP, per-subtype estimates (descriptive)
di <- model.matrix(~ label * condition, info)
fi <- eBayes(lmFit(m$v["PSAP", ], di, block = m$block, correlation = m$cor))
ic <- grep(":conditionMS", colnames(di))
ft <- topTable(fi, coef = ic, number = 1)
register_test("B_PSAP_interaction_F", 3, "B", "PSAP", NA, NA, ft$P.Value, nb["control"], nb["ms"], paste0("subtype x condition interaction F-test, ", length(ic), " df"))
per <- do.call(rbind, lapply(sort(unique(info$label)), function(s) {
  u <- which(info$label == s)
  gi <- which(rownames(m$v$E) == "PSAP")  # v$weights has no dimnames
  d <- data.frame(E = m$v$E[gi, u], w = m$v$weights[gi, u], condition = info$condition[u])
  f <- summary(lm(E ~ condition, d, weights = w))
  est <- f$coefficients["conditionMS", ]; tq <- qt(0.975, f$df[2])
  data.frame(subtype = s, n_control = sum(d$condition == "Control"), n_ms = sum(d$condition == "MS"), log2FC = est[1], se = est[2], ci_low = est[1] - tq * est[2], ci_high = est[1] + tq * est[2])
}))
write.csv(per, file.path(out_root, "B_psap_per_subtype.csv"), row.names = FALSE)
log_line("B: PSAP subtype x condition interaction F p ", signif(ft$P.Value, 3))

# sensitivities (each one registry row); long (about 25 refits), so only with B_SENS=1
if (Sys.getenv("B_SENS", "0") == "1") {
sens <- list()
add_sens <- function(id, label, mm) {
  rr <- p_of(mm$tab, "PSAP"); nn <- n_by_cond(mm$info)
  register_test(id, 3, "B", "PSAP", rr$log2FC, rr$se, rr$p, nn["control"], nn["ms"], label)
  sens[[id]] <<- data.frame(sensitivity = id, description = label, n_units = nrow(mm$info), n_subtypes = length(unique(mm$info$label)), n_control = nn["control"], n_ms = nn["ms"], log2FC = rr$log2FC, se = rr$se, p = rr$p, padj = rr$padj)
  log_line("B sens ", id, ": log2FC ", round(rr$log2FC, 4), " SE ", round(rr$se, 4), " p ", signif(rr$p, 4))
}
add_sens("B_sens_a_no_quality_weights", "without sample quality weights (voom only)", run_B(info_all, quality_weights = FALSE))
excl <- c("RGC4", "RGC16", "RGC23", "RGC26", "RGC13")
add_sens("B_sens_b_excl_odd_clusters", "excluding units from RGC4, RGC16, RGC23, RGC26, RGC13", run_B(info_all[!info_all$label %in% excl, ]))
small <- pb$retina$retina[pb$retina$n_cells < 1000]
log_line("B sens c: retinas with < 1000 RGCs: ", paste(small, collapse = ","))
add_sens("B_sens_c_excl_small_retinas", paste("excluding retinas with < 1000 RGCs:", paste(small, collapse = ",")), run_B(info_all[!info_all$retina %in% small, ]))
loo <- do.call(rbind, lapply(pb$retina$retina, function(rt) {
  mm <- run_B(info_all[info_all$retina != rt, ])
  rr <- p_of(mm$tab, "PSAP")
  data.frame(left_out = rt, condition = pb$retina$condition[pb$retina$retina == rt], log2FC = rr$log2FC, se = rr$se, p = rr$p)
}))
write.csv(loo, file.path(out_root, "B_leave_one_retina_out.csv"), row.names = FALSE)
register_test("B_sens_d_leave_one_retina_out", 3, "B", "PSAP", median(loo$log2FC), NA, NA, 10, 12, paste0("median of 22 leave-one-retina-out PSAP log2FC; min ", round(min(loo$log2FC), 4), " (left out ", loo$left_out[which.min(loo$log2FC)], "), max ", round(max(loo$log2FC), 4), " (left out ", loo$left_out[which.max(loo$log2FC)], "); max p ", signif(max(loo$p), 3)))
sens[["loo"]] <- data.frame(sensitivity = "B_sens_d_leave_one_retina_out", description = "22 refits", n_units = NA, n_subtypes = NA, n_control = NA, n_ms = NA, log2FC = median(loo$log2FC), se = NA, p = NA, padj = NA)
write.csv(do.call(rbind, sens), file.path(out_root, "B_sensitivity.csv"), row.names = FALSE)
log_line("B: leave-one-out PSAP log2FC range ", round(min(loo$log2FC), 4), " to ", round(max(loo$log2FC), 4), "; p range ", signif(min(loo$p), 3), " to ", signif(max(loo$p), 3))
}
