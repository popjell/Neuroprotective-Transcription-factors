# Ambient diagnostic: rod, glia and microglia marker levels in RGC nuclei per retina, by condition; rod score vs pooled PSAP; search of deposited text for ambient-correction statements.
source("../mayer_rgc_v2/common_v2.R")

v2_res <- "../mayer_rgc_v2/results"
cnt <- load_counts_ro()
x <- cnt$x; genes <- cnt$genes
meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
stopifnot(identical(cnt$cells, meta$cell_barcode))
retina <- paste0("S", meta$Sample)
cond_of <- tapply(ifelse(meta$Condition == "Ctrl", "Control", "MS"), retina, unique)

panel <- c(rod = "RHO", rod = "PDE6A", rod = "SAG", rod = "GNAT1", glia = "GFAP", glia = "GLUL", glia = "RLBP1", microglia = "P2RY12", microglia = "CSF1R")
present <- panel[panel %in% genes]
log_line("ambient: panel genes present ", length(present), " of ", length(panel), "; missing ", ifelse(all(panel %in% genes), "none", paste(setdiff(panel, genes), collapse = ",")))
cn <- as.matrix(x[match(present, genes), ])
rownames(cn) <- present
ln <- log1p(t(t(cn) / colSums(x) * 1e4))  # log1p counts per 10,000

# retina-level mean log-CPM and detection fraction
rows <- do.call(rbind, lapply(present, function(g) {
  data.frame(gene = g, class = names(panel)[panel == g], retina = names(cond_of), condition = as.character(cond_of),
             mean_lognorm = as.numeric(tapply(ln[g, ], retina, mean)[names(cond_of)]), detect_frac = as.numeric(tapply(cn[g, ] > 0, retina, mean)[names(cond_of)]))
}))
rod <- present[names(present) == "rod"]
rodscore <- data.frame(gene = "ROD_SCORE", class = "rod", retina = names(cond_of), condition = as.character(cond_of),
                       mean_lognorm = as.numeric(tapply(colMeans(ln[rod, , drop = FALSE]), retina, mean)[names(cond_of)]), detect_frac = NA)
rows <- rbind(rows, rodscore)
write.csv(rows, file.path(out_root, "ambient_retina_values.csv"), row.names = FALSE)

# by condition, exact permutation over retinas
tests <- list()
for (g in unique(rows$gene)) for (met in c("mean_lognorm", "detect_frac")) {
  if (g == "ROD_SCORE" && met == "detect_frac") next
  d <- rows[rows$gene == g, ]
  y <- d[[met]]; is_ms <- d$condition == "MS"
  pp <- perm_p(y, is_ms)
  tests[[paste(g, met)]] <- data.frame(gene = g, class = d$class[1], metric = met, mean_control = mean(y[!is_ms]), mean_ms = mean(y[is_ms]), diff_ms_minus_ctrl = as.numeric(pp["obs"]), perm_p = as.numeric(pp["p"]), mode = pp["mode"])
  register_test(paste0("AMB_", g, "_", met), 3, "ambient", g, as.numeric(pp["obs"]), NA, as.numeric(pp["p"]), 10, 12, paste("retina-level", met, "MS - Control; exact permutation"))
}
tt <- do.call(rbind, tests); rownames(tt) <- NULL
write.csv(tt, file.path(out_root, "ambient_by_condition.csv"), row.names = FALSE)
print(tt, digits = 3)

# rod-marker score vs pooled PSAP log-CPM across the 22 retinas
pa <- read.csv(file.path(v2_res, "A_psap_pooled_lcpm_per_retina.csv"), stringsAsFactors = FALSE)
rs <- rodscore$mean_lognorm[match(pa$retina, rodscore$retina)]
stopifnot(!anyNA(rs), nrow(pa) == 22)
ct <- cor.test(rs, pa$psap); cs <- suppressWarnings(cor.test(rs, pa$psap, method = "spearman"))
write.csv(data.frame(n_retinas = 22, pearson_r = ct$estimate, pearson_p = ct$p.value, spearman_rho = cs$estimate, spearman_p = cs$p.value), file.path(out_root, "ambient_rod_vs_psap.csv"), row.names = FALSE)
register_test("AMB_rod_vs_PSAP_pearson", 3, "ambient", "ROD_SCORE vs PSAP", unname(ct$estimate), NA, ct$p.value, 10, 12, "Pearson r across 22 retinas, rod score vs pooled PSAP log-CPM")
register_test("AMB_rod_vs_PSAP_spearman", 3, "ambient", "ROD_SCORE vs PSAP", unname(cs$estimate), NA, cs$p.value, 10, 12, "Spearman rho across 22 retinas")
log_line("ambient: rod score vs pooled PSAP Pearson r ", round(ct$estimate, 3), " p ", signif(ct$p.value, 3), "; Spearman rho ", round(cs$estimate, 3), " p ", signif(cs$p.value, 3))

# text search of the deposited files
files <- list.files(in_dir, full.names = TRUE, recursive = FALSE)
files <- files[!dir.exists(files) & grepl("\\.(json|csv|tsv|txt|md)$", files)]
pat <- c("CellBender", "ambient", "SoupX", "corrected")
out <- c(paste("files searched:", paste(basename(files), collapse = ", ")))
for (p in pat) {
  hit <- unlist(lapply(files, function(f) { l <- readLines(f, warn = FALSE); h <- grep(p, l, ignore.case = TRUE, value = TRUE); if (length(h)) paste0(basename(f), ": ", substr(h, 1, 300)) }))
  out <- c(out, paste0("[", p, "] ", if (length(hit)) paste(hit, collapse = " || ") else "not stated"))
}
writeLines(out, file.path(out_root, "ambient_text_search.txt"))
cat(out, sep = "\n")
