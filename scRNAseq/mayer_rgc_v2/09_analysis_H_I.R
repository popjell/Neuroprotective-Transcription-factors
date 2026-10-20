# Analysis H: minimum detectable effects and meta-analysis inputs. Analysis I: look for the authors' module A-F gene lists in the deposited files (read-only).
source("common_v2.R")

pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
rd <- function(f) read.csv(file.path(out_root, f), stringsAsFactors = FALSE)
pick <- function(tab, g, se_col = "se") tab[tab$gene == g, ]
rows <- list()
add <- function(analysis, gene, tab, nc, nm, se_col = "se", fc_col = "log2FC") {
  r <- tab[tab$gene == gene, ]
  rows[[length(rows) + 1]] <<- data.frame(analysis = analysis, gene = gene, log2FC = r[[fc_col]], SE = r[[se_col]], n_control = nc, n_ms = nm, unit = "retina")
}
pooled <- rd("pooled_deseq2_rebuilt.csv")
bu <- rd("B_units_used.csv"); nb_c <- length(unique(bu$retina[bu$condition == "Control"])); nb_m <- length(unique(bu$retina[bu$condition == "MS"]))
for (g in c("PSAP", "ITM2C")) {
  add("pooled_DESeq2", g, pooled, 10, 12, se_col = "lfcSE")
  add("B_subtype_adjusted", g, rd("B_all_genes.csv"), nb_c, nb_m)
  add("C_composition_standardized", g, rd("C_all_genes.csv"), 10, 12)
  add("D_resilient_MS_effect", g, rd("D_resilient_MS_effect_all_genes.csv"), 10, 12)   # n retinas from D unit table below
  add("D_susceptible_MS_effect", g, rd("D_susceptible_MS_effect_all_genes.csv"), 10, 12)
  add("D_interaction_susc_minus_resil", g, rd("D_interaction_all_genes.csv"), 10, 12)
}
meta_in <- do.call(rbind, rows)
meta_in$min_detectable_abs_log2FC_80pct_power <- 2.8 * meta_in$SE
write.csv(meta_in[, c("analysis", "gene", "log2FC", "SE", "n_control", "n_ms", "unit")], file.path(out_root, "meta_inputs.csv"), row.names = FALSE)
pw <- meta_in[meta_in$gene == "PSAP", ]
pw$trobisch_psap_range <- "0.80 to 1.59 (not recomputed)"
pw$mde_below_0.80 <- pw$min_detectable_abs_log2FC_80pct_power < 0.80
write.csv(pw, file.path(out_root, "H_power_psap.csv"), row.names = FALSE)
print(pw[, c("analysis", "log2FC", "SE", "min_detectable_abs_log2FC_80pct_power", "mde_below_0.80")], digits = 3)
log_line("H: minimum detectable |log2FC| (2.8 x SE) for PSAP: ", paste(paste0(pw$analysis, " ", round(pw$min_detectable_abs_log2FC_80pct_power, 3)), collapse = "; "))

# Analysis I: search the deposited text files for module lists
files <- list.files(in_dir, full.names = TRUE)
files <- files[!dir.exists(files) & grepl("\\.(json|csv|txt|tsv|md)$", files)]
hits <- unlist(lapply(files, function(f) { l <- readLines(f, warn = FALSE); h <- grep("module", l, ignore.case = TRUE, value = TRUE); if (length(h)) paste(basename(f), ":", substr(h, 1, 200)) }))
txt <- if (length(hits)) c("lines mentioning 'module' in deposited files:", hits, "", "check whether these are gene lists before use") else "module lists not available in the deposited files"
writeLines(c(paste("files searched:", paste(basename(files), collapse = ", ")), txt), file.path(out_root, "I_author_modules.txt"))
log_line("I: ", txt[1])
