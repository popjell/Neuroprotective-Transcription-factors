# Analysis E (Tier 3): PSAP and the 31-gene set against expression-matched random genes, lysosomal comparators, and the transcriptome-wide shift, in Analyses B and C.
source("common_v2.R")

n_draw <- 10000  # random draws
n_boot <- 2000   # bootstrap resamples across genes
pb <- readRDS(file.path(out_root, "pseudobulk.rds"))
g31 <- read_gene_list(pb$genes)
srcs <- list(B = read.csv(file.path(out_root, "B_all_genes.csv"), stringsAsFactors = FALSE),
             C = read.csv(file.path(out_root, "C_all_genes.csv"), stringsAsFactors = FALSE))

null_rows <- list(); lys_rows <- list(); shift_rows <- list()
for (a in names(srcs)) {
  tab <- srcs[[a]]
  tab$decile <- cut(tab$ave_expr, quantile(tab$ave_expr, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  pool <- split(tab, tab$decile)
  nc <- 10; nm <- 12
  # PSAP against random single genes from its own decile (excluding PSAP)
  ps <- tab[tab$gene == "PSAP", ]
  cand <- pool[[ps$decile]]; cand <- cand[cand$gene != "PSAP", ]
  for (stat in c("t", "log2FC")) {
    draw <- sample(cand[[stat]], n_draw, replace = TRUE)
    pct <- mean(draw <= ps[[stat]])
    pup <- (1 + sum(draw >= ps[[stat]])) / (n_draw + 1)
    null_rows[[paste(a, "PSAP", stat)]] <- data.frame(analysis = a, target = "PSAP", statistic = stat, observed = ps[[stat]], decile = ps$decile, null_median = median(draw), percentile = pct, empirical_p_upper = pup)
    register_test(paste0("E_", a, "_PSAP_", stat), 3, "E", "PSAP", ps[[stat]], NA, pup, nc, nm, paste0("vs ", n_draw, " expression-matched random genes (decile ", ps$decile, "); p = upper-tail empirical; percentile ", round(100 * pct, 1)))
  }
  # 31-gene set: median log2FC against sets where each gene is replaced by a random gene of the same decile
  gs <- tab[tab$gene %in% g31, ]
  obs <- median(gs$log2FC)
  nul <- replicate(n_draw, median(sapply(gs$decile, function(d) { cc <- pool[[d]]; cc <- cc[!cc$gene %in% g31, ]; sample(cc$log2FC, 1) })))
  pup <- (1 + sum(nul >= obs)) / (n_draw + 1)
  null_rows[[paste(a, "set31")]] <- data.frame(analysis = a, target = paste0("set31 (", nrow(gs), " genes)"), statistic = "median_log2FC", observed = obs, decile = NA, null_median = median(nul), percentile = mean(nul <= obs), empirical_p_upper = pup)
  register_test(paste0("E_", a, "_set31_median_log2FC"), 3, "E", "set31", obs, NA, pup, nc, nm, paste0(nrow(gs), " genes vs ", n_draw, " expression-matched random sets; upper-tail empirical p"))
  # lysosomal comparators
  ly <- tab[tab$gene %in% c(lysosomal, "PSAP"), c("gene", "log2FC", "t", "p")]  # PSAP ranked among the comparators plus itself
  ly$rank_by_log2FC <- rank(-ly$log2FC); ly$rank_by_t <- rank(-ly$t); ly$n_ranked <- nrow(ly)
  ly$analysis <- a
  lys_rows[[a]] <- ly[order(ly$rank_by_log2FC), ]
  log_line("E ", a, ": PSAP t percentile among matched genes ", round(100 * null_rows[[paste(a, "PSAP t")]]$percentile, 1), "; log2FC percentile ", round(100 * null_rows[[paste(a, "PSAP log2FC")]]$percentile, 1),
           "; set31 median log2FC ", round(obs, 4), " null median ", round(median(nul), 4), " p ", signif(pup, 3), "; lysosomal comparators present ", nrow(ly) - 1, " (missing: ", paste(setdiff(lysosomal, tab$gene), collapse = ","), "); PSAP rank by log2FC ",
           ly$rank_by_log2FC[ly$gene == "PSAP"], " of ", nrow(ly), ", by t ", ly$rank_by_t[ly$gene == "PSAP"], " of ", nrow(ly))
  # transcriptome-wide shift with a bootstrap over genes
  bs_med <- replicate(n_boot, median(sample(tab$log2FC, replace = TRUE))); bs_mean <- replicate(n_boot, mean(sample(tab$log2FC, replace = TRUE)))
  shift_rows[[a]] <- data.frame(analysis = a, n_genes = nrow(tab), median_log2FC = median(tab$log2FC), median_ci_low = quantile(bs_med, 0.025), median_ci_high = quantile(bs_med, 0.975),
                                mean_log2FC = mean(tab$log2FC), mean_ci_low = quantile(bs_mean, 0.025), mean_ci_high = quantile(bs_mean, 0.975), frac_positive = mean(tab$log2FC > 0))
}
write.csv(do.call(rbind, null_rows), file.path(out_root, "E_matched_null.csv"), row.names = FALSE)
write.csv(do.call(rbind, lys_rows), file.path(out_root, "E_lysosomal_comparators.csv"), row.names = FALSE)
write.csv(do.call(rbind, shift_rows), file.path(out_root, "E_global_shift.csv"), row.names = FALSE)
print(do.call(rbind, shift_rows))
