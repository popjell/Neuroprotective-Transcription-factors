# Final figure export - TWO PDFs: one with TFs as units, one with clusters.
#
# Each PDF has 4 pages:
#   p1: G vs U coverage, color = Cliff's delta (n>=3, shrunk)  [2x2 grid, all groups]
#   p2: log2 specificity bars                                 [2x2 grid, all groups]
#   p3: x = Cliff's delta (n>=3, shrunk), y = log2 specificity, static purple
#                                                             [2x2 grid, all groups]
#   p4: size = specificity, color = Cliff's delta, all groups [own tall page]
#
# Color scale: soft RdBu diverging, grey pinned at 0, stops at the quartiles of
# the negative/positive tails (even colour usage, top positives highlighted).
# p4: size capped to avoid circle overlap, explicit specificity legend breaks.
#
# Run from the repo ROOT (script 15 uses relative paths).

suppressMessages({
  library(dplyr); library(ggplot2); library(gridExtra); library(RColorBrewer)
})

groups <- c("GO_cell_death", "GO_metabolism", "LLM_cell_death", "LLM_metabolism")

# ---- inputs needed for expressible-n ----
expr <- read.table("Expression_data.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
tgt <- expr %>% dplyr::select(symbol, EAE_Combined_Score)
tgt$EAE_Combined_Score <- as.numeric(tgt$EAE_Combined_Score)

expressible_n <- function(target_str, G_genes) {
  tg <- strsplit(target_str, ",\\s*")[[1]]
  tg <- tg[tg != "" & !is.na(tg)]
  sum(tgt$symbol %in% intersect(tg, G_genes) & !is.na(tgt$EAE_Combined_Score))
}

# ---- soft diverging scale (grey at 0, even colour usage) ----
make_cliff_scale <- function(v) {
  v <- v[!is.na(v)]
  L <- max(abs(v))
  absv <- abs(v[v != 0])
  qs <- if (length(absv) > 0) as.numeric(quantile(absv, c(0.25, 0.5, 0.75))) else numeric(0)
  stops <- unique(sort(c(-L, -rev(qs), 0, qs, L)))
  positions <- (stops + L) / (2 * L)
  cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(stops))
  scale_color_gradientn(
    colors = cols, values = positions, limits = c(-L, L),
    oob = scales::squish, na.value = "grey85",
    name = "Cliff's delta\n(n>=3, shrunk)",
    breaks = signif(pretty(c(-L, L)), 2)
  )
}

# ---- build cl + dominant_tf_per_cluster (Jaccard, h = 1.0, as in pipeline) ----
build_clusters <- function() {
  source("utils/helpers.R")
  net <- load_tf_network("reformatted_edges.csv")
  all_tfs <- net$all_tfs; tf_targets <- net$tf_targets
  all_genes <- unique(unlist(tf_targets))
  gene_idx <- setNames(seq_along(all_genes), all_genes)
  bin_mat <- matrix(0L, length(all_tfs), length(all_genes),
                    dimnames = list(all_tfs, all_genes))
  for (tf in all_tfs) bin_mat[tf, gene_idx[tf_targets[[tf]]]] <- 1L
  inter_mat <- bin_mat %*% t(bin_mat)
  union_mat <- outer(rowSums(bin_mat), rowSums(bin_mat), "+") - inter_mat
  jac_mat <- inter_mat / union_mat; diag(jac_mat) <- 1
  hc <- hclust(as.dist(1 - jac_mat), method = "ward.D2")
  cl <- cutree(hc, h = 1.0)
  cluster_ids <- sort(unique(cl))
  dominant_tf_per_cluster <- data.frame(
    cluster = cluster_ids, dominant_tf = NA_character_,
    label = NA_character_, color = rainbow(length(cluster_ids)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(dominant_tf_per_cluster))) {
    cs <- dominant_tf_per_cluster$cluster[i]
    members <- names(cl[cl == cs])
    tf_target_counts <- sapply(members, function(tf) length(tf_targets[[tf]]))
    dom_tf <- names(which.max(tf_target_counts))
    dominant_tf_per_cluster$dominant_tf[i] <- dom_tf
    dominant_tf_per_cluster$label[i] <- paste0(dom_tf, "_cluster")
  }
  list(cl = cl, dominant_tf_per_cluster = dominant_tf_per_cluster)
}

# ---- render one full 4-page PDF for a unit type ----
render_pdf <- function(use_clusters, outfile) {
  assign("readline",
         function(prompt = "") if (use_clusters) "y" else "n",
         envir = .GlobalEnv)
  if (use_clusters) {
    built <- build_clusters()
    assign("cl", built$cl, envir = .GlobalEnv)
    assign("dominant_tf_per_cluster", built$dominant_tf_per_cluster, envir = .GlobalEnv)
  }
  pdf(NULL)
  invisible(capture.output(source("15.coverage_comparasion.R", local = FALSE)))
  dev.off()

  pdata <- list(); comb_rows <- list()
  for (nm in groups) {
    pd <- get(paste0(nm, "_p3"), envir = .GlobalEnv)$data
    G_genes <- read.table(paste0("filtered_genelists/", nm, "_genes.txt"), header = FALSE)$V1
    pd$n_expr <- sapply(pd$targets_in_G, expressible_n, G_genes = G_genes)
    pd$cliff_clean <- ifelse(pd$n_expr < 3, NA_real_,
                             pd$cliff_delta * pd$n_expr / (pd$n_expr + 5))
    pdata[[nm]] <- pd
    comb_rows[[length(comb_rows) + 1]] <- data.frame(
      functional_group = nm, TF = pd$TF,
      specificity = pd$specificity, n_expr = pd$n_expr,
      cliff_clean = pd$cliff_clean, stringsAsFactors = FALSE
    )
  }
  comb <- do.call(rbind, comb_rows)
  comb$functional_group <- factor(comb$functional_group, levels = groups)
  cliff_scale <- make_cliff_scale(comb$cliff_clean)

  cov_xlim   <- range(unlist(lapply(pdata, function(p) p$G_coverage_pct)))
  cov_ylim   <- range(unlist(lapply(pdata, function(p) p$U_coverage_pct)))
  spec_ylim  <- range(unlist(lapply(pdata, function(p) p$log2spec)), na.rm = TRUE)
  cliff_xlim <- range(unlist(lapply(pdata, function(p) p$cliff_clean)), na.rm = TRUE)

  size_range  <- if (use_clusters) c(4, 12) else c(0.5, 6.5)
  size_breaks <- if (use_clusters) c(0.8, 1.1, 1.4, 1.75) else c(0.9, 1.1, 1.3, 1.5, 1.7)

  p1_list <- list(); p2_list <- list(); p3_list <- list()
  for (nm in groups) {
    pd <- pdata[[nm]]

    p1_list[[nm]] <- ggplot(pd, aes(x = G_coverage_pct, y = U_coverage_pct, color = cliff_clean)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7) +
      geom_point(size = 3, alpha = 0.85) +
      geom_text(aes(label = TF), vjust = -0.8, size = 2.6, show.legend = FALSE) +
      cliff_scale +
      labs(title = paste0(nm, ": ", unit_label, " Target Coverage"),
           x = paste0("% of Functional Group G Covered by ", unit_label),
           y = paste0("% of U Covered by ", unit_label)) +
      coord_cartesian(xlim = cov_xlim, ylim = cov_ylim) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(), axis.text = element_text(color = "black"))
    print(p1_list[[nm]])

    p2_list[[nm]] <- pd %>%
      dplyr::filter(is.finite(log2spec)) %>%
      mutate(TF = reorder(TF, log2spec)) %>%
      ggplot(aes(x = TF, y = log2spec)) +
      geom_col() +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
      coord_flip(ylim = spec_ylim) +
      labs(title = paste0(nm, ": ", unit_label, " Functional Specificity"),
           x = unit_label, y = expression(log[2]~"Specificity")) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    print(p2_list[[nm]])

    p3_list[[nm]] <- ggplot(pd, aes(x = cliff_clean, y = log2spec)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
      geom_point(size = 3.5, alpha = 0.9, color = "purple") +
      geom_text(aes(label = TF), vjust = -0.8, size = 2.6, color = "purple") +
      labs(title = paste0(nm, ": ", unit_label, " Specificity vs Effect"),
           x = "Cliff's delta (n>=3, shrunk)",
           y = expression(log[2]~"Specificity")) +
      coord_cartesian(xlim = cliff_xlim, ylim = spec_ylim) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())
    print(p3_list[[nm]])
  }

  p4 <- ggplot(comb, aes(x = functional_group, y = TF,
                         size = specificity^1.5, color = cliff_clean)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = size_range,
                          breaks = size_breaks,
                          name = "Specificity") +
    cliff_scale +
    labs(title = paste0(unit_label, " Functional Specificity Across Disease Programs"),
         x = "Functional Group", y = unit_label) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  print(p4)

  # pages 1-3: landscape 2x2 grids; p4: own tall portrait page; then merge
  grid_file <- tempfile(fileext = ".pdf")
  p4_file   <- tempfile(fileext = ".pdf")
  pdf(grid_file, width = 11, height = 8.5, onefile = TRUE)
  for (pl in list(p1_list, p2_list, p3_list)) {
    gridExtra::grid.arrange(grobs = pl, ncol = 2, nrow = 2)
  }
  dev.off()
  pdf(p4_file, width = 8.5, height = 11, onefile = TRUE)
  print(p4)
  dev.off()
  qpdf::pdf_combine(c(grid_file, p4_file), outfile)

  write.csv(comb, sub("\\.pdf$", "_data.csv", outfile), row.names = FALSE)
  cat("Wrote", outfile, "(pages 1-3: p1/p2/p3 2x2, page 4: p4 portrait)\n")
}

render_pdf(FALSE, "temp_eda/tf_plots.pdf")
render_pdf(TRUE,  "temp_eda/cluster_plots.pdf")
