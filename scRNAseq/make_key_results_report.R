# Make a short, separate report from the saved Schirmer and Trobisch results.
# Run from the scRNAseq folder: Rscript make_key_results_report.R "Rank 1"
library(ggplot2)

target_rank <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "Rank 2"
stopifnot(target_rank %in% c("Rank 1", "Rank 2"))
slug <- if (target_rank == "Rank 1") "rank1" else "rank2"
report_file <- if (target_rank == "Rank 1") "scRNAseq_rank1_results.html" else "scRNAseq_key_results.html"
tro_plot <- if (target_rank == "Rank 1") "key_l2_l3_rank1_volcano.png" else "key_l2_l3_volcano.png"
sch_plot <- if (target_rank == "Rank 1") "key_en_l2_3_rank1_volcano.png" else "key_en_l2_3_volcano.png"

tro_sets <- read.csv("Trobisch_2022/results/full_atlas_gene_set_tests.csv")
sch_sets <- read.csv("Schirmer_2019/results/rank_set_method_comparison.csv")
tro_scores <- read.csv("Trobisch_2022/results/single_sample_rank_scores.csv")
tro_tests <- read.csv("Trobisch_2022/results/single_sample_rank_tests.csv")
tro_genes <- read.csv("Trobisch_2022/results/cortex_integrated_rank_genes.csv")

tro <- subset(tro_sets, subtype == "L2-L3-EN" & rank == target_rank & method == "ROAST")
sch <- subset(sch_sets, subtype == "EN-L2-3" & rank == target_rank & method == "roast")
sing <- subset(tro_tests, subtype == "L2-L3-EN" & rank == target_rank & method == "singscore")
donors <- subset(tro_scores, subtype == "L2-L3-EN" & rank == target_rank & method == "singscore")
donors$condition <- factor(donors$condition, levels = c("Control", "MS"))
stopifnot(nrow(tro) == 1, nrow(sch) == 1, nrow(sing) == 1, nrow(donors) == 8)
tro_hit <- sum(tro_sets$rank == target_rank & tro_sets$method == "ROAST" & tro_sets$q_method < 0.05, na.rm = TRUE)
tro_tested <- sum(tro_sets$rank == target_rank & tro_sets$method == "ROAST" & is.finite(tro_sets$p_value))
sch_hit <- sum(sch_sets$rank == target_rank & sch_sets$method == "roast" & sch_sets$q_within_method < 0.05,
               na.rm = TRUE)
sch_tested <- sum(sch_sets$rank == target_rank & sch_sets$method == "roast" & is.finite(sch_sets$p_value))

# Show every cortical donor, not just a group average.
p <- ggplot(donors, aes(condition, score, color = condition)) +
  geom_point(size = 3.5, position = position_jitter(width = 0.07, height = 0)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, color = "#333333") +
  scale_color_manual(values = c(Control = "#4274A5", MS = "#C46B42")) +
  labs(title = paste("Trobisch L2-L3-EN:", target_rank, "score by donor"),
       subtitle = "Four control and four MS cortical donors; line is the group median",
       x = NULL, y = "singscore") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")
ggsave(paste0("Trobisch_2022/plots/key_l2_l3_", slug, "_donors.png"), p,
       width = 7, height = 4.6, dpi = 180)

# Compare gene-level MS-control results in the closest cortical labels.
tro_de <- subset(read.csv("Trobisch_2022/results/cortex_integrated_de.csv"),
                 subtype == "L2-L3-EN")
tro_de$effect <- tro_de$logFC
tro_de$p <- tro_de$PValue
tro_de$q <- tro_de$q_global
tro_de$symbol <- tro_de$gene
sch_de <- read.csv("Schirmer_2019/results/deseq2_results.csv")
sch_de$q_global <- p.adjust(sch_de$pvalue, method = "BH")
sch_de <- subset(sch_de, subtype == "EN-L2-3")
sch_de$effect <- sch_de$log2FoldChange
sch_de$p <- sch_de$pvalue
sch_de$q <- sch_de$q_global
sch_de$symbol <- sch_de$gene

volcano <- function(d, title, path) {
  d <- d[is.finite(d$effect) & is.finite(d$p), ]
  d$height <- -log10(pmax(d$p, 1e-300))
  d$hit <- is.finite(d$q) & d$q < 0.05
  labels <- if (target_rank == "Rank 1") "PSAP" else c("PSAP", "ITM2C")
  highlight <- d[d$symbol %in% labels, ]
  p <- ggplot(d, aes(effect, height)) +
    geom_point(aes(color = hit), size = 0.8, alpha = 0.65) +
    geom_point(data = highlight, size = 2.7, color = "#174E78") +
    geom_text(data = highlight, aes(label = symbol), vjust = -0.7,
              color = "#174E78", check_overlap = TRUE, size = 3.7) +
    scale_color_manual(values = c(`FALSE` = "#BFC4C8", `TRUE` = "#C46B42")) +
    coord_cartesian(xlim = c(-3, 3), ylim = c(0, 10)) +
    labs(title = title, subtitle = "Orange: q < 0.05 across subtype × gene tests",
         x = "MS − control log2 fold change", y = "−log10(p value)") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  ggsave(path, p, width = 7, height = 4.8, dpi = 180)
}
volcano(tro_de, "Trobisch cortical L2-L3-EN",
        paste0("Trobisch_2022/plots/", tro_plot))
volcano(sch_de, "Schirmer EN-L2-3",
        paste0("Schirmer_2019/plots/", sch_plot))

# The exact donor-label check shows how much four-versus-four tests depend on assumptions.
score_rows <- subset(tro_scores, method == "singscore" & region == "Cortex")
permutation <- lapply(split(score_rows,
                            interaction(score_rows$subtype, score_rows$rank, drop = TRUE)),
  function(d) {
    n_ms <- sum(d$condition == "MS")
    if (n_ms < 2 || nrow(d) - n_ms < 2) return(NULL)
    observed <- mean(d$score[d$condition == "MS"]) -
      mean(d$score[d$condition == "Control"])
    groups <- combn(nrow(d), n_ms)
    possible <- apply(groups, 2, function(i)
      mean(d$score[i]) - mean(d$score[-i]))
    data.frame(subtype = d$subtype[1], rank = d$rank[1],
               p = mean(abs(possible) >= abs(observed) - 1e-12))
  })
permutation <- do.call(rbind, permutation)
permutation$q <- p.adjust(permutation$p, method = "BH")
exact <- subset(permutation, subtype == "L2-L3-EN" & rank == target_rank)

fmt <- function(x) formatC(x, format = "g", digits = 3)
lead_genes <- subset(tro_genes, subtype == "L2-L3-EN" & rank == target_rank &
                       q_global < 0.05 & logFC > 0)
lead_genes <- paste(lead_genes$human_symbol, collapse = " and ")

# Show donor scores and tests for every neuronal label.
{
  sch_scores <- read.csv("Schirmer_2019/results/single_sample_rank_scores.csv")
  sch_tests <- read.csv("Schirmer_2019/results/single_sample_rank_tests.csv")
  all_tro <- subset(tro_tests, rank == target_rank & method == "singscore")
  all_sch <- subset(sch_tests, rank == target_rank & method == "singscore")
  tro_points <- subset(tro_scores, rank == target_rank & method == "singscore")
  sch_points <- subset(sch_scores, rank == target_rank & method == "singscore")
  en_labels <- all_tro$subtype[grepl("EN$", all_tro$subtype)]
  other_labels <- setdiff(all_tro$subtype, en_labels)

  score_plot <- function(d, labels, title, path, show_region = FALSE) {
    d <- d[d$subtype %in% labels, ]
    d$subtype <- factor(d$subtype, levels = labels)
    d$condition <- factor(d$condition, levels = c("Control", "MS"))
    absent <- setdiff(labels, as.character(unique(d$subtype)))
    missing <- data.frame(subtype = factor(absent, levels = labels),
                          condition = rep("Control", length(absent)),
                          score = rep(0, length(absent)),
                          label = rep("Not tested", length(absent)))
    p <- ggplot(d, aes(condition, score, color = condition)) +
      geom_point(aes(shape = region), size = 1.9, alpha = 0.85,
                 position = position_jitter(width = 0.12, height = 0)) +
      facet_wrap(~subtype, ncol = 3, scales = "free_y", drop = FALSE) +
      scale_color_manual(values = c(Control = "#4274A5", MS = "#C46B42")) +
      labs(title = title, x = NULL, y = "singscore", shape = "Region") +
      guides(color = "none") +
      theme_minimal(base_size = 11) +
      theme(legend.position = if (show_region) "bottom" else "none",
            strip.text = element_text(face = "bold"))
    if (nrow(missing)) p <- p + geom_text(
      data = missing, aes(condition, score, label = label),
      inherit.aes = FALSE, color = "#56636e", size = 3.5)
    ggsave(path, p, width = 10, height = if (length(labels) > 9) 10 else 8,
           dpi = 180)
  }
  score_plot(tro_points, en_labels, paste("Trobisch:", target_rank, "excitatory subtypes"),
             paste0("Trobisch_2022/plots/key_", slug, "_singscore_en.png"), TRUE)
  score_plot(tro_points, other_labels, paste("Trobisch:", target_rank, "other neuronal subtypes"),
             paste0("Trobisch_2022/plots/key_", slug, "_singscore_other.png"), TRUE)
  score_plot(sch_points, all_sch$subtype, paste("Schirmer:", target_rank, "neuronal subtypes"),
             paste0("Schirmer_2019/plots/key_", slug, "_singscore_all.png"))

  show_value <- function(x) ifelse(is.finite(x), fmt(x), "—")
  score_table <- function(d) paste(vapply(seq_len(nrow(d)), function(i) {
    x <- d[i, ]
    paste0("<tr><td>", x$subtype, "</td><td>", x$control_donors, " / ",
           x$ms_donors, "</td><td>", show_value(x$genes_tested),
           "</td><td>", show_value(x$ms_minus_control), "</td><td>",
           show_value(x$p_value), "</td><td>", show_value(x$q_method),
           "</td><td>", x$status, "</td></tr>")
  }, character(1)), collapse = "")
  table_head <- "<thead><tr><th>Subtype</th><th>Donors, control / MS</th><th>Genes</th><th>MS − control</th><th>p</th><th>q</th><th>Status</th></tr></thead>"
  all_scores_html <- c(
    "<h2>Singscore across all neuronal subtypes</h2>",
    "<p>Each point is one donor-subtype profile. Trobisch includes cortex and cerebellum; tests adjust for region. Panel y-axes vary by subtype. Empty panels have no testable score.</p>",
    "<h3>Trobisch</h3>",
    paste0("<figure><img src='Trobisch_2022/plots/key_", slug, "_singscore_en.png' alt='Trobisch ", target_rank, " singscores across excitatory subtypes'><figcaption>Trobisch excitatory-neuron donor scores. Point shape denotes region.</figcaption></figure>"),
    paste0("<figure><img src='Trobisch_2022/plots/key_", slug, "_singscore_other.png' alt='Trobisch ", target_rank, " singscores across other neuronal subtypes'><figcaption>Trobisch other neuronal donor scores. Point shape denotes region.</figcaption></figure>"),
    paste0("<table class='singscore-table'><caption>Trobisch ", target_rank, " singscore tests</caption>"), table_head,
    paste0("<tbody>", score_table(all_tro), "</tbody></table>"),
    "<h3>Schirmer</h3>",
    paste0("<figure><img src='Schirmer_2019/plots/key_", slug, "_singscore_all.png' alt='Schirmer ", target_rank, " singscores across all neuronal subtypes'><figcaption>Schirmer neuronal donor scores.</figcaption></figure>"),
    paste0("<table class='singscore-table'><caption>Schirmer ", target_rank, " singscore tests</caption>"), table_head,
    paste0("<tbody>", score_table(all_sch), "</tbody></table>"),
    "<p>q values are Benjamini-Hochberg corrected within singscore across both Rank 1 and Rank 2 subtype tests in each dataset.</p>")
}

html <- c(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  paste0("<title>Neuronal ", target_rank, " gene-set results</title>"),
  "<style>body{font:15px/1.5 Arial,sans-serif;max-width:1040px;margin:36px auto;padding:0 24px;color:#242b30}h1{font-size:27px;font-weight:600;margin:0 0 4px}h2{font-size:19px;font-weight:600;margin:34px 0 10px;border-bottom:1px solid #b9c2c9;padding-bottom:5px}h3{font-size:16px;margin:24px 0 8px}.subtitle{color:#56636e;margin:0 0 22px}figure{margin:16px 0 26px}figure img{width:100%;height:auto;border:1px solid #d7dde1}figcaption{font-size:13px;color:#4c5962;margin-top:6px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}.grid figure{margin:0 0 18px}table{border-collapse:collapse;width:100%;margin:12px 0 20px}th,td{border-bottom:1px solid #d7dde1;padding:8px;text-align:left}th{font-weight:600;background:#f4f6f7}caption{text-align:left;font-weight:600;margin-bottom:6px}.singscore-table{font-size:13px}.methods{font-size:13px;color:#4c5962;border-top:1px solid #b9c2c9;padding-top:12px}a{color:#14598b}@media(max-width:760px){.grid{grid-template-columns:1fr}table{display:block;overflow-x:auto}}</style></head><body>",
  paste0("<h1>Neuronal ", target_rank, " gene-set results in human MS</h1>"),
  "<p class='subtitle'>Schirmer et al. 2019 and Trobisch et al. 2022 · MS versus control · donor-level analysis</p>",
  "<h2>Gene-set tests</h2>",
  "<table><thead><tr><th>Dataset and neuronal label</th><th>Direction</th><th>Genes tested</th><th>Donors, control / MS</th><th>ROAST p</th><th>ROAST q</th></tr></thead><tbody>",
  paste0("<tr><td>Trobisch · L2-L3-EN · cortex</td><td>MS-up</td><td>", tro$genes_tested,
         " / 31</td><td>", tro$control_donors, " / ", tro$ms_donors,
         "</td><td>", fmt(tro$p_value), "</td><td>", fmt(tro$q_method), "</td></tr>"),
  paste0("<tr><td>Schirmer · EN-L2-3</td><td>MS-up</td><td>", sch$genes_tested,
         " / 31</td><td>9 / 9</td><td>", fmt(sch$p_value),
         "</td><td>", fmt(sch$q_within_method), "</td></tr>"),
  "</tbody></table>",
  paste0("<p>", target_rank, " ROAST tests with q &lt; 0.05 across neuronal subtypes: Trobisch ",
         tro_hit, " / ", tro_tested, " (all MS-up); Schirmer ", sch_hit,
         " / ", sch_tested, ".</p>"),
  paste0("<h2>Donor-level ", target_rank, " scores</h2>"),
  paste0("<figure><img src='Trobisch_2022/plots/key_l2_l3_", slug, "_donors.png' alt='Eight cortical donor scores for ", target_rank, "'><figcaption>Trobisch L2-L3-EN; one point per donor. Singscore: p ", fmt(sing$p_value), ", q ", fmt(sing$q_method), ". Exact donor-label permutation: p ", fmt(exact$p), ", q ", fmt(exact$q), " across cortical subtype/list tests.</figcaption></figure>"),
  "<h2>Cell annotations</h2>",
  "<div class='grid'>",
  "<figure><img src='Trobisch_2022/plots/full_atlas_umap.png' alt='Trobisch full integrated atlas UMAP'><figcaption>Trobisch integrated atlas UMAP with author annotations.</figcaption></figure>",
  "<figure><img src='Schirmer_2019/plots/umap.png' alt='Schirmer neuronal UMAP'><figcaption>Schirmer neuronal UMAP with subtype annotations.</figcaption></figure>",
  "</div>",
  "<h2>Gene-level differential expression</h2>",
  "<div class='grid'>",
  paste0("<figure><img src='Trobisch_2022/plots/", tro_plot, "' alt='Trobisch L2-L3-EN volcano plot'><figcaption>Trobisch L2-L3-EN. ", lead_genes, if (target_rank == "Rank 1") " passes" else " pass", " global gene-level FDR (q &lt; 0.05).</figcaption></figure>"),
  paste0("<figure><img src='Schirmer_2019/plots/", sch_plot, "' alt='Schirmer EN-L2-3 volcano plot'><figcaption>Schirmer EN-L2-3. ", if (target_rank == "Rank 1") "PSAP is" else "PSAP and ITM2C are", " highlighted for comparison.</figcaption></figure>"),
  "</div>",
  all_scores_html,
  "<h2>Methods</h2>",
  "<p class='methods'>Raw counts were summed by donor and neuronal subtype. ROAST tested MS versus control within each subtype; q values use Benjamini-Hochberg correction within the method across testable subtype/list comparisons. Singscore used ranked, normalized donor expression; donor-level tests compare MS versus control and adjust for region in Trobisch when both regions are represented. Singscore q values are corrected within the method across testable subtype/list comparisons. Volcano plots use gene-level q values across subtype-by-gene tests and show −3 to 3 log2 fold change and 0 to 10 on −log10(p). Trobisch cortical samples overlap Schirmer samples; subtype labels and expressed-gene sets differ. <a href='scRNAseq_report.html'>Detailed report</a> · <a href='Trobisch_2022/results/full_atlas_gene_set_tests.csv'>Complete Trobisch gene-set results</a>.</p>",
  "</body></html>")
writeLines(html, report_file)
