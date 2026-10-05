# Rebuild the compact HTML report from the saved result tables.
rank_list <- read.csv("Schirmer_2019/data/rank_genes.csv")
rank_list <- unique(rank_list[rank_list$use_primary_analysis == "Y",
                              c("rank", "human_symbol")])
schirmer <- read.csv("Schirmer_2019/results/rank_gene_global.csv")
new <- read.csv("Trobisch_2022/results/integrated_neuron_rank_genes.csv")
cortex_file <- "Trobisch_2022/results/cortex_integrated_rank_genes.csv"
cortex <- if (file.exists(cortex_file)) read.csv(cortex_file) else NULL
status <- read.csv("Trobisch_2022/results/integrated_neuron_test_status.csv")
cortex_status_file <- "Trobisch_2022/results/cortex_integrated_test_status.csv"
if (file.exists(cortex_status_file)) {
  cs <- read.csv(cortex_status_file)
  status <- status[status$region != "Cortex", ]
  cs$region <- "Cortex"
  status <- rbind(status, cs[, names(status)])
}

number <- function(x) ifelse(is.na(x), "—", formatC(x, format = "g", digits = 3))
best <- function(data, gene, rank_name, effect, p, q) {
  if (is.null(data)) return("<td>—</td>")
  d <- data[data$human_symbol == gene & data$rank == rank_name &
              is.finite(data[[p]]), ]
  if (!nrow(d)) return("<td>—</td>")
  order_by <- ifelse(is.finite(d[[q]]), d[[q]], Inf)
  i <- if (all(is.infinite(order_by))) which.min(d[[p]]) else which.min(order_by)
  d <- d[i, ]
  direction <- if (d[[effect]] > 0) "↑" else "↓"
  label <- paste0(direction, " ", d$subtype, " · p ", number(d[[p]]),
                  " · q ", number(d[[q]]))
  class <- if (is.finite(d[[q]]) && d[[q]] < 0.05) " class='hit'" else ""
  paste0("<td", class, ">", label, "</td>")
}

rank_table <- function(rank_name) {
  genes <- rank_list$human_symbol[rank_list$rank == rank_name]
  rows <- vapply(genes, function(g) paste0(
    "<tr><td>", g, "</td>",
    best(schirmer, g, rank_name, "log2FoldChange", "pvalue", "q_global"),
    best(new, g, rank_name, "logFC", "PValue", "q_global"),
    best(cortex, g, rank_name, "logFC", "PValue", "q_global"), "</tr>"), "")
  c(paste0("<h3>", rank_name, "</h3>"), "<div class='table-wrap'><table>",
    "<thead><tr><th>Gene</th><th>Schirmer</th><th>Trobisch cerebellum/spinal</th><th>Trobisch cortex*</th></tr></thead>",
    "<tbody>", rows, "</tbody></table></div>")
}

labels <- sort(unique(status$subtype))
coverage_rows <- vapply(labels, function(s) {
  d <- status[status$subtype == s, ]
  cells <- function(region) {
    z <- d[d$region == region, ]
    if (!nrow(z)) return("—")
    paste0(z$cells, " · ", z$status)
  }
  paste0("<tr><td>", s, "</td><td>", cells("Cortex"), "</td><td>",
         cells("Cerebellum"), "</td><td>", cells("Spinal_cord"), "</td></tr>")
}, "")

sch_broad <- read.csv("Schirmer_2019/results/broad_neuron_rank_tests.csv")
new_broad <- read.csv("Trobisch_2022/results/integrated_neuron_broad_ucell.csv")
cortex_broad_file <- "Trobisch_2022/results/cortex_integrated_broad_ucell.csv"
cortex_broad <- if (file.exists(cortex_broad_file)) read.csv(cortex_broad_file) else NULL
cortex_sensitivity_file <- "Trobisch_2022/results/cortex_broad_ucell_sensitivity.csv"
cortex_sensitivity <- if (file.exists(cortex_sensitivity_file))
  read.csv(cortex_sensitivity_file) else NULL
cortex_loo_file <- "Trobisch_2022/results/cortex_rank_set_leave_one_out.csv"
cortex_loo <- if (file.exists(cortex_loo_file)) read.csv(cortex_loo_file) else NULL
driver_file <- "Trobisch_2022/results/cortex_gene_driver_tests.csv"
drivers <- if (file.exists(driver_file)) read.csv(driver_file) else NULL
camera <- read.csv("Trobisch_2022/results/integrated_neuron_camera.csv")
full_sets <- read.csv("Trobisch_2022/results/full_atlas_gene_set_tests.csv")
regional_sets <- read.csv("Trobisch_2022/results/full_atlas_region_gene_set_tests.csv")
sch_scores <- read.csv("Schirmer_2019/results/single_sample_rank_tests.csv")
tro_scores <- read.csv("Trobisch_2022/results/single_sample_rank_tests.csv")
tro_score_regions <- read.csv("Trobisch_2022/results/single_sample_rank_regional_tests.csv")
score_methods <- c("ssGSEA", "GSVA", "singscore", "AUCell", "PLAGE", "zscore", "minmax")
score_count <- function(d, method, direction) sum(
  d$method == method & is.finite(d$q_method) & d$q_method < 0.05 &
    if (direction == "up") d$ms_minus_control > 0 else d$ms_minus_control < 0)
score_rows <- vapply(score_methods, function(method) paste0(
  "<tr><td>", method, "</td><td>",
  sum(sch_scores$method == method & is.finite(sch_scores$p_value)),
  "</td><td>", score_count(sch_scores, method, "up"), " / ",
  score_count(sch_scores, method, "down"), "</td><td>",
  sum(tro_scores$method == method & is.finite(tro_scores$p_value)),
  "</td><td>", score_count(tro_scores, method, "up"), " / ",
  score_count(tro_scores, method, "down"), "</td></tr>"), "")
sch_tests <- read.csv("Schirmer_2019/results/donor_rank_tests.csv")
vip <- sch_tests[sch_tests$subtype == "IN-VIP" & sch_tests$rank == "Rank 1", ]
sch_hit <- unique(schirmer$human_symbol[is.finite(schirmer$q_global) &
                                           schirmer$q_global < 0.05])
new_hit <- unique(new$human_symbol[is.finite(new$q_global) & new$q_global < 0.05])
cx_hit <- if (is.null(cortex)) character() else unique(cortex$human_symbol[
  is.finite(cortex$q_global) & cortex$q_global < 0.05])
tested_new <- sum(status$region != "Cortex" & status$status == "Tested")
tested_cx <- sum(status$region == "Cortex" & status$status == "Tested")
full_hit <- function(method) sum(full_sets$method == method &
  is.finite(full_sets$q_method) & full_sets$q_method < 0.05)
full_tested <- sum(full_sets$method == "GSEA" & full_sets$status == "Tested")

set_cell <- function(d) {
  if (is.na(d$q_method)) return(paste0("<td>", d$status, "</td>"))
  class <- if (d$q_method < 0.05) " class='hit'" else ""
  paste0("<td", class, ">", d$direction, " · p ", number(d$p_value),
         " · q ", number(d$q_method), "</td>")
}
set_pairs <- unique(full_sets[, c("subtype", "rank")])
set_rows <- vapply(seq_len(nrow(set_pairs)), function(i) {
  d <- full_sets[full_sets$subtype == set_pairs$subtype[i] &
                   full_sets$rank == set_pairs$rank[i], ]
  g <- d[d$method == "GSEA", ]
  r <- d[d$method == "ROAST", ]
  paste0("<tr><td>", d$subtype[1], "</td><td>", d$rank[1],
         "</td><td>", d$control_donors[1], "/", d$ms_donors[1],
         "</td><td>", if (nzchar(d$regions[1])) d$regions[1] else "—",
         "</td><td>", number(d$genes_tested[1]), "</td>",
         set_cell(g), set_cell(r), "</tr>")
}, "")

html <- c(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>Schirmer and Trobisch neuronal MS analysis</title>",
  "<style>body{font:16px/1.5 Arial,sans-serif;max-width:1100px;margin:30px auto;padding:0 18px;color:#222}h1,h2,h3{line-height:1.2}h2{margin-top:38px}figure{margin:18px 0 28px}img{display:block;width:100%;height:auto;border:1px solid #ddd}figcaption,.note{color:#555;font-size:14px}.note{background:#f5f6f7;padding:12px 16px}.table-wrap{max-height:430px;overflow:auto;border:1px solid #ddd;margin-bottom:24px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:7px 10px;text-align:left;border-bottom:1px solid #e5e5e5}th{position:sticky;top:0;background:#f5f6f7}td:first-child{font-weight:600}tr:nth-child(even){background:#fafafa}.hit{background:#e4f4e7!important;color:#176b32;font-weight:700}a{color:#14598b}</style></head><body>",
  "<h1>Neuronal MS analysis: Schirmer and Trobisch</h1>",
  "<p>Schirmer and Trobisch were compared at donor level within neuronal subtypes. Trobisch is analyzed as one dataset across cortex, cerebellum, and spinal cord; region-by-region results show where a combined signal arises. Its full atlas includes 20 neuronal labels.</p>",
  "<h2>Summary</h2>",
  paste0("<p>In Trobisch's combined, region-adjusted analysis, ", full_hit("GSEA"), " of ", full_tested, " testable subtype/list comparisons pass FDR with preranked GSEA, and ", full_hit("ROAST"), " of ", full_tested, " pass with ROAST. Every significant comparison is MS-up. However, most testable groups are cortex-only; a combined-dataset result does not mean all three regions support it.</p>"),
  paste0("<p>In the region-by-region check, cerebellar GRAN-EN Rank 2 is nominally MS-up by ROAST (p ", number(regional_sets$p_value[regional_sets$scope == "Cerebellum" & regional_sets$subtype == "GRAN-EN" & regional_sets$rank == "Rank 2" & regional_sets$method == "ROAST"]), "), but does not pass correction among cerebellar tests (q ", number(regional_sets$q_region[regional_sets$scope == "Cerebellum" & regional_sets$subtype == "GRAN-EN" & regional_sets$rank == "Rank 2" & regional_sets$method == "ROAST"]), "). Spinal-cord donor-subtype samples are too sparse for these tests. Schirmer broad Rank 1 remains nonsignificant (q ", number(sch_broad$q_two_lists[sch_broad$rank == "Rank 1"]), ").</p>"),
  if (!is.null(cortex_broad)) paste0("<p>Trobisch's finer cortical labels give a separate view of the reused Schirmer donors. The broad UCell scores are MS-up: Rank 1 p ", number(cortex_broad$p_value[cortex_broad$rank == "Rank 1"]), " (q ", number(cortex_broad$q_two_lists[cortex_broad$rank == "Rank 1"]), "); Rank 2 p ", number(cortex_broad$p_value[cortex_broad$rank == "Rank 2"]), " (q ", number(cortex_broad$q_two_lists[cortex_broad$rank == "Rank 2"]), "). This is supportive within the cortical cohort, not independent replication.</p>") else "",
  if (!is.null(cortex_sensitivity)) paste0("<p>The broad cortical scores ", if (all(cortex_sensitivity$q_two_lists < 0.05)) "remain" else "do not all remain", " FDR-significant after separate adjustments for age, RNA integrity, post-mortem interval, median UMIs, detected genes, or mitochondrial percentage. ", if (!is.null(cortex_loo)) paste0(sum(cortex_loo$same_direction == cortex_loo$donors), " of ", nrow(cortex_loo), " FDR-significant cortical subtype/list scores retain their MS-up direction after removing any one donor. ") else "", "With only eight donors, these checks cannot remove unmeasured batch or tissue differences.</p>") else "",
  if (!is.null(drivers)) paste0("<p>The rank lists are not behaving as uniform 31-gene programs. Removing PSAP and ITM2C together removes the broad Rank 2 difference (p ", number(drivers$p_value[drivers$set == "R2_no_PSAP_ITM2C"]), "). Rank 1 remains higher after removing PSAP and PREX1 (p ", number(drivers$p_value[drivers$set == "R1_no_PSAP_PREX1"]), "), but its effect is much smaller.</p>") else "",
  paste0("<p>A competitive gene-set test also finds MS-up enrichment in ", sum(camera$region == "Cortex" & camera$direction == "Up" & camera$q_global < 0.05, na.rm = TRUE), " cortical subtype/list comparisons after correction across all testable Trobisch comparisons. Only a few rank-list genes pass the expression filter in each subtype, so this supports a focused gene module rather than every gene on the lists.</p>"),
  paste0("<p>Schirmer IN-VIP Rank 1 was higher in MS at p ", number(vip$p_value), " (q ", number(vip$q_global), ") in the donor-level mean-expression test. It is a lead, not a confirmed signature. Rank-list genes passing global gene-level FDR were ", if (length(sch_hit)) paste(sch_hit, collapse = ", ") else "none", " in Schirmer, ", if (length(new_hit)) paste(new_hit, collapse = ", ") else "none", " in Trobisch's new regions, and ", if (length(cx_hit)) paste(cx_hit, collapse = ", ") else "none", " in Trobisch cortex. Direction matters: an MS-down gene is not evidence that the rank list is up.</p>"),
  "<h2>Single-sample Rank 1/2 scores</h2>",
  "<p>Seven scores were calculated from the same normalized donor-by-subtype pseudobulk profiles: ssGSEA, GSVA, singscore, AUCell, PLAGE, z-score, and gene-wise 0–1 min-max averaging. Each score was tested for MS versus control within every eligible neuronal label. Trobisch used cortex, cerebellum, and spinal cord together, adjusting the disease effect for region when a label had usable donors in multiple regions. The table counts q &lt; 0.05 within each method; the two numbers are MS-up / MS-down. Methods reuse the same donors and genes, so agreement is a robustness check, not independent replication.</p>",
  "<div class='table-wrap'><table><thead><tr><th>Method</th><th>Schirmer testable</th><th>Schirmer up / down</th><th>Trobisch testable</th><th>Trobisch up / down</th></tr></thead><tbody>",
  score_rows, "</tbody></table></div>",
  "<figure><img src='Schirmer_2019/plots/single_sample_rank_tests.png' alt='Schirmer rank-list score tests'><figcaption>All 9 Schirmer neuronal labels. The single GSVA hit is EN-PYR Rank 2, MS-down, and is not supported by the other methods.</figcaption></figure>",
  "<figure><img src='Trobisch_2022/plots/single_sample_rank_tests.png' alt='Trobisch rank-list score tests'><figcaption>All 20 Trobisch neuronal labels are included; light tiles are untestable, not negative results. The significant combined scores are mostly cortical.</figcaption></figure>",
  paste0("<p>Among Trobisch's non-cortical donors, only cerebellar GRAN-EN and VIP-IN had score comparisons with at least two donors per condition; spinal cord had none. For cerebellar GRAN-EN Rank 2, AUCell (q ", number(tro_score_regions$q_region_method[tro_score_regions$region == "Cerebellum" & tro_score_regions$subtype == "GRAN-EN" & tro_score_regions$rank == "Rank 2" & tro_score_regions$method == "AUCell"]), "), z-score (q ", number(tro_score_regions$q_region_method[tro_score_regions$region == "Cerebellum" & tro_score_regions$subtype == "GRAN-EN" & tro_score_regions$rank == "Rank 2" & tro_score_regions$method == "zscore"]), "), and min-max (q ", number(tro_score_regions$q_region_method[tro_score_regions$region == "Cerebellum" & tro_score_regions$subtype == "GRAN-EN" & tro_score_regions$rank == "Rank 2" & tro_score_regions$method == "minmax"]), ") were MS-up after correction within each region and method, but each group has only three control and three MS donors; ssGSEA, GSVA, singscore, and PLAGE were not significant. This is a small, method-dependent regional lead, not broad validation.</p>"),
  "<p class='note'>Complete <a href='Schirmer_2019/results/single_sample_rank_tests.csv'>Schirmer tests</a>, <a href='Trobisch_2022/results/single_sample_rank_tests.csv'>combined Trobisch tests</a>, and <a href='Trobisch_2022/results/single_sample_rank_regional_tests.csv'>Trobisch region-by-region checks</a> include every label, list, method, donor count, direction, p value, q value, and untestable status. Individual donor scores are saved beside them.</p>",
  "<h2>Schirmer 2019</h2>",
  "<figure><img src='Schirmer_2019/plots/umap.png' alt='Schirmer neuronal subtype UMAP'><figcaption>Neuronal subtypes used for donor-level pseudobulk.</figcaption></figure>",
  "<figure><img src='Schirmer_2019/plots/volcano.png' alt='Schirmer pseudobulk volcano plots'><figcaption>MS versus control gene-level differential expression within subtype.</figcaption></figure>",
  "<figure><img src='Schirmer_2019/plots/donor_rank_scores.png' alt='Schirmer donor-level rank scores'><figcaption>Rank-list mean normalized expression; one point per patient. The score uses genes passing the expression filter, with the count given in the results table.</figcaption></figure>",
  "<h2>Trobisch 2022</h2>",
  "<figure><img src='Trobisch_2022/plots/full_atlas_umap.png' alt='Full Trobisch atlas UMAP with 38 labels'><figcaption>The authors' full atlas: 38 annotations, including all 20 integrated neuronal labels.</figcaption></figure>",
  "<h3>Full-dataset Rank 1 and Rank 2 tests</h3>",
  "<figure><img src='Trobisch_2022/plots/full_atlas_gene_sets.png' alt='GSEA and ROAST results for all 20 Trobisch neuronal labels'><figcaption>Combined Trobisch donors; MS effect adjusted for region. Dark blue means q &lt; 0.05 within the method. Grey means not significant; near-white means untestable.</figcaption></figure>",
  "<p>Raw counts were summed per donor and neuronal label. Each donor-label sample needed at least 10 nuclei. Genes were retained if expressed in an adequately sampled region, then TMM/voom and a region-adjusted MS-control model supplied the ranked statistic for GSEA and the model for ROAST. Both methods were corrected across all testable combined subtype/list comparisons. The table shows every label, including untestable ones; C/MS is the number of usable control/MS donor samples. Cortex-only rows use only the eight donors shared with Schirmer.</p>",
  "<div class='table-wrap'><table><thead><tr><th>Label</th><th>List</th><th>C/MS</th><th>Regions tested</th><th>Genes</th><th>GSEA</th><th>ROAST</th></tr></thead><tbody>",
  set_rows, "</tbody></table></div>",
  "<p class='note'>Regional sensitivity is in <a href='Trobisch_2022/results/full_atlas_region_gene_set_tests.csv'>the region-by-region results</a>. These tests are still within the same Trobisch dataset, not three separate datasets. Only 77 integrated neuronal nuclei came from spinal cord, and none of its donor-label groups met the 10-nucleus minimum. The significant combined results are therefore largely cortical.</p>",
  "<h3>Regional context within Trobisch</h3>",
  paste0("<p>All 20 neuronal labels appear in the coverage table. ", tested_new, " new-region subtype-by-region groups and ", tested_cx, " overlapping cortical groups had enough donors in both conditions for a gene-level test. Small groups remain visible but are unreliable; untestable groups are not treated as negative results.</p>"),
  "<figure><img src='Trobisch_2022/plots/integrated_neuron_de_counts.png' alt='Gene-level hits across new-region integrated neuronal groups'><figcaption>New-region integrated groups. Zero and untestable are distinguished in the table below.</figcaption></figure>",
  "<figure><img src='Trobisch_2022/plots/integrated_neuron_volcano.png' alt='New-region integrated neuronal volcano plots'><figcaption>New cerebellum/spinal-cord donors; global FDR across tested subtype-gene comparisons.</figcaption></figure>",
  "<figure><img src='Trobisch_2022/plots/integrated_neuron_ucell.png' alt='New-region UCell rank scores'><figcaption>UCell scores for every integrated neuronal label with new-region cells; each point is a donor. None of the subtype comparisons passed FDR.</figcaption></figure>")
if (file.exists("Trobisch_2022/plots/cortex_integrated_de_counts.png")) html <- c(html,
  "<figure><img src='Trobisch_2022/plots/cortex_integrated_de_counts.png' alt='Cortical integrated neuronal gene hits'><figcaption>Cortex is shown for label coverage, but these donors overlap Schirmer.</figcaption></figure>",
  "<figure><img src='Trobisch_2022/plots/cortex_integrated_volcano.png' alt='Cortical integrated neuronal volcano plots'><figcaption>Overlapping cortical samples; not independent replication of Schirmer.</figcaption></figure>")
if (file.exists("Trobisch_2022/plots/cortex_donor_qc.png")) html <- c(html,
  "<figure><img src='Trobisch_2022/plots/cortex_donor_qc.png' alt='Cortical donor quality measures'><figcaption>MS cortical nuclei tend to have higher detected genes and UMIs. This is a potential technical or tissue-composition confound for cell-level scores.</figcaption></figure>")
if (file.exists("Trobisch_2022/plots/cortex_integrated_ucell.png")) html <- c(html,
  "<figure><img src='Trobisch_2022/plots/cortex_broad_ucell.png' alt='Balanced cortical rank scores by donor'><figcaption>Broad cortical UCell scores: 18 labels represented in all eight overlapping donors receive equal weight.</figcaption></figure>",
  "<figure><img src='Trobisch_2022/plots/cortex_vip_ucell.png' alt='Cortical VIP-IN rank scores'><figcaption>VIP-IN, the strongest cortical neuronal subtype lead; each point is a donor shared with Schirmer. The table below retains every integrated label.</figcaption></figure>")
html <- c(html,
  "<h3>Every integrated neuronal label</h3>",
  "<div class='table-wrap'><table><thead><tr><th>Label</th><th>Cortex*</th><th>Cerebellum</th><th>Spinal cord</th></tr></thead><tbody>",
  coverage_rows, "</tbody></table></div>",
  "<p class='note'>*Trobisch cortex reuses Schirmer material. Its finer labels can refine the cortical analysis, but they do not provide an independent cohort.</p>",
  "<h2>Individual Rank 1 and Rank 2 genes</h2>",
  "<p>These older gene-level results are regional sensitivity analyses within Trobisch, not one combined gene-level model. Each cell shows the best tested subtype for that gene, with raw p and q corrected across the tested subtype-gene comparisons in that column's analysis. Green means q &lt; 0.05, whether the gene is MS-up (↑) or MS-down (↓). A dash means no usable test, not proof of zero expression.</p>",
  rank_table("Rank 1"), rank_table("Rank 2"),
  "<h2>Interpretation</h2>",
  "<p>Combined Trobisch GSEA, ROAST, and single-sample scores find MS-up rank-list signals, but most testable observations are cortical. The cerebellar GRAN-EN Rank 2 score signal depends on scoring method and only six donors; spinal cord is too sparse. The cortical cells overlap Schirmer and show higher detected genes and UMIs in MS, so this is a focused cortical lead, not independent validation or a verified pan-neuronal human-MS signature.</p>",
  "<p>Methods: raw counts were summed by donor and subtype; gene-level comparisons used DESeq2 for Schirmer and edgeR quasi-likelihood for Trobisch. Trobisch's primary gene-set test combines all usable donors with region in the model; region-specific results are sensitivity checks within that dataset. Rank-list scores used normalized pseudobulk expression and donor-level tests; UCell and camera were supporting checks. Sources: <a href='https://www.nature.com/articles/s41586-019-1404-z'>Schirmer et al. 2019</a> and <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC9547805/'>Trobisch et al. 2022</a>.</p>",
  "</body></html>")
writeLines(html, "scRNAseq_report.html")
