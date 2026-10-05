# Trobisch 2022

This folder analyzes all 20 integrated neuronal labels in Trobisch's atlas as one dataset. Script 18 is the combined gene-set analysis: it uses all usable donor-subtype samples and adjusts for brain region. The same script also writes region-by-region sensitivity results. The cached author Seurat objects contain raw RNA counts. The Cell Browser expression matrix is SCT-corrected and must not be used for pseudobulk DE.

For the integrated-label analysis, run from this folder:

```bash
Rscript 05_author_full_atlas_umap.R
Rscript 09_integrated_neurons.R
Rscript 10_integrated_ucell.R
Rscript 12_cortex_integrated_neurons.R
Rscript 13_integrated_camera.R
Rscript 14_cortex_ucell.R
Rscript 15_cortex_sensitivity.R
Rscript 16_cortex_gene_drivers.R
Rscript 17_cortex_qc.R
Rscript 18_full_atlas_gene_sets.R
# From the parent scRNAseq folder:
Rscript single_sample_rank_scores.R Trobisch
```

Run `prlimit --as=6500000000 -- Rscript 11_extract_cortex.R` once, before script 12, only if `data/cortex_neurons.rds` is missing. The compact cortical cache is now present, so normal reruns skip this large-object extraction. If `data/new_regions.rds` is missing, run `Rscript 01_extract_raw_counts.R` once. None of the analysis scripts redownload data. The main inputs are:

- `data/meta.tsv` and `data/UMAP.coords.tsv.gz` from [the integrated Cell Browser atlas](https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional/integrated/)
- `data/cer_seurat.rds` from [the cerebellum atlas](https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional/cer/seurat.rds)
- `data/sc_seurat.rds` from [the spinal-cord atlas](https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional/sc/seurat.rds)
- `data/cor_seurat.rds` from [the cortical atlas](https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional/cor/seurat.rds)

The full-atlas UMAP uses the authors' 38 integrated cell labels. Scripts 03–08 use six *regional* neuronal labels and are retained as a separate, older analysis. Scripts 09–18 use the 20 *integrated* neuronal labels. Script 18 joins the existing cortical and new-region pseudobulk caches, filters donor-subtype samples with fewer than 10 nuclei, retains every label in the output, uses TMM/voom, ranks genes by the region-adjusted MS-control moderated t statistic, and runs preranked GSEA plus ROAST. A gene passes the expression filter if it passes in any adequately sampled region. BH q-values are calculated across all testable subtype/list comparisons within each method. Groups without enough donors or expressed list genes are marked untestable. The 9,675 cortical neuronal raw-count columns match the integrated atlas cell IDs and count totals. Cortex overlaps Schirmer, so Trobisch is not an independent replication cohort. The older cell-level `FindMarkers` p-values treat nuclei as independent and are exploratory.

The integrated-label outputs are `results/integrated_neuron_coverage.csv`, `results/integrated_neuron_test_status.csv`, `results/integrated_neuron_de.csv`, `results/integrated_neuron_rank_sets.csv`, `results/integrated_neuron_ucell_tests.csv`, and the matching `cortex_integrated_*` tables. The cortical object is extracted only once; the smaller `data/cortex_neurons.rds` is reused. Very small groups are shown, not hidden, but their estimates are unstable. Run `Rscript make_scrnaseq_report.R` from the parent folder to refresh the HTML report.

The combined gene-set output is `results/full_atlas_gene_set_tests.csv`; `results/full_atlas_region_gene_set_tests.csv` is a within-dataset regional sensitivity check. In the regional table, `q_method` corrects across all regional tests of one method, while `q_region` corrects within one region and method. `plots/full_atlas_gene_sets.png` shows the combined tests. In the integrated annotations, spinal-cord neurons are too sparse for donor-subtype gene-set testing after the 10-nucleus minimum; spinal labels still appear as untestable rather than disappearing.

`single_sample_rank_scores.R` scores Rank 1/2 by ssGSEA, GSVA, singscore, AUCell, PLAGE, z-score, and gene-wise 0–1 min-max averaging. It joins the cached cortex and new-region donor pseudobulks, keeps all 20 neuronal labels in `results/single_sample_rank_tests.csv`, and saves the donor scores and `plots/single_sample_rank_tests.png`. The primary tests adjust for region where applicable; `results/single_sample_rank_regional_tests.csv` checks the same scores within each region. q-values are corrected within method across testable label/list comparisons; scoring methods are not independent cohorts. The project-local `.r_lib` contains the needed extra packages.

The older six-regional-label outputs are `results/deseq2_results.csv`, `results/rank_gene_results.csv`, `results/cell_level_markers.csv`, and `results/granule_sensitivity.csv`.

`07_ucell_signature_scores.R` scores the Rank 1 and Rank 2 gene sets with UCell in the six author-labelled neuronal subtypes from the new cerebellum and spinal-cord samples. Scores are averaged per donor before MS-control tests. The gene names used to restore count-matrix row labels are cached in `data/gene_names.txt` and are regenerated by `01_extract_raw_counts.R`.

`08_rank_set_tests.R` runs pseudobulk preranked GSEA, limma `camera` and `roast`, and donor-level mean-expression tests for both lists. Tests with fewer than three list genes passing the expression filter are marked not testable. Results are saved in `results/rank_set_method_comparison.csv`.
