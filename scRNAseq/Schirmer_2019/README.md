# Schirmer 2019

This folder contains a simple neuronal single-nucleus RNA-seq analysis.

Run from this folder:

```bash
Rscript 01_qc_cluster.R
Rscript 02_pseudobulk_deseq2.R
Rscript 03_rank_gene_plot.R
Rscript 04_plot_subtype_umap.R
Rscript 05_cell_level_markers.R
Rscript 06_ucell_signature_scores.R
Rscript 07_rank_set_tests.R
Rscript 08_donor_rank_scores.R
# From the parent scRNAseq folder:
Rscript single_sample_rank_scores.R Schirmer
```

- `01_qc_cluster.R`: QC, standard Seurat clustering, and a UMAP labelled with the same neuronal subtypes used for pseudobulk.
- `02_pseudobulk_deseq2.R`: sums raw counts per patient and supplied neuronal subtype, runs DESeq2 for MS versus control, and makes the pseudobulk PCA, volcano, and Rank 1/2 bar chart.
- `03_rank_gene_plot.R`: plots the significant Rank 1 and Rank 2 genes and their effect sizes.
- `06_ucell_signature_scores.R`: uses UCell to score Rank 1 and Rank 2 from the sparse count matrix, averages cell scores per patient and subtype, then compares donor-level scores between MS and control. The supplied lists do not specify up/down direction.
- `07_rank_set_tests.R`: tests the two lists using pseudobulk preranked GSEA, limma `camera` and `roast`, and patient-level mean logCPM scores. Results are saved in `results/rank_set_method_comparison.csv`.
- `08_donor_rank_scores.R`: compares mean normalized Rank 1 and Rank 2 expression per patient and subtype, then gives each subtype equal weight in a broad neuronal comparison. The results are in `results/donor_rank_tests.csv` and `results/broad_neuron_rank_tests.csv`.
- `single_sample_rank_scores.R`: compares seven Rank 1/2 scores across all nine neuronal subtypes using one pseudobulk profile per patient and subtype. It reuses cached counts and saves `results/single_sample_rank_tests.csv`, individual scores, and a plot. The script needs GSVA, singscore, and AUCell (already installed in the project-local `.r_lib`).
- `05_cell_level_markers.R`: exploratory Seurat `FindMarkers` comparison of MS and control cells within each neuronal subtype. It uses Seurat's default Wilcoxon test and Bonferroni-adjusted p values; cells are not independent donors, so use the patient-level pseudobulk results as the main evidence.
- `data/`: cached raw neuronal counts and the Rank 1/2 gene list. No download is needed.
- `results/`: pseudobulk counts and CSV result tables.
- `plots/`: requested figures.

The condition-derived layer 2/3 labels are combined as one neutral `EN-L2-3` subtype before pseudobulk.
