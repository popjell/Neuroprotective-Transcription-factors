# Schirmer 2019

This folder contains a simple neuronal single-nucleus RNA-seq analysis.

Run from this folder:

```bash
Rscript 01_qc_cluster.R
Rscript 02_pseudobulk_deseq2.R
```

- `01_qc_cluster.R`: QC, standard Seurat clustering, and a UMAP labelled with the same neuronal subtypes used for pseudobulk.
- `02_pseudobulk_deseq2.R`: sums raw counts per patient and supplied neuronal subtype, runs DESeq2 for MS versus control, and makes the pseudobulk PCA, volcano, and Rank 1/2 bar chart.
- `data/`: cached raw neuronal counts and the Rank 1/2 gene list. No download is needed.
- `results/`: pseudobulk counts and CSV result tables.
- `plots/`: requested figures.

The condition-derived layer 2/3 labels are combined as one neutral `EN-L2-3` subtype before pseudobulk.
