# Fagiani 2025

This folder contains a simple neuronal single-nucleus RNA-seq analysis.

Run from this folder:

```bash
Rscript 01_qc_cluster.R
Rscript 02_pseudobulk_deseq2.R
```

- `01_qc_cluster.R`: QC and a full metadata table using the published neuronal subtype labels. It does not run the memory-heavy neuron-only UMAP.
- `02_pseudobulk_deseq2.R`: sums raw counts per patient and published neuronal subtype, runs DESeq2 for MS versus control, and makes the pseudobulk PCA, volcano, and Rank 1/2 bar chart.
- A neuronal cluster is tested only when it has at least 20 cells per pseudobulk and at least two control and two MS patients.
- `data/`: cached raw neuronal counts, metadata, and the Rank 1/2 gene list. No download is needed.
- `results/`: pseudobulk counts and CSV result tables.
- `plots/`: requested figures.
