# Jakel 2019

This folder contains a simple neuronal single-nucleus RNA-seq analysis.

Run from this folder:

```bash
Rscript 01_qc_cluster.R
Rscript 02_pseudobulk_deseq2.R
```

- `01_qc_cluster.R`: QC and standard Seurat clustering of all cell types, followed by separate reclustering of the neuronal cells.
- `02_pseudobulk_deseq2.R`: uses the new neuronal Seurat clusters, sums raw counts per patient and cluster, runs DESeq2 for MS versus control, and makes the pseudobulk PCA, volcano, and Rank 1/2 bar chart.
- A neuronal cluster is tested only when it has at least 20 cells per pseudobulk and at least two control and two MS patients.
- `data/`: cached full-cell and neuronal count matrices plus the Rank 1/2 gene list. No download is needed.
- `results/`: pseudobulk counts and CSV result tables.
- `plots/`: requested figures.
