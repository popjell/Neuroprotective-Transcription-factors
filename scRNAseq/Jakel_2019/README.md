# Jakel 2019

This folder contains a simple neuronal single-nucleus RNA-seq analysis.

Run from this folder:

```bash
Rscript 01_qc_cluster.R
Rscript 03_prepare_celltypist.R
python3 04_celltypist.py
Rscript 05_plot_celltypist.R
Rscript 02_pseudobulk_deseq2.R
```

- `01_qc_cluster.R`: QC and standard Seurat clustering of the full atlas.
- `03_prepare_celltypist.R` and `04_celltypist.py`: export the sparse counts and annotate them with the adult human prefrontal-cortex CellTypist model.
- `05_plot_celltypist.R`: makes the broad cell-class and neuronal-subtype UMAPs.
- `02_pseudobulk_deseq2.R`: sums raw counts by donor and CellTypist neuronal subtype, runs DESeq2 for MS versus control, and makes the pseudobulk PCA, volcano, and Rank 1/2 bar chart.
- A neuronal subtype is tested only when it has at least 20 cells per pseudobulk and at least two control and two MS donors.
- `data/`: cached full-cell and neuronal count matrices plus the Rank 1/2 gene list. No download is needed.
- `results/`: pseudobulk counts and CSV result tables.
- `plots/`: requested figures.
