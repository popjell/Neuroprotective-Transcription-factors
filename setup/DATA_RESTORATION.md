# Restoring scRNA-seq inputs on WSL

From the repository root, run:

```bash
bash setup/restore_scrnaseq.sh
```

The restoration runs sequentially to limit memory use. It resumes partial downloads, skips existing input/cache files, and does not rerun the statistical analyses or overwrite result CSVs/HTML reports. Downloads and reconstructed data are Git-ignored.

Sources:

- Trobisch: author RNA-count Seurat objects at https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional/{cer,sc,cor}/seurat.rds; integrated coordinates at the same site's integrated/UMAP.coords.tsv.gz. The tracked integrated metadata and new-region metadata are reused. Existing extraction scripts verify cell IDs and published count totals.
- Schirmer: https://cells.ucsc.edu/ms/rawMatrix.zip and https://cells.ucsc.edu/ms/meta.tsv. The archive contains the raw MatrixMarket counts, gene map, and barcodes. The browser expression TSV is normalized and is not used. Neurons are selected by EN-/IN- author labels; browser sample identifiers supply patient IDs. Published QC values remain unchanged. `data/restored_count_totals.csv` records differences between archive totals and browser UMIs; a mismatch is not silently treated as agreement.
- Jakel: GSE118257 expression and annotation files from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/. The dense TSV is converted one row at a time into sparse integer counts and matched to the author annotation by cell ID. The released gene matrix omits some genes counted by the original QC annotation. Detected-gene totals are compared to the annotation and recorded in `data/restored_detected_genes.csv`; the analysis scripts recompute their QC fields from the count matrix.
- Fagiani: GSE277430 metadata and normalized sparse matrix from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277430/suppl/. GEO calls this matrix normalized expression. The outer gzip wrapper is removed before R reads the compressed RDS. Count reconstruction attempts `expm1(value) * nCount_RNA / 10000`, and aborts unless every value is within 1e-5 of an integer and every reconstructed neuron has exactly the published UMI and detected-gene totals. The neuronal cache uses the published QC thresholds and the cell list in the committed prior-analysis metadata, retaining the same analysis selection. A successful reconstruction records its provenance beside the cache. This is a mathematically reconstructed count cache, not an independent raw-count download.
- Human NicheNet inputs: https://zenodo.org/records/7074291. Both downloaded files are checked against the record's MD5 checksums.

Raw-cache formats follow the current analysis scripts. Existing reports can be inspected without rebuilding downstream results. The restore command also recreates missing donor/subtype pseudobulk RDS files using the existing analysis scripts, stopping before statistical tests. To refresh the combined tests, activate `setup/activate.sh`, run scripts 09 and 12 from `scRNAseq/Trobisch_2022`, then follow that directory's README. Schirmer/Jakel/Fagiani run orders remain in their READMEs.

Restoration and verification logs are under `setup/logs/`. The studies' older statements that no download is needed describe the previous laptop's caches; use this recipe on a fresh clone.
