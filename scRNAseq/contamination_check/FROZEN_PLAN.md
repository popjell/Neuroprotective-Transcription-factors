# FROZEN PLAN: ambient RNA / glial contamination check of the Trobisch PSAP result

Written 2026-10-19 before any analysis was run. Parameters and rules below are copied from the task text and are not changed after seeing results.

## Frozen parameters (identical to trobisch.R unless stated)

- Donor mapping: CR_ctrl1 = CO39, CR_ctrl2 = CO14, CR_ctrl3 = CO36, CR_ctrl4 = CO25, CR_ms1 = MS60, CR_ms2 = MS229, CR_ms3 and CR_ms4 = MS371 (pooled).
- Eligibility: at least 20 cells per donor per subtype; fit a subtype only if at least 3 control and 3 MS donors remain.
- Pseudobulk: integer counts summed per donor within author subtype (cell_subtypes_intg).
- Filter: counts >= 10 in at least 2 pseudobulks.
- DESeq2: design ~condition, reference Control, contrast MS vs Control. padj is BH across the transcriptome within subtype.
- Glial marker panel (negative controls), fixed:
  - Microglia/macrophage: C1QA, C1QB, C1QC, CSF1R, P2RY12
  - Astrocyte: GFAP, AQP4, ALDH1L1, GJA1, SLC1A3
  - Oligodendrocyte: MBP, PLP1, MOBP, MOG, MAG
  - If any marker is in the 31-gene Rank 2 list, drop it and report that. Do not add replacements.
- set.seed(1010). OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1.

## Part A (existing results only)

1. For each fitted subtype, read results/<subtype>/deseq2.csv (read only) and extract the 15 markers. Report log2FC, padj, or "removed by count filter".
2. Rule: a subtype shows a contamination signal for a glial type if, among that type's markers that passed the filter, the median log2FC is > 0 AND at least 2 have padj < 0.05. Applied per subtype and per glial type.
3. Donor-level, descriptive only: per donor per fitted subtype, (a) glial marker CPM (sum of each glial type's markers / library size x 1e6) and (b) PSAP CPM. Report Spearman rho between microglial-marker CPM and PSAP CPM across donors per subtype. No p-values.

## Part B (full Trobisch object with non-neuronal cells)

1. Locate and report every candidate full cortex object (path, size, contents).
2. decontX: celda::decontX with default parameters, counts = raw integer counts of the full object, z = broad cell-type labels, batch = sample_id. No tuning.
3. Report contamination per cell, summarised per donor and per broad cell type (median and IQR); whether neurons from MS donors have higher estimated contamination than control (per-donor medians side by side).
4. Take decontaminated counts, round to integers, keep only neurons in cortex_neurons.rds (match on cell ID), and rerun the pipeline (eligibility, pseudobulk, filter, DESeq2).
5. Rule: PSAP survives in a subtype if decontaminated padj < 0.05 AND decontaminated log2FC >= 0.5 x the original log2FC. Report per subtype original log2FC and padj, decontaminated log2FC and padj, survives yes/no. Summary: number surviving out of the originally significant subtypes.
6. Rerun the Part A step 2 rule on the decontaminated results.

## Part C (availability check only; no download, no SoupX/CellBender/velocyto)

Report whether unfiltered/raw droplet matrices, BAM files or FASTQ files exist locally for Trobisch, with paths and sizes.

## Decisions made because the inputs named in the task are not on disk (recorded before any run)

1. trobisch.R, plot_trobisch.R, _keep_input/ and scRNAseq_old_20261016/ do not exist anywhere in WSL. The baseline is the existing scRNAseq/results/Cortex/<subtype>/deseq2.csv, which uses author labels `cell_subtypes` (not cell_subtypes_intg): 12 fitted subtypes, PSAP padj < 0.05 in 10 of 12 (not significant: L4-GABRG1-EN, SST-IN), log2FC positive in 12/12. "Subtype" below means `cell_subtypes`; the summary is out of these 10, not "10 of 11".
2. The DE pipeline is rebuilt from the frozen parameters above (min_cells 20, min_donors 3 per condition, CR_ms3+CR_ms4 pooled as MS371, integer sum pseudobulk, counts >= 10 in >= 2 pseudobulks, DESeq(~condition), MS vs Control, unshrunk log2FoldChange). Reproduction check: run on the uncorrected cortex_neurons.rds counts, it must give the same 12 fitted subtypes and PSAP log2FoldChange and padj within 1e-6 of results/Cortex. If not, stop and report; no tuning.
3. No broad-type column exists. decontX z is a fixed map of cell_subtypes: Neuron = all *-EN, *-IN, motoEN; Microglia = MG; Astrocyte = GM-AS, WM-AS; Oligodendrocyte = OL-A, OL-Stress, WM-OL; OPC = OPC; Vascular = END, PER; Immune = B cells, T cells. Any unmapped label stops the script.
4. The cell ID column of cortex_neurons.rds is `cell` (matched to `Cell` in cor/meta.tsv), not `cellId`.
5. The 31-gene Rank 2 list is rank_genes.csv with rank == "Rank 2" and use_primary_analysis == "Y". None of the 15 glial markers are in it, so none are dropped.
6. Layout follows the user's style guide: numbered scripts in contamination_check/ (not a single contamination_check.R), outputs in results/trobisch_cortex/ and figures/trobisch_cortex/ with the task's file names.
7. Figure panel G (violins of the lysosomal lipid-processing score and PSAP, Welch t on donor means) is an added request and plays no part in any rule above. The UMAP panels E and F were dropped at the user request; their files are in _removed_umap_figures/.
