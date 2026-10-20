# FROZEN PLAN: Schirmer replication of PSAP with harmonized cell types and pooled neurons

Written after steps 1-2 (inputs, label transfer; no PSAP result seen) and before steps 3-5. Not to be changed after PSAP results.
Run with set.seed(1010), OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.

## Documented deviations from the task text (all decided with the user before any analysis)
1. **No independent-donor exclusion.** The user chose to use all 19 Schirmer donors as replication and separate shared from independent donors later.
   The 7 donors also present in Trobisch are listed below and flagged in `results/shared_donors.csv` and `results/label_transfer_counts.csv`.
   Results are therefore NOT independent of the Trobisch discovery cohort.
2. **Shared-donor list source.** `scRNAseq_old_20261016/schirmer_slimmed.R` does not exist. The list is the Trobisch sample -> donor mapping in `scRNAseq/cell_death/00_prepare_cortex.R`.
3. **Reference.** `_keep_input/cortex_neurons.rds` and the label column `cell_subtypes_intg` do not exist. The reference is `ms-cross-regional/cor/seurat.rds` + `meta.tsv`
   (RNA counts), label column `cell_subtypes`, restricted to the 18 neuronal labels that exist (not 20).
4. **Trobisch comparator** uses the same cells and the donor mapping from item 2 (trobisch.R does not exist).
5. **Rank 2 genes** are read from `Neuroprotective-Transcription-factors/rank_genes.csv`; it already lists TBC1D13 (the TBC13 -> TBC1D13 mapping is a no-op).
6. Code is split into numbered scripts with `common.R` (user style guide) instead of a single `schirmer_harmonized.R`.
7. SingleR with `de.method = "wilcox"` needs scran; Bioconductor 3.20 packages (matching the installed DESeq2 1.46.0) were installed into `schirmer_harmonized/Rlib`.

## Shared donors (in both Schirmer and Trobisch; flagged, not excluded)
CO39, CO14, CO36, CO25 (controls); MS60, MS229, MS371 (MS).
Trobisch sample -> donor: CR_ctrl1=CO39, CR_ctrl2=CO14, CR_ctrl3=CO36, CR_ctrl4=CO25, CR_ms1=MS60, CR_ms2=MS229, CR_ms3=MS371, CR_ms4=MS371.

## Inputs
- Schirmer 2019: `scRNAseq/data/Schirmer_2019/raw/matrix.mtx` (raw integer counts), meta.tsv, sample_donor.csv. 21 samples -> 19 donors; samples from one donor are summed (MS371 = MS6+MS7, MS377 = MS2+MS8).
- Neurons by the authors' labels: EN-L2-3-A, EN-L2-3-B, EN-L4, EN-L5-6, EN-MIX, EN-PYR, IN-PVALB, IN-SST, IN-SV2C, IN-VIP (23,705 cells).

## Step 2: harmonize cell types (label transfer from Trobisch) - frozen
- Method: SingleR with the reference aggregated to pseudobulk per label (aggregateReference, default settings), log-normalized counts, de.method = "wilcox", default pruning. Pruned (NA) labels are "unassigned" and excluded.
- Coarse groups, assigned by this name rule only:
  - upper_EN: label ends in "-EN" and its first layer token is L2 or L3 (e.g. L2-EN, L2-L3-EN, L3-EN)
  - deep_EN: any other label ending in "-EN"
  - IN: label ends in "-IN"
  - anything else: excluded
- Assignment (`results/coarse_group_map.csv`):
  - upper_EN: L2-EN, L2-L3-EN, L3-EN
  - deep_EN: L4-GABRG1-EN, L4-LRRK1-EN, L4-MME-EN, L4-NPFFR2-EN, L5-L6-ADAMTSL1-EN, L5-L6-FRMD6-EN, L5-L6-HTR2C-EN, L5-L6-NTNG2-EN
  - IN: KIT-IN, PVALB-CA8-IN, PVALB-IN, RELN-IN, SST-IN, VIP-IN
  - excluded: motoEN
- Result (seen before freezing): 446 of 23,705 Schirmer neurons unassigned (1.9%); upper_EN 7,408; deep_EN 7,006; IN 6,875; excluded 1,970.

## Step 3: primary test (one pre-specified test)
- Unit: upper_EN, Schirmer donors (all 19, see deviation 1).
- Pseudobulk: integer counts summed per donor across all upper_EN cells.
- Eligibility: at least 20 upper_EN cells per donor; require at least 3 control and 3 MS donors, else report "not testable" and stop the primary.
- Filter: counts >= 10 in at least 2 pseudobulks.
- DESeq2: design ~condition, reference Control, contrast MS vs Control.
- Primary endpoint: PSAP log2FC, its 95% interval (lfcSE x 1.96), and the two-sided Wald p-value. Single pre-specified test, so p < 0.05 is the threshold with no further correction. Report direction regardless of significance.
- Comparator: the identical upper_EN pooled model on Trobisch (author labels, same donor mapping), reported with the same three numbers, labeled "Trobisch discovery cohort, for effect-size comparison only".

## Step 4: secondary pooled model (all neurons)
- Pseudobulk per donor x coarse group (upper_EN, deep_EN, IN), eligibility 20 cells per donor per group; drop ineligible donor-group pseudobulks.
- Filter: edgeR::filterByExpr with group = condition, then TMM normalization.
- limma-voom with design ~ coarse_group + condition, block = donor via duplicateCorrelation (two rounds, as in the limma user guide). Test the condition coefficient (MS vs Control).
- Report PSAP logFC, 95% CI, p, and the consensus within-donor correlation. Run the same model on Trobisch as comparator.

## Step 5: secondary descriptive outputs
- Per Trobisch fine label with enough Schirmer cells (same 20-cell / 3-and-3 rule), run the trobisch.R-style DESeq2 and report PSAP only. Descriptive; not part of the primary claim.
- For the 31 Rank 2 genes (human_symbol from rank_genes.csv, rank == "Rank 2", use_primary_analysis == "Y", TBC13 -> TBC1D13), report the upper_EN primary-model results with BH applied within the 31 genes. Descriptive.

## Rules
No additional subtypes, covariates, thresholds or methods. If the primary is not significant, report it as not significant; do not search for a configuration where it is.
