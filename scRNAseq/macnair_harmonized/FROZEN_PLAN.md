# FROZEN PLAN: Macnair re-analysis of PSAP with harmonized labels

Written after Stage 1 (inputs, label transfer, eligibility tables; no PSAP or Rank 2 result seen) and before any Stage 2 test. Not to be changed after PSAP results.
Run with set.seed(1010), OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1. No age, sex, PMI or brain-bank terms in any model.

## Inputs and Stage 1 decisions
- Macnair et al., Zenodo 8338963; grey matter (control GM, NAGM, GML) neurons (type_broad Excitatory + Inhibitory); 175,237 neurons, 64 samples, 35 donors.
- Labels: SingleR transfer of the 18 Trobisch neuronal subtypes (31 Rank 2 genes removed from reference and query), coarse groups upper_EN / deep_EN / IN by the Schirmer name rule. 1.6% unassigned.
- Exclusions: none (no donor or sample ID matches the 7 shared Trobisch/Schirmer donors or the 19 Schirmer donors; IDs are anonymised, so overlap cannot be ruled out by ID). `results/excluded_donors.csv` is empty.
- Eligible donors: at least 20 pooled upper_EN cells. 34 donors: 14 control, 20 MS (Ind83, a control, has 19 and is out).
- Rule A (bank adjustment) is dropped by decision of the user. Control donors are all Edinburgh BB and MS donors all UK MS TB, so MS vs control is confounded with brain bank and the primary cannot separate them.

## Primary: PSAP in upper_EN, MS vs Control
- One pseudobulk per donor: all of the donor's GM samples summed (GML and NAGM together for MS). 34 donors, 14 control vs 20 MS.
- Gene filter: counts >= 10 in at least 2 pseudobulks.
- DESeq2 ~condition, reference Control, contrast MS vs Control.
- Report PSAP log2FC, 95% interval (1.96 x lfcSE) and two-sided Wald p. Single pre-specified test, p < 0.05. Report direction regardless.
- Note carried in this file and in psap.csv: control and MS donors come from different brain banks (Edinburgh BB vs UK MS TB); the comparison is confounded with bank.

## Secondary 1: pooled limma, all neurons
- As schirmer_harmonized/03_pooled_limma.R: pseudobulk per donor x coarse group (upper_EN, deep_EN, IN) with at least 20 cells, filterByExpr (group = condition), TMM, voom.
- Design ~ coarse_group + condition (reference Control), block = donor via duplicateCorrelation, two rounds.
- Report PSAP logFC, 95% CI, p, and the consensus correlation. Same bank caveat.

## Secondary 2: within-donor GML vs NAGM (local pathology, not MS vs control)
- The 10 MS donors with at least one GML and one NAGM sample that each have at least 20 upper_EN cells: Ind66, Ind67, Ind68, Ind70, Ind75, Ind80, Ind81, Ind93, Ind96, Ind99. Need at least 3.
- Samples with at least 20 upper_EN cells are summed within donor x tissue, so each donor contributes one GML and one NAGM pseudobulk.
- Gene filter as above. DESeq2 ~ donor + tissue, reference NAGM, contrast GML vs NAGM.
- Report PSAP log2FC, 95% interval, p. samples.csv lists every summed sample with its seq_pool.

## Secondary 3: tissue gradient (control GM, NAGM, GML)
- One upper_EN pseudobulk per donor x tissue (samples summed within donor x tissue); at least 20 cells per pseudobulk. Controls give a control GM pseudobulk; MS donors give NAGM and/or GML.
- limma-voom after filterByExpr (group = tissue) and TMM; design ~tissue (reference control GM); block = donor via duplicateCorrelation, two rounds.
- Contrasts: NAGM vs control GM, GML vs control GM, GML vs NAGM. Report PSAP logFC, 95% CI, p for each. Control GM contrasts carry the same bank caveat.

## Descriptive only
- The 31 Rank 2 genes from the primary fit, BH within the 31.
- Per-donor normalized upper_EN PSAP (primary-fit size factors) with condition and tissue composition: descriptive/donor_psap.csv.

## Rules
No additional subsets, covariates, thresholds or methods. A null is reported as a null.
