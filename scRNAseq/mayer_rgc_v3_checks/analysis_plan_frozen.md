# FROZEN PLAN (mayer_rgc_v3_checks): pipeline validity checks. Do not edit after the first analysis is run.
Header for every output: unit = retina (Sample); 10 control, 12 MS; two MS donors contribute two retinas each; condition is confounded with tissue source (controls: Johns Hopkins; MS: Netherlands Brain Bank). All tests are Tier 3 and do not count toward the earlier families. These are sensitivity checks, not independent evidence (the authors derived these effects from the same data).

## PURPOSE
Test whether this pipeline can detect effects the authors report in the same data.

## 1. POSITIVE CONTROLS (frozen expectations, from Mayer et al.)
   a. CFH, pooled MS vs control: expected up. Run it in the baseline pooled DESeq2 and in the Analysis B subtype-adjusted model. Report log2FC, SE, p, padj.
   b. In the resilient group (author-defined: PG-ON, RGC26, RGC13), MS vs control, directional gene-set tests (ROAST/fry) for: induced stress/inflammatory genes JUN, JUND, FOSB, ACSL4, SQSTM1, NFKBIA, NFKBIZ, TGFB2 (expected up); neuronal-identity/excitability genes TBR1, POU4F1, EBF2, DLGAP1, SORCS1, GRIN2C, GPT2, PSAT1 (expected down).
   c. In the susceptible group (MG-ON4, MG-OFF2, ipRGC), MS vs control: DCN, SCUBE1, BMP5, MT-ATP8, MT-ND3, MT-ND4, MT-ND4L (expected up).
   d. For each set, report how many genes were testable, how many matched the expected direction, the set p, and the estimate. State plainly whether the pipeline reproduces the authors' directions. Use only genes present after filtering; list the missing ones.
   e. Do not tune anything to make controls pass. If they fail, say so and stop that section.

## 2. AMBIENT DIAGNOSTIC
   For each retina, in RGC nuclei: mean log-CPM and detection fraction of RHO, PDE6A, SAG, GNAT1 (rod), GFAP, GLUL, RLBP1 (glia), and P2RY12, CSF1R (microglia). Report by condition (donor-aware: retina-level means, exact permutation test). Report the correlation across the 22 retinas between the rod-marker score and pooled PSAP log-CPM. Also check the figshare readme or any text in _input/ for the strings "CellBender", "ambient", "SoupX", "corrected" and quote what is stated, or write "not stated".

## 3. DEPTH-ADJUSTED CELL-LEVEL MODELS
   Per-nucleus covariates: log total counts and number of genes detected. Models (lme4/lmerTest): PSAP log-CPM ~ condition + subtype + log_total_counts + (1|Sample); PSAP detection (binomial) ~ condition + subtype + log_total_counts + (1|Sample). Report the condition estimate, SE, p, with and without the depth covariate, and the per-retina mean total counts by condition (retina-level permutation test).

## 4. SENSITIVITY BOUND
   From the Tier 1 model in Analysis B, report the 95% and 90% intervals for the PSAP log2FC and the smallest effect detectable at 80% power. Add the number of retinas that would be needed for 80% power if the true effect were +0.14, +0.23 and +0.40 (assuming SE scales as 1/sqrt(n)).

## FINAL REPORT
one table of control results with a one-line verdict on whether the pipeline reproduces the authors' effects; the ambient table; the depth-adjusted estimates next to the unadjusted ones; the bound and the sample-size table; deviations. Do not interpret beyond the numbers.

---
# APPENDIX: interpretation decisions made BEFORE any analysis was run
1. Group fits for 1b and 1c reuse the Analysis D model exactly: limma-voom ~ group*condition on retina x author-group pseudobulks (resilient, susceptible; 27-label groups), block = retina, two-pass duplicateCorrelation, the same code as mayer_rgc_v2/05_analysis_D.R. The resilient MS effect is the conditionMS coefficient; the susceptible MS effect is the contrast conditionMS + interaction. The fit must reproduce the v2 D PSAP estimates (resilient 0.1388, susceptible 0.2926) or the script stops.
2. Gene-set tests: ROAST with 9999 rotations, seed 1010 before each set, block and correlation from the fit. The registered p is the DIRECTIONAL p (Up p for sets expected up, Down p for sets expected down). fry (two-sided, and one-sided in the expected direction) and the ROAST mixed and two-sided p are also reported. "Matched expected direction" = sign of the gene log2FC in the fitted model. The estimate is the mean log2FC of the testable genes (the median is also reported).
3. Genes are "testable" if they pass filterByExpr in the group-fit model (the same filter as Analysis D). Missing genes are listed with the reason (not in matrix, or removed by the filter).
4. 1a: baseline pooled DESeq2 values are read from mayer_rgc/results/all_RGC/deseq2.csv; Analysis B values from mayer_rgc_v2/results/B_all_genes.csv. Both are read-only; each is registered once as Tier 3 with notes.
5. Ambient: "mean log-CPM" = mean over a retina's nuclei of log1p(counts per 10,000); detection = fraction of the retina's nuclei with count > 0. The rod-marker score = mean of the retina-level mean log-CPM of RHO, PDE6A, SAG, GNAT1 (those present). Pooled PSAP log-CPM per retina is read from mayer_rgc_v2/results/A_psap_pooled_lcpm_per_retina.csv. Retina-level exact permutation test on the difference of means (all 646,646 labelings). Correlation: Pearson and Spearman (both registered).
6. Text search: all text files in mayer_rgc/_input/ (json, csv, tsv, txt, md); case-insensitive; quote matching lines; "not stated" otherwise.
7. Depth models: log_total_counts = log(nCount_RNA) (from the metadata), centered at its mean; n genes = nFeature_RNA, log-transformed and centered. Models: M0 (no depth covariate, equals mayer_rgc_v2 Analysis F), M1 (+ log_total_counts), M2 (+ log_total_counts + log_n_genes; an extra variant). PSAP log-CPM = log1p(counts per 10,000). Subtype = the 22-label celltypes. Per-retina mean total counts and mean n genes by condition: exact permutation test on retina means.
8. Sensitivity bound: normal-approximation intervals (estimate +/- 1.96 SE and +/- 1.645 SE) using the Analysis B PSAP log2FC and SE from mayer_rgc_v2/results/B_all_genes.csv; minimum detectable effect = 2.8 x SE; retinas needed n = 22 x (2.8 x SE / delta)^2, rounded up, for delta = 0.14, 0.23, 0.40. Because SE scales as 1/sqrt(n) this treats the current 22 retinas as the reference and ignores df changes.
9. Every hypothesis test is registered Tier 3 with register_test (the existing helper, sourced from mayer_rgc_v2/common_v2.R, writing to mayer_rgc_v3_checks/results/tests_registry.csv).
