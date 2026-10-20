# FROZEN ANALYSIS PLAN (mayer_rgc_v2). Do not edit after the first model is run.

Header for every output: unit = retina (Sample); 10 control, 12 MS; two MS donors contribute two retinas each (donor map not available unless _input/donor_map.csv exists); condition is confounded with tissue source (controls: Johns Hopkins; MS: Netherlands Brain Bank).

## TIERS (multiplicity control; frozen)
- Tier 1 (the single primary test, alpha 0.05, two-sided, pre-specified direction MS > control): PSAP condition effect in the subtype-adjusted pseudobulk model (Analysis B).
- Tier 2 (secondary family, BH across these tests only): (i) ITM2C in Analysis B; (ii) PSAP condition-by-group interaction (resilient vs susceptible) in Analysis D; (iii) PSAP in the composition-standardized pooled test (Analysis C); (iv) PSAP in Analysis D within the resilient group; (v) the 31-gene module (ROAST, mixed direction) in Analysis B.
- Tier 3 (exploratory, descriptive): everything else. Report estimates and intervals. No claims from Tier 3 results, no uncorrected p-value screening. Report how many Tier 3 tests were run.
The final report must state the total number of registered tests.

## STEP 1: INSPECT AND REPRODUCE THE BASELINE (checkpoint; STOP if any check fails)
1. Confirm the metadata columns: Sample, condition, celltypes (22 labels), cellsubtypes (27 labels). Record their exact label sets.
2. Map the six author-defined groups: susceptible = MG-ON4, MG-OFF2, ipRGC; resilient = PG-ON, RGC26, RGC13. Match by exact name after normalizing case and "-" vs "_". These need the 27-label column (cellsubtypes) because MG-ON4 and MG-OFF2 are merged in the 22-label column. Report the mapping. STOP if any of the six cannot be matched.
3. Rebuild the pooled retina pseudobulk (sum of raw integer counts per Sample), run DESeq2 exactly as the baseline (design ~condition, reference Control, gene filter counts >= 10 in >= 2 pseudobulks), and confirm PSAP log2FC matches the baseline results to 3 decimals (read the baseline file read-only). STOP if not.
4. Report the cells per Sample per label (27-label), the retinas by condition, and, if _input/donor_map.csv exists, the donor map. Otherwise state "no donor map; donor collapse skipped".

## ANALYSIS A: COMPOSITION DIAGNOSTICS (Tier 3)
- Per-retina proportions of each of the 27 subtypes; Mann-Whitney U (MS vs control) per subtype, with BH over the 27, to confirm the authors' composition shift. Register each test.
- Per-retina fraction of cells in the resilient group and in the susceptible group; report by condition.
- Variance decomposition for the pooled log-CPM PSAP value: how much between-retina variance is explained by condition, by resilient fraction, and by both (simple linear models, descriptive only). Also report the correlation of PSAP pooled log-CPM with resilient fraction across the 22 retinas.

## ANALYSIS B: SUBTYPE-ADJUSTED MODEL (TIER 1 AND PART OF TIER 2)
- Units: pseudobulk per retina x 22-label subtype (celltypes), summing raw counts. Eligible units: at least 20 cells. Keep a subtype only if it has eligible units in at least 3 control and 3 MS retinas.
- Gene filter: edgeR::filterByExpr with group = condition, using the eligible units.
- Model: limma-voom (voomWithQualityWeights), design ~ subtype + condition, with block = Sample and duplicateCorrelation to account for several subtypes per retina. Contrast: condition MS vs Control. Report log2FC, SE, t, p, padj (BH across the transcriptome) for all genes.
- Tier 1 = PSAP from this model. Tier 2 = ITM2C and the ROAST (or fry if ROAST fails) test on the 31-gene set, and also on the 30 genes excluding PSAP (Tier 3).
- Heterogeneity (Tier 3): subtype x condition interaction F-test for PSAP, and the per-subtype condition estimates for PSAP with intervals (descriptive; no per-subtype significance claims).
- Sensitivity (Tier 3), each as one extra row in the registry: (a) without sample quality weights; (b) excluding units from RGC4, RGC16, RGC23, RGC26, RGC13; (c) excluding retinas with fewer than 1,000 RGCs (S46, S136 plus any other); (d) leave-one-retina-out (22 refits, PSAP estimate only).

## ANALYSIS C: COMPOSITION-STANDARDIZED POOLED TEST (TIER 2)
- For each retina, compute log2 CPM per subtype from the subtype pseudobulks, then combine subtypes with FIXED weights equal to each subtype's mean proportion across all 22 retinas (renormalized over the eligible subtypes in that retina). This gives one composition-standardized expression value per retina per gene.
- Test condition with limma (design ~condition), and with an exact or 100,000-draw permutation test of the condition labels across retinas (seed 1010) for PSAP only. Report both.
- Also report the same estimate for ITM2C and the 31-gene set (Tier 3).

## ANALYSIS D: RESILIENT VS SUSCEPTIBLE GROUPS (TIER 2)
- Units: pseudobulk per retina x group, summing raw counts over the member subtypes (27-label column). Eligible: at least 20 cells per unit; fit a group only if at least 3 control and 3 MS retinas remain.
- Model: limma-voom, design ~group * condition, block = Sample, duplicateCorrelation. Pre-specified hypothesis (adaptive-response prediction): PSAP is more strongly increased in MS in the resilient group than in the susceptible group (interaction contrast, one-sided in this direction, stated in the registry).
- Also run DESeq2 per group (as the authors did): design ~condition, resilient and susceptible separately.
- Sensitivity (Tier 3): resilient group excluding RGC13 and RGC26 (i.e. PG-ON only), if eligible; susceptible group per member subtype.

## ANALYSIS E: SPECIFICITY AND MODULE NULLS (TIER 3)
- For Analysis B and Analysis C: compare the PSAP t-statistic and log2FC with 10,000 random gene sets of size 1, drawn from genes matched on mean log-CPM (deciles), and report the percentile of PSAP. Do the same for the median log2FC of the 31-gene set against 10,000 size-matched, expression-matched random sets.
- Comparator lysosomal genes (frozen list, use those present): LAMP1, LAMP2, CTSB, CTSD, CTSZ, GBA1, HEXA, HEXB, GLB1, NPC1, NPC2, GRN, SCARB2, CD63, TMEM106B, ASAH1, SMPD1, GALC, LIPA. Report PSAP's rank by log2FC and by t among these.
- Transcriptome-wide shift: median and mean log2FC of ALL tested genes in Analysis B and C, with a bootstrap interval across genes, to say whether MS retinas are globally shifted (which would give the "28 of 30 positive" pattern a technical explanation).

## ANALYSIS F: CELL-LEVEL, SUBTYPE-ADJUSTED (TIER 3, SUPPORTIVE)
- Reuse the baseline UCell scores if saved in scRNAseq/mayer_rgc/results/ (read-only). If not saved, recompute with UCell for the 15 module genes passing the 5% detection filter, with and without PSAP.
- Models: score ~ condition + subtype + (1|Sample) with lme4/lmerTest, for log-CPM PSAP, PSAP detection (binomial), module score, leave-PSAP-out score. Report the condition estimate, SE and p for each, with n retinas.

## ANALYSIS G: METHOD CONCORDANCE (TIER 3)
For PSAP, ITM2C, CNDP2, report the pooled condition log2FC and SE from: baseline DESeq2; edgeR quasi-likelihood (design ~condition); limma-trend on pooled log-CPM; and the Analysis B and C estimates. One concordance table.

## ANALYSIS H: POWER AND META-ANALYSIS INPUTS
- For each of the Analysis B, C and D PSAP estimates, report SE and the minimum detectable absolute log2FC at 80% power, alpha 0.05 (approximately 2.8 x SE). State how that compares with the Trobisch PSAP range 0.80 to 1.59 (do not recompute Trobisch).
- Write results/meta_inputs.csv: analysis, gene (PSAP, ITM2C), log2FC, SE, n_control, n_ms, unit (retina), for the pooled DESeq2, B, C and D estimates, ready for a later random-effects meta-analysis. No meta-analysis in this task.

## ANALYSIS I: AUTHOR MODULE LOOKUP (read-only)
Search the files already in _input/ (including article.json) for the authors' module A to F gene lists. If present, report which of the 31 genes and PSAP fall in each module. If not present, write "module lists not available in the deposited files" and do NOT guess.

## OPTIONAL ANALYSIS J: DONOR COLLAPSE
Only if _input/donor_map.csv exists. Collapse the retinas of the two doubled MS donors by summing counts, then repeat the pooled DESeq2 and Analysis B (20 units). Otherwise skip.

## VALIDATION (before reporting)
- Rerun the Tier 1 model and the Analysis C permutation test from scratch into a temp name inside mayer_rgc_v2/ and confirm identical output (checksums). Never delete the temp folder; list it for me to remove.
- Confirm that all registry rows have a tier, and that Tier 2 padj is BH across exactly the Tier 2 tests.
- Confirm the sha256 check of the baseline folder (no change).
- Write sessionInfo.txt and run.log (every filter count, every skipped analysis and why).

---
# APPENDIX: interpretation decisions made BEFORE any model was run (the text above is verbatim; ambiguities resolved here)
1. Tier 1 and Tier 2 p-values are two-sided, EXCEPT Tier 2 (ii) (interaction), which is one-sided in the direction resilient > susceptible as the plan states. Tier 1 additionally reports the one-sided (MS > control) p in the notes column, for information only.
2. Tier 2 (iii) uses the limma p-value from Analysis C. The permutation p-value for PSAP is registered separately as Tier 3.
3. Tier 2 (iv) is the PSAP condition effect (MS vs Control) within the resilient group from the Analysis D interaction model (two-sided).
4. Tier 2 (v) is the ROAST mixed-direction p (9,999 rotations, seed 1010) on the 31-gene set restricted to genes that pass the filterByExpr gene filter in Analysis B.
5. voom: TMM normalization; two-pass duplicateCorrelation (voomWithQualityWeights -> duplicateCorrelation -> voomWithQualityWeights with block/correlation -> duplicateCorrelation -> lmFit). Same recipe in B and D.
6. Log-CPM in Analysis C, E, G: edgeR cpm(log = TRUE, prior.count = 1) with TMM-normalized library sizes.
7. Analysis E: decile of mean log-CPM (Analysis B average log-CPM; Analysis C mean of the retina values). Random single genes are drawn with replacement from the same decile excluding PSAP. Percentile = fraction of random draws with value <= PSAP value; one-sided empirical p for "PSAP at or above" = fraction >= PSAP value. Gene-set null: each of the 31 genes (testable) is replaced by a random gene from the same decile; statistic = median log2FC. Transcriptome-wide bootstrap: 2,000 resamples of genes.
8. Analysis F: log-CPM PSAP = log1p(counts per 10,000) (as in the baseline); detection by glmer binomial (default Laplace); subtype = 22-label celltypes; module = the 15 genes passing 5% detection in the baseline module_genes.csv.
9. Analysis G: edgeR QL SE derived as |logFC| / sqrt(F) (1 df contrast). The baseline DESeq2 row is read from the baseline files and is not a new test.
10. Analysis H meta inputs include the D estimates for the resilient group, the susceptible group, and the interaction.
11. Every new hypothesis test is registered; descriptive items (ranks, transcriptome-wide shift, variance fractions, per-subtype PSAP estimates, minimum detectable effects) are not tests.
12. Leave-one-out row (B sensitivity d) registers the median PSAP estimate across the 22 refits (se and p NA) with the min and max in the notes.
