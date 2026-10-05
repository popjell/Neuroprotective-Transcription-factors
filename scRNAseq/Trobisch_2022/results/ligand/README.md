# Exploratory ligand analysis of the Trobisch atlas

This analysis asks which microglial or astrocytic ligands could explain MS-up neuronal expression. It uses all three Trobisch regions and all 20 author-labelled neuronal subtypes; a region–subtype is scored only when it has enough donor samples and MS-up genes. Ligand activity is a model-based ranking, **not** a statistical test or evidence of physical signaling.

## Main results

- The strict analysis used donor-level pseudobulk MS-vs-control results and required at least 20 FDR-significant MS-up target genes. Sixteen cortical neuronal subtypes passed. No cerebellar or spinal subtype passed; see [subtype_status.csv](subtype_status.csv) for every region–subtype and its reason.
- In cortical L2-L3-EN, the leading glia-enriched candidates for the **broader** MS-up gene program were ANOS1, MFGE8, LGALS3, HLA-DRA, and SDC4. The strongest plausible pairings include astrocytic ANOS1–neuronal SDC2, astrocytic MFGE8–neuronal ITGAV, and astrocytic LGALS3–neuronal PTPRK. These are predictions from expression plus the reference network, not validated interactions. [All ligand rankings](ligand_activity.csv) and [ligand–receptor support](ligand_receptor_support.csv) are provided.
- The leading ligands do **not** provide a compelling direct explanation for the Rank 1/2 genes. In L2-L3-EN, PSAP and ITM2C are MS-up at FDR < 0.05, but neither appears among the top 250 predicted targets of the top 20 ligand candidates. [All significant rank-gene links](significant_rank_target_links.csv) show the target-specific scores and ranks.
- Among all eligible L2-L3 ligands, IL33 has the strongest predicted link to ITM2C (target rank 283), but the required IL33 primary receptor IL1RL1/ST2 was undetected in these neuronal nuclei; IL1RAP alone is not sufficient evidence for signaling. SPP1 has the strongest link to PSAP (target rank 777), which is not a strong target-specific prediction. These are not confirmed mechanisms.
- In the glial donor check, median detection of ANOS1, MFGE8, LGALS3, and microglial HLA-DRA was higher in MS than control; SDC4 was not. These are descriptive comparisons with four donors per condition, not independent confirmation.
- A separate, **nominal-p-only** sensitivity run scored cerebellar GRAN-EN and VIP-IN because neither had 20 FDR-significant MS-up genes. It yielded ANOS1 among its candidates, but these cerebellar rankings are less reliable and should not be combined with the strict cortical rankings. [Sensitivity results](cerebellum_nominal_activity.csv).

## Figures

![Top cortical ligand candidates](top_ligands.png)

![Donor-level detection of leading L2-L3 candidates](top_ligand_donors.png)

![Exploratory cerebellar candidates](cerebellum_nominal_ligands.png)

## Method and limits

The strict analysis uses the existing author-annotation pseudobulk comparisons. Potential ligands must be detected in at least 10% of MS glial nuclei and at least two MS donors, and be at least twice as detectable in glia as in pooled neurons from the same region. Receptors must be detected in at least 5% of MS nuclei within the receiving subtype. NicheNet then ranks ligands for the FDR-significant MS-up neuronal genes. Rank 1/2 target links are checked separately. The cerebellar sensitivity run substitutes nominal p < 0.05 for the FDR target requirement and is explicitly exploratory.

Nuclear RNA detection does not prove secretion or surface protein expression. The atlas does not establish physical proximity, ligand release, receptor-complex assembly, or causality. Multiple cortical subtypes share the same donors, so their similar ligand rankings are not independent replication. NicheNet's corrected AUPR is a ranking score, not a p-value or FDR.

The expression filter excludes lowly detected cytokines such as TNF and IL1B in these nuclei. Their exclusion is a limitation of this RNA-based screen, not evidence that their proteins are absent in MS tissue.

To rerun from `Trobisch_2022/`, run `Rscript 19_extract_cortex_senders.R` once to cache the cortical sender counts, then `Rscript 20_ligand_analysis.R`, `Rscript 21_ligand_donor_check.R`, `Rscript 22_rank_target_links.R`, and optionally `Rscript 23_cerebellum_sensitivity.R`. The two NicheNet reference files are cached in `data/` and are not downloaded again.

Reference: [NicheNet analysis guide](https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper.md); [NicheNet-v2 human reference files](https://zenodo.org/records/7074291); [IL33 receptor-complex study](https://pubmed.ncbi.nlm.nih.gov/17675517/).
