# Neuroprotective Transcription Factors in Multiple Sclerosis

## Table of Contents
- [1. Project Overview](#1-project-overview)
- [2. Biological Background](#2-biological-background)
  - [2.1 Multiple Sclerosis and EAE](#21-multiple-sclerosis-and-eae)
  - [2.2 Neuronal Immune Signatures](#22-neuronal-immune-signatures)
  - [2.3 Transcription Factors as Therapeutic Targets](#23-transcription-factors-as-therapeutic-targets)
- [3. Reference Literature](#3-reference-literature)
  - [3.1 Biorxiv Paper: Mohammadnia et al., 2024](#31-biorxiv-paper-mohammadnia-et-al-2024)
  - [3.2 Gianfelice Thesis, 2025](#32-gianfelice-thesis-2025)
- [4. Project Goals](#4-project-goals)
- [5. Data Sources](#5-data-sources)
  - [5.1 Motor Neuron Dataset (GSE104897)](#51-motor-neuron-dataset-gse104897)
  - [5.2 Retinal Ganglion Cell Dataset (GSE247173)](#52-retinal-ganglion-cell-dataset-gse247173)
  - [5.3 Published Gene Lists](#53-published-gene-lists)
- [6. Computational Pipeline](#6-computational-pipeline)
  - [6.1 Pipeline Overview](#61-pipeline-overview)
  - [6.2 Step 1: Data Acquisition](#62-step-1-data-acquisition)
  - [6.3 Step 2: Data Wrangling (DESeq2)](#63-step-2-data-wrangling-deseq2)
  - [6.4 Step 3: Exploratory Data Analysis](#64-step-3-exploratory-data-analysis)
  - [6.5 Step 4: Functional Enrichment (ssGSEA)](#65-step-4-functional-enrichment-ssgsea)
  - [6.6 Step 5: Gene Intersection](#66-step-5-gene-intersection)
  - [6.7 Step 6: GO Enrichment and Term Reduction](#67-step-6-go-enrichment-and-term-reduction)
  - [6.8 Step 7: Gene Regulatory Network Inference (iRegulon)](#68-step-7-gene-regulatory-network-inference-iregulon)
  - [6.9 Step 8: Network Centrality Analysis (Cytoscape)](#69-step-8-network-centrality-analysis-cytoscape)
  - [6.10 Step 9: Subregulatory Network Generation](#610-step-9-subregulatory-network-generation)
  - [6.11 Step 10: Visualization and Reporting](#611-step-10-visualization-and-reporting)
- [7. Supplementary Analyses](#7-supplementary-analyses)
  - [7.1 Comparison with Published 97-Gene Signature](#71-comparison-with-published-97-gene-signature)
  - [7.2 Transcription Factor Clustering (Jaccard Similarity)](#72-transcription-factor-clustering-jaccard-similarity)
  - [7.3 Community Detection (Leiden Algorithm)](#73-community-detection-leiden-algorithm)
  - [7.4 Neurite Outgrowth Plate Analysis](#74-neurite-outgrowth-plate-analysis)
  - [7.5 Combinatorial GO Reduction Sweep](#75-combinatorial-go-reduction-sweep)
- [8. Key Results and Findings](#8-key-results-and-findings)
  - [8.1 Differentially Expressed Genes](#81-differentially-expressed-genes)
  - [8.2 Shared Immune Gene Signature](#82-shared-immune-gene-signature)
  - [8.3 Identified Transcription Factors](#83-identified-transcription-factors)
  - [8.4 Gene Regulatory Network Structure](#84-gene-regulatory-network-structure)
  - [8.5 GO Biological Process Categories](#85-go-biological-process-categories)
  - [8.6 TF Activity Across GO Categories](#86-tf-activity-across-go-categories)
  - [8.7 TF Functional Clusters](#87-tf-functional-clusters)
  - [8.8 Community Structure](#88-community-structure)
  - [8.9 Neurite Outgrowth Experimental Results](#89-neurite-outgrowth-experimental-results)
- [9. Directory Structure](#9-directory-structure)
- [10. Dependencies](#10-dependencies)
- [11. Running the Pipeline](#11-running-the-pipeline)
- [12. Current Status and Next Steps](#12-current-status-and-next-steps)

---

## 1. Project Overview

This project is a **bioinformatics and computational neuroscience** investigation into **neuroprotective transcription factors (TFs)** in the context of **Multiple Sclerosis (MS)**, using the **Experimental Autoimmune Encephalomyelitis (EAE)** mouse model. The work is conducted in the **Fournier Lab at the Montreal Neurological Institute, McGill University**.

The central discovery motivating this project is that neurons in MS and EAE are not merely passive victims of inflammation — they actively upregulate immune-related genes, including MHC class I molecules and interferon response genes. This neuronal immune signature is regulated by transcription factors, particularly from the **IRF** and **STAT** families, and inhibition of these TFs (notably **IRF5**) can protect neurons from inflammation-mediated cell death.

This computational pipeline re-analyzes publicly available RNA-sequencing data from two neuronal cell types — **motor neurons (MNs)** and **retinal ganglion cells (RGCs)** — to:
1. Identify shared upregulated genes between these cell types in EAE
2. Characterize the biological processes these genes are involved in
3. Infer the gene regulatory networks (GRNs) controlling these genes
4. Identify the key transcription factor regulators
5. Build subregulatory networks organized by biological process
6. Rank and cluster TFs by their regulatory roles

The experimental companion to this computational work (Gianfelice, 2025) tests transcription factor inhibitors (TFIs) for neuroprotective effects in vitro, providing biological validation of the computational predictions.

---

## 2. Biological Background

### 2.1 Multiple Sclerosis and EAE

**Multiple Sclerosis (MS)** is a chronic, immune-mediated disease of the central nervous system (CNS) characterized by:
- **Inflammation**: Infiltration of peripheral immune cells (T-cells, B-cells) across the blood-brain barrier
- **Demyelination**: Loss of myelin sheaths around axons, impairing signal conduction
- **Neurodegeneration**: Progressive neuronal and axonal damage leading to permanent disability

MS is broadly classified into:
- **Relapsing-Remitting MS (RRMS)**: ~85% of initial diagnoses, with alternating relapses and remission
- **Secondary Progressive MS (SPMS)**: Gradual worsening after initial RRMS phase
- **Primary Progressive MS (PPMS)**: ~10-15% of cases, steady decline from onset

Current disease-modifying therapies (DMTs) primarily target the **immune system** (e.g., B-cell depletion, natalizumab, fingolimod) but do not address the underlying **neuronal pathology** that drives progressive disability. This is the critical unmet need this project aims to address.

**Experimental Autoimmune Encephalomyelitis (EAE)** is the most widely used animal model for MS. It is induced by immunizing mice with myelin antigens (e.g., MOG35-55) and mimics key features of MS including:
- Immune cell infiltration into the CNS
- Demyelination
- Neuronal damage
- Paralysis (clinical scoring)

EAE allows researchers to study disease mechanisms in a controlled, reproducible system and test potential therapeutic interventions before moving to human trials.

### 2.2 Neuronal Immune Signatures

Historically, research on immune processes in MS has focused on **glial cells** (astrocytes, microglia, oligodendrocytes) and **infiltrating immune cells**. Neurons were considered passive targets of inflammation. However, recent evidence shows that neurons themselves actively participate in the inflammatory response:

- **MHC Class I expression**: Neurons can present antigens via MHC-I molecules when stimulated with interferon-gamma, potentially making them targets for T-cell-mediated attack
- **Interferon signaling**: Neurons express interferon receptors and can activate downstream STAT/IRF signaling cascades
- **Immune gene upregulation**: In both MS patient neurons and EAE mouse neurons, there is broad upregulation of immune-related genes

This project characterizes this **neuronal immune signature** in two different neuron types to find conserved patterns and identify the transcription factors driving this response.

### 2.3 Transcription Factors as Therapeutic Targets

Transcription factors (TFs) are proteins that bind to DNA regulatory regions to control gene expression. They are particularly attractive therapeutic targets because:
- A single TF can regulate dozens to hundreds of downstream genes
- Targeting one TF can modulate an entire pathological gene program
- Small molecule inhibitors of TFs can cross the blood-brain barrier

The key TF families identified in this project include:
- **IRF family** (IRF1, IRF2, IRF5, IRF7, IRF8, IRF9): Regulate interferon signaling and immune responses. **IRF5** is particularly notable as it overlaps with MS GWAS susceptibility loci.
- **STAT family** (STAT1, STAT2, STAT3, STAT6): Signal transducers activated by interferon receptor signaling
- **NF-κB family** (RelA/p65, NFKB1): Master regulators of inflammation
- **ETS family** (ETS1, ETS2, ELF4, FLI1): Regulate immune cell development and function
- **PU.1 (SPI1)**: A key myeloid and immune TF, also implicated in neuronal immune gene regulation

---

## 3. Reference Literature

### 3.1 Biorxiv Paper: Mohammadnia et al., 2024

**"Uncovering immune characteristics of neurons in the context of multiple sclerosis and experimental autoimmune encephalomyelitis"**
*Mohammadnia A, Drake SS, Zandee SEJ, Gosselin D, Antel JP, Prat A, Fournier AE*
bioRxiv preprint, August 2024

This is the foundational paper for this project. The key findings are:

1. **Shared immune gene expression in EAE neurons**: RNA-seq of MNs and RGCs from EAE mice identified numerous genes significantly upregulated in both cell types (p < 0.05, logFC > 0). GO analysis of these shared genes revealed enrichment for immune-related processes.

2. **ssGSEA confirms immune pathway enrichment**: Single-sample GSEA on both MN and RGC datasets showed extensive enrichment of immune pathways, including "cellular response to type II interferon", "type I interferon production", "response to type I interferon", and "antigen processing and presentation".

3. **ATAC-seq confirms chromatin accessibility**: ATAC-seq data from RGCs showed that immune pathway genes were in the top quartiles of chromatin accessibility, indicating these genes are in open, transcriptionally active chromatin states.

4. **97-gene common signature with human MS**: By comparing EAE neuron data with snRNA-seq from MS patient cortical neurons, the authors identified **97 genes** commonly upregulated across all datasets. This signature was enriched in neurons from chronic active MS lesions.

5. **Transcription factor prediction**: Using **iRegulon** on the top 48 accessible genes from the 97-gene signature, the authors predicted upstream TF regulators. The top regulators were from the **IRF** and **STAT** families. **Leiden clustering** of TF protein-protein interactions (via STRING DB) revealed a highly interconnected subnetwork of IRF and STAT family members.

6. **IRF5 as a key regulator**: Cross-referencing predicted TFs with MS GWAS susceptibility genes identified **IRF5** as the only overlapping transcription factor. IRF5 was predicted to regulate numerous members of the 48-gene signature.

7. **IRF5 inhibition is neuroprotective**: Treatment of primary mouse cortical neurons with the IRF5 inhibitor **YE6144** (1 µM) significantly promoted neuronal survival in inflammatory conditions (PBMC-conditioned media) and reduced enrichment of the 48-gene immune signature.

**Relevance to this project**: This paper provides the biological rationale and analytical framework. Our computational pipeline reproduces and extends this analysis with additional GO term reduction methods, TF clustering approaches, subregulatory network construction, and combinatorial parameter optimization.

### 3.2 Gianfelice Thesis, 2025

**"Investigating Transcription Factor Inhibitors for Neuroprotection in Multiple Sclerosis"**
*Christine Gianfelice*, Master's Thesis, McGill University, 2025

This thesis represents the **experimental validation** arm of the project, conducted in parallel with this computational work. Only **Aim 1** is relevant here:

**Aim 1: Investigate Transcription Factor Inhibitors for Neuroprotection**

Building on the biorxiv paper's identification of IRF5 as a key regulator, this aim tested four additional TFIs for neuroprotective effects:

| Inhibitor | Target TF | Key Finding |
|-----------|-----------|-------------|
| **YE6144** | IRF5 | Previously shown to protect neurons (Mohammadnia et al., 2024) |
| **DB2313** | PU.1 (SPI1) | **Significantly protected cortical neurons** from PBMC-CM-induced cell death (p < 0.05) |
| **Fludarabine** | STAT1 | Showed **dose-dependent trend** toward neuroprotection (10 nM – 1 µM), though PBMC-CM did not induce significant baseline death in this experiment |
| **SC75741** | NF-κB (RelA) | Showed **trend toward neuroprotection** (10 nM – 1 µM), similar pattern to Fludarabine |
| **Auranofin** | IRF3 | **Toxic to neurons** at all concentrations tested (even 1 pM), ruled out |

**Methodology**: Cryopreserved rat embryonic cortical neurons (E17-18) were cultured for 7 days, then subjected to 80% media swap with PBMC-conditioned media (from PMA/ionomycin-stimulated PBMCs) supplemented with TFIs or DMSO control. After 72 hours, neurons were fixed, stained (β-III tubulin + Hoechst), imaged, and quantified for cell survival.

**Key results from neurite outgrowth plate analysis** (`neurite_outgrowth_plate_analysis/`):
- Two-way ANOVA was performed on % β-III tubulin positive cells and total cell counts
- YE6144 (IRF5 inhibitor) was tested at 0.1, 1, and 10 µM in Control vs PBMC media
- Results are visualized as bar plots with t-test significance annotations

**Relevance to this project**: The experimental findings validate the computational predictions. DB2313 (PU.1 inhibitor) showing significant neuroprotection is particularly interesting because PU.1/SPI1 was identified by our pipeline as a key TF in the regulatory network.

---

## 4. Project Goals

The specific goals of this computational project are:

1. **Reproduce the published analysis** using publicly available data, with transparent, reproducible code
2. **Identify shared upregulated genes** between motor neurons and retinal ganglion cells in EAE
3. **Characterize biological processes** through GO enrichment analysis with semantic similarity reduction
4. **Infer gene regulatory networks** using iRegulon to predict TF → target gene relationships
5. **Rank transcription factors** by centrality metrics (in-degree, out-degree, betweenness, closeness)
6. **Build subregulatory networks** organized by GO biological process parent categories
7. **Assess TF specificity** to different biological processes through enrichment scoring
8. **Cluster TFs** by target gene overlap (Jaccard similarity) and community structure (Leiden algorithm)
9. **Compare results** with the published 97-gene signature from the Mohammadnia et al. paper
10. **Optimize parameters** across multiple GO reduction methods, similarity thresholds, and gene assignment modes

---

## 5. Data Sources

### 5.1 Motor Neuron Dataset (GSE104897)

- **Source**: NCBI Gene Expression Omnibus (GEO)
- **Organism**: Mouse (Mus musculus)
- **Cell type**: Motor neurons (translating ribosome affinity purification, TRAP-seq)
- **Condition**: EAE (MOG35-55 immunization) vs. Naive/Healthy controls
- **Samples**: 20 per condition (40 total)
- **Raw data**: HTSeq count files (`.txt.gz`) extracted from `GSE104897_RAW.tar`
- **Processing**: DESeq2 normalization and differential expression (EAE vs. Healthy)
- **Key parameters**: p-value < 0.05, log2FC > 0 for upregulated genes
- **Reference**: Schattling et al. (2019), Drake et al. (2023)

### 5.2 Retinal Ganglion Cell Dataset (GSE247173)

- **Source**: NCBI GEO
- **Organism**: Mouse (Mus musculus)
- **Cell type**: Retinal ganglion cells (bulk RNA-seq)
- **Condition**: EAE vs. Naive/Healthy controls
- **Samples**: Specific naive, onset, and peak timepoints selected from metadata
- **Metadata**: Sample information from `GSE247173_RNA_Metadata_EAE_RGC.xlsx`
- **Processing**: DESeq2 normalization and differential expression
- **ID mapping**: Gene symbols mapped to Ensembl IDs via biomaRt (with mirror fallbacks)
- **Reference**: Woo et al. (2024)

### 5.3 Published Gene Lists

- **97 common genes** (`List of 97 common genes.xlsx`): The gene signature from Mohammadnia et al. (2024) representing genes commonly upregulated across EAE MNs, EAE RGCs, and human MS patient neurons. Used for validation via Venn diagram overlap and correlation analysis.

---

## 6. Computational Pipeline

### 6.1 Pipeline Overview

The pipeline is orchestrated by `0.Pipeline_execution.R`, which runs all steps sequentially:

```
install_packages.R
        |
iregulon_reformat.R (if iRegulon results exist)
        |
   ┌────┴────┐
   v         v
Motor Neuron    RGC Analysis
Analysis        (Scripts 1-4)
(Scripts 1-4)       |
   │         └──────┘
   v         v
5.gene_intersection.R → 215 shared genes
        |
6.GO_analysis.R → GO enrichment + reduction
        |
7.cytoscape_output_analysis.R → Centrality analysis
        |
8.generate_subregulatory_networks.R → Per-GO subnetworks
        |
9.visualize_subnetworks.R → PDF report + bubble plots
```

### 6.2 Step 1: Data Acquisition

**Motor Neurons** (`Motor Neuron Analysis/1.data_aquisiton.R`):
- Downloads `GSE104897_RAW.tar` from GEO via `GEOquery`
- Extracts 40 HTSeq count files (20 EAE + 20 Control)
- Builds a count matrix with Ensembl gene IDs as row names
- Fetches phenotype metadata (disease state, tissue info)

**RGCs** (`RGC analysis/1.data_aquisition.R`):
- Downloads count data and metadata from GSE247173
- Reads sample metadata from Excel spreadsheet
- Selects specific timepoints (naive, onset, peak) for analysis
- Builds count matrix from selected samples

### 6.3 Step 2: Data Wrangling (DESeq2)

**Motor Neurons** (`Motor Neuron Analysis/2.data_wrangling.R`):
- Constructs DESeqDataSet from count matrix
- Runs full DESeq2 pipeline: estimation, testing, shrinkage
- Contrast: EAE vs. Healthy
- Maps Ensembl IDs to gene symbols using `EnsDb.Mmusculus.v79`
- Exports results: baseMean, log2FC, p-value, adjusted p-value, gene symbol

**RGCs** (`RGC analysis/2.data_wrangling.R`):
- Same DESeq2 workflow
- Additional ID mapping step: Ensembl IDs → gene symbols via `biomaRt` (with mirror fallbacks for bioconductor mirror issues)
- Handles duplicate mappings and missing IDs

### 6.4 Step 3: Exploratory Data Analysis

Both cell types undergo identical EDA (`3.data_analysis.R`):
- **PCA plots**: Visualize sample clustering by condition (EAE vs. Healthy)
- **Volcano plots**: Show distribution of log2FC vs. -log10(p-value), highlighting significant genes
- **Heatmaps**: Top upregulated genes (MN: log2FC > 10; RGC: log2FC > 2) with hierarchical clustering
- **Summary statistics**: Number of significant upregulated/downregulated genes

### 6.5 Step 4: Functional Enrichment (ssGSEA)

Both cell types undergo ssGSEA (`4.functional_enrichment_analysis.R`):
- Uses MSigDB **GO: Biological Process** gene sets
- Runs single-sample GSEA on normalized count matrices
- Identifies top 50 enriched pathways by enrichment score delta (EAE - Healthy)
- Generates clustered heatmaps of pathway enrichment across samples
- Key finding: Immune/interferon pathways are top enriched in both cell types

### 6.6 Step 5: Gene Intersection

`5.gene_intersection.R` performs:
- Extracts significantly upregulated genes from each cell type (p < 0.05, positive log2FC)
- Computes intersection: **215 shared upregulated genes**
- Generates **Venn diagram** showing overlap
- Exports:
  - `intersected_upregulated_genes.txt` (215 gene symbols)
  - `intersected_upregulated_ensembl_ids.txt` (Ensembl IDs)
  - `intersected_upregulated_IDs.txt` (formatted for g:Profiler)
  - `Expression_data.txt` (combined DESeq2 results for all genes, both cell types)

### 6.7 Step 6: GO Enrichment and Term Reduction

`6.GO_analysis.R` performs two stages:

**Stage 1: GO enrichment** (run externally via g:Profiler web tool):
- Input: 215 shared genes
- Output: `GO_results.csv` (687 enriched terms across GO:BP, GO:MF, GO:CC)
- Plots top 20 GO:BP terms as bar plots

**Stage 2: GO term reduction** (two methods compared):
1. **ReviGO** (via `rrvgo` R package): Uses semantic similarity to collapse redundant GO terms into representative "parent" categories
2. **simplifyEnrichment**: Uses hierarchical clustering of GO term overlap to group similar terms

Both methods produce:
- Reduced term sets with parent categories
- Gene-to-parent-category mappings
- Similarity matrices at multiple thresholds (0.70 – 0.95)

The three major parent categories that emerge:
- **cellular_process**
- **immune_system_process**
- **response_to_external_stimulus**

### 6.8 Step 7: Gene Regulatory Network Inference (iRegulon)

`iregulon_reformat.R` processes raw iRegulon output:
- **iRegulon** was run externally in Cytoscape on the 215 shared genes
- Raw results (`iregulon_results.csv`): 119 motif-target associations with NES scores
- The reformatting function:
  1. Parses the iRegulon TSV, extracting TF, target genes, motif ID, and NES
  2. Separates multi-gene/multi-TF entries into individual edges
  3. Filters TFs: must be significantly upregulated in at least one dataset AND not in opposite directions between MN and RGC
  4. Classifies each gene's specificity: MN_Only, RGC_Only, Both, or Opposite
  5. Computes an `EAE_Combined_Score` (average log2FC across cell types)
- Outputs:
  - `reformatted_edges.csv`: **3,739 TF → target gene edges**
  - `reformatted_nodes.csv`: **50 TF nodes** with expression data
  - `regulators.txt`: List of 50 TF names

### 6.9 Step 8: Network Centrality Analysis (Cytoscape)

`7.cytoscape_output_analysis.R` analyzes the Cytoscape-enriched network:
- Imports node table with STRING DB protein interaction data and tissue expression annotations
- Calculates centrality metrics:
  - **Out-degree**: Number of genes a TF regulates (regulatory influence)
  - **In-degree**: Number of TFs regulating a gene (degree of co-regulation)
- Visualizations:
  - Bar chart of TF regulatory impact (out-degree ranking)
  - Bar chart of top co-regulated targets (in-degree ranking)
  - Scatter plot: TF expression (MN log2FC vs. RGC log2FC), sized by out-degree

### 6.10 Step 9: Subregulatory Network Generation

`8.generate_subregulatory_networks.R` is the core network analysis script:

1. **Map genes to GO parent categories**: Each of the 215 genes is assigned to one or more parent GO terms from the reduced term set
2. **Optional unique mode**: Each gene assigned to only its most significant parent category (reduces overlap)
3. **Build master network**: From reformatted iRegulon edges
4. **For each GO parent category**:
   - Extract subnetwork (edges where both TF and target are in that category)
   - Calculate **igraph centrality metrics**: in-degree, out-degree, betweenness, closeness
   - Run **fgsea** enrichment for each TF's target set against the MN+RGC combined expression ranking
   - Compute NES (Normalized Enrichment Score) for each TF in each category
   - Export per-category edge lists, node metadata, and TF regulator lists
5. **Output**: `tf_target_enrichment.csv` (40 rows: TF × GO category enrichment scores)

### 6.11 Step 10: Visualization and Reporting

`9.visualize_subnetworks.R` generates:
- **Bubble plot (dotplot)**: TF activity across GO parent categories
  - X-axis: GO parent categories
  - Y-axis: Transcription factors
  - Size: % of subnetwork genes regulated by that TF (out-degree / total targets)
  - Color: NES enrichment score from fgsea
- **PDF report** (`subnetwork_report.pdf`):
  - Title page
  - Network legend
  - One page per GO category: Cytoscape network image + top TFs bar chart + top targets bar chart

---

## 7. Supplementary Analyses

### 7.1 Comparison with Published 97-Gene Signature

`comparing_to_old_analysis.R` validates our results against the published gene list:
- Loads the 97 genes from Mohammadnia et al.
- Computes Venn diagram overlap with our 215 shared genes
- Calculates overlap percentages in both directions
- Generates scatter plots comparing:
  - p-values between our analysis and the published supplementary data
  - log2 fold changes between the two analyses
- Identifies genes unique to each analysis

### 7.2 Transcription Factor Clustering (Jaccard Similarity)

Located in `tf_clustering/`, this analysis:
1. Computes a **Jaccard similarity matrix** between all 50 TFs based on overlap of their target gene sets
2. Performs **hierarchical clustering** (Ward D2 method)
3. Tests multiple cut heights (0.5, 1.0, 1.5) to define clusters
4. Labels each cluster via `enrichGO()` on the union of target genes in that cluster
5. Runs **fgsea** to score TF activity per cluster
6. Generates:
   - Dendrograms showing TF relationships
   - Dotplots showing TF cluster membership and enrichment scores
   - Cluster enrichment CSV files

**Key TF clusters identified** (5 functional clusters):
- **C1**: Myeloid leukocyte activation (red)
- **C2**: Antigen processing and presentation (yellow)
- **C3**: Leukocyte migration (green)
- **C4**: Regulation of immune effector process (blue)
- **C5**: Response to molecule of bacterial origin (purple)

### 7.3 Community Detection (Leiden Algorithm)

Located in `community_detection/`, this analysis:
1. Builds an **undirected graph** from TF-target edges
2. Runs **Leiden community detection** at multiple resolutions (0.5 – 2.0)
3. Identifies communities of densely connected genes
4. Labels communities via GO enrichment of member genes
5. Computes TF activity per community using fgsea
6. Generates dotplots showing TF specificity across communities

**Communities identified** include:
- Antigen processing and presentation
- Positive regulation of chemokine production
- Various immune signaling modules

### 7.4 Neurite Outgrowth Plate Analysis

Located in `neurite_outgrowth_plate_analysis/`, this analyzes experimental plate reader data:
- Input: `JP1.xlsx` (raw plate data from high-content screening)
- Plate layout: Control vs PBMC media × DMSO/YE6144 (0.1, 1, 10 µM)
- Processing: Parses well positions, collapses ROIs per well, maps treatments
- Statistics: **Two-way ANOVA** (Media × Drug interaction)
- Outputs: Bar plots of % β-III tubulin positive cells and total cell counts with t-test annotations

### 7.5 Combinatorial GO Reduction Sweep

`run_all.R` systematically tests **24 parameter combinations**:
- **2 gene modes**: All genes (genes can belong to multiple categories) vs. Unique genes (one category per gene)
- **2 GO reduction methods**: ReviGO (rrvgo) vs. simplifyEnrichment
- **6 similarity thresholds**: 0.70, 0.75, 0.80, 0.85, 0.90, 0.95

Each combination produces:
- Dotplot visualizations
- GO reduction results
- Per-category gene lists
- TF enrichment scores
- Summary comparison images

Results are organized in:
- `subsetted options/` (all-genes mode)
- `subsetted options_unique/` (unique-genes mode)

---

## 8. Key Results and Findings

### 8.1 Differentially Expressed Genes

| Cell Type | Upregulated | Downregulated | Total Significant |
|-----------|-------------|---------------|-------------------|
| Motor Neurons (EAE vs. Healthy) | ~2,500+ | ~2,000+ | ~4,500+ |
| RGCs (EAE vs. Healthy) | ~2,000+ | ~1,500+ | ~3,500+ |

(DESeq2: p < 0.05, log2FC > 0 for upregulated)

### 8.2 Shared Immune Gene Signature

- **215 genes** are significantly upregulated in both MNs and RGCs in EAE
- These include immune/inflammation-related genes: `S100a8`, `S100a9`, `Cd68`, `Trem2`, `Tmem119`, `Tyrobp`, `Il1b`, `Tnf`, `Ccl2`, `Cxcl10`
- And transcription factors: `Stat1`, `Runx1`, `Spi1`, `Atf3`, `Cebpa`, `Nfkbia`, `Nfkbiz`
- This shared signature represents a conserved neuronal immune response across different CNS neuron types

### 8.3 Identified Transcription Factors

**50 transcription factors** were predicted by iRegulon as regulators of the shared gene set:

| Family | TFs | Key Role |
|--------|-----|----------|
| **IRF** | Irf1, Irf2, Irf5, Irf7, Irf8, Irf9 | Interferon signaling, immune gene regulation |
| **STAT** | Stat1, Stat2, Stat3, Stat6 | JAK-STAT signaling downstream of cytokine receptors |
| **NF-κB** | Rela, Nfkb1, Bcl3 | Master inflammatory regulators |
| **ETS** | Ets1, Ets2, Elf4, Ehf, Fli1 | Immune cell development, neuronal function |
| **PU.1/SPI** | Spi1 | Myeloid/immune TF, also active in neurons |
| **RUNX** | Runx1, Runx2, Runx3 | Hematopoiesis, immune cell differentiation |
| **CEBP** | Cebpa, Cebpb, Cebpd | Myeloid differentiation, acute phase response |
| **NFI** | Nfia, Nfib, Nfix | Neurodevelopment, also immune functions |
| **Other** | Sox10, Ltf, Nr1h2/3, Nr3c1, Yy1, Tbp, Bptf, Smad3, Jun, Tfeb, Myf6, Mxd4, Usf1, Pura | Diverse regulatory roles |

**Top TFs by expression fold change** (from reformatted_nodes.csv):
- Irf7: log2FC = 10.08 (MN)
- Irf5: log2FC = 9.21 (MN)
- Irf8: log2FC = 8.78 (MN)
- Bcl3: log2FC = 7.32 (MN)
- Elf4: log2FC = 8.57 (MN)

**TFs significant in both cell types**: Cebpa, Nfia, Runx1, Elf4, Ltf

### 8.4 Gene Regulatory Network Structure

The inferred regulatory network contains:
- **50 transcription factor nodes** (regulators)
- **~207 target gene nodes** (regulated genes)
- **3,739 directed edges** (TF → target relationships)
- **119 unique motif-TF associations** from iRegulon
- Top motifs: PU.1/IRF4/GABPB1 (NES=7.91), SPIB/PU.1 (NES=7.52), NF-κB family (NES=7.48)

### 8.5 GO Biological Process Categories

**443 GO:BP terms** were identified from the enrichment analysis. After reduction, **3 major parent categories** emerged:

1. **Immune system process** (strongest enrichment, p = 10^-56)
   - Includes: immune response, innate immune response, antigen processing and presentation, interferon signaling
   - Most TFs enriched here: IRF family, STAT family

2. **Response to external stimulus** (p = 10^-51)
   - Includes: response to bacterium, inflammatory response, defense response
   - Key TFs: NF-κB family, C/EBP family

3. **Cellular process** (broad category)
   - Includes: cell death, stress response, cell-cell communication
   - Key TFs: YY1, NFATC3, TBP, RUNX2, SMAD3

### 8.6 TF Activity Across GO Categories

From `tf_target_enrichment.csv` (fgsea NES scores across categories):

**Immune system process TFs** (highest enrichment):
- IRF family members (Irf1, Irf5, Irf7, Irf8, Irf9)
- STAT family (Stat1, Stat2, Stat3)
- Pura (purine-rich element binding protein A)

**Cellular process TFs** (highest enrichment):
- Yy1, Nfatc3, Tbp, Runx2, Bptf, Runx3, Smad3

**Response to external stimulus TFs**:
- Mix of IRF/STAT and NF-κB family members

### 8.7 TF Functional Clusters

Jaccard similarity clustering of TFs by target gene overlap identified **5 functional clusters**:

| Cluster | GO Enrichment | Key TFs |
|---------|---------------|---------|
| C1 (red) | Myeloid leukocyte activation | Cebpa, Cebpb, Spi1 |
| C2 (yellow) | Antigen processing and presentation | Irf1, Irf8, Stat1 |
| C3 (green) | Leukocyte migration | Ets1, Ets2, Fli1 |
| C4 (blue) | Regulation of immune effector process | Irf5, Irf7, Stat2 |
| C5 (purple) | Response to molecule of bacterial origin | Nfkb1, Rela, Bcl3 |

### 8.8 Community Structure

Leiden community detection identified distinct functional modules in the TF-target network:
- **Antigen processing and presentation community**: Strongly enriched for MHC class I pathway genes
- **Chemokine regulation community**: Genes involved in immune cell recruitment
- Additional communities corresponding to interferon signaling, cell death regulation, and inflammatory response

### 8.9 Neurite Outgrowth Experimental Results

From the plate analysis (`neurite_outgrowth_plate_analysis/`):
- **Two-way ANOVA** tested the interaction between Media (Control vs PBMC) and Drug (DMSO vs YE6144 at 3 concentrations)
- **% β-III tubulin positive cells**: Measures neuronal survival (β-III tubulin is a neuronal marker)
- **Total DAPI cell counts**: Measures overall cell survival
- YE6144 (IRF5 inhibitor) was tested as a proof-of-concept for the computational predictions

---

## 9. Directory Structure

```
Neuroprotective-Transcription-factors/
|
|-- PROJECT_OVERVIEW.md                     # This file
|-- .gitignore
|
|-- ========== PIPELINE SCRIPTS ==========
|-- 0.Pipeline_execution.R                  # Master orchestrator
|-- install_packages.R                      # Dependency installation
|-- iregulon_reformat.R                     # iRegulon output → edge/node tables
|-- 5.gene_intersection.R                   # Shared gene identification + Venn diagram
|-- 6.GO_analysis.R                         # GO enrichment + ReviGO/simplifyEnrichment
|-- 7.cytoscape_output_analysis.R           # Cytoscape centrality analysis
|-- 8.generate_subregulatory_networks.R     # Per-GO-category subnetworks + fgsea
|-- 9.visualize_subnetworks.R               # PDF report + bubble plots
|-- run_all.R                               # Combinatorial parameter sweep
|-- a.R                                     # Earlier iRegulon reformatting version
|-- comparing_to_old_analysis.R             # Validation against published 97-gene set
|
|-- ========== UTILITIES ==========
|-- utils/
|   |-- helpers.R                           # build_ranked(), load_expression(), load_tf_network()
|
|-- ========== MOTOR NEURON ANALYSIS ==========
|-- Motor Neuron Analysis/
|   |-- 1.data_aquisiton.R                  # Download GSE104897, build count matrix
|   |-- 2.data_wrangling.R                  # DESeq2 differential expression
|   |-- 3.data_analysis.R                   # PCA, volcano plots, heatmaps
|   |-- 4.functional_enrichment_analysis.R  # ssGSEA with MSigDB GO:BP
|   |-- GSE104897/                          # Raw + normalized count files (40 samples)
|
|-- ========== RGC ANALYSIS ==========
|-- RGC analysis/
|   |-- 1.data_aquisition.R                 # Download GSE247173, read metadata
|   |-- 2.data_wrangling.R                  # DESeq2 + biomaRt ID mapping
|   |-- 3.data_analysis.R                   # PCA, volcano plots, heatmaps
|   |-- 4.functional_enrichment_analysis.R  # ssGSEA with MSigDB GO:BP
|   |-- GSE247173/                          # Sample metadata Excel file
|
|-- ========== DATA FILES ==========
|-- Expression_data.txt                     # Combined DESeq2 results (25,436 genes)
|-- GO_results.csv                          # g:Profiler GO enrichment (687 terms)
|-- all_BP_terms.txt                        # GO:BP terms ranked by -log10(adj. p-val)
|-- iregulon_results.csv                    # Raw iRegulon output (119 entries)
|-- reformatted_edges.csv                   # TF→Target edge list (3,739 edges)
|-- reformatted_nodes.csv                   # TF node attributes (50 nodes)
|-- metadata.txt                            # Full node metadata (257 genes)
|-- regulators.txt                          # 50 identified TF regulators
|-- intersected_upregulated_genes.txt       # 215 shared gene symbols
|-- intersected_upregulated_ensembl_ids.txt # 215 genes as Ensembl IDs
|-- intersected_upregulated_IDs.txt         # Formatted for g:Profiler
|-- tf_target_enrichment.csv                # TF enrichment per GO category (40 rows)
|-- List of 97 common genes.xlsx            # Published reference gene set
|
|-- ========== VISUALIZATIONS ==========
|-- GRN_legend.gif / .png                   # Network style legend
|-- dotplot_allgenes_s0.85.png              # Dotplot, all-genes mode, threshold 0.85
|-- dotplot_unique.png                      # Dotplot, unique-genes mode
|-- subnetwork_report.pdf                   # Multi-page subnetwork PDF report
|-- subnetwork_report_revigo.pdf            # Report using ReviGO method
|-- cytoscape_analysis_fixed.cys            # Cytoscape workspace
|-- Parent_Subnetworks_Workspace.cys        # Cytoscape subnetwork workspace
|
|-- ========== GO REDUCTION SWEEP ==========
|-- subsetted options/                      # All-genes mode
|   |-- rrvgo/                              # ReviGO reduction
|   |   |-- 0.70/ ... 0.95/                # Threshold-specific results
|   |-- simplifyEnrichment/                 # simplifyEnrichment reduction
|       |-- 0.70/ ... 0.95/
|-- subsetted options_unique/               # Unique-genes mode
|   |-- rrvgo/
|   |-- simplifyEnrichment/
|
|-- ========== TF CLUSTERING ==========
|-- tf_clustering/
|   |-- README.txt                          # Documentation
|   |-- run_clustering.R                    # Jaccard similarity + hierarchical clustering
|   |-- dendrograms_labeled.pdf             # Labeled dendrogram
|   |-- dendrogram_h*.png                   # Dendrograms at various cut heights
|   |-- dotplot_h*.png                      # Dotplots at various cut heights
|   |-- tf_cluster_enrichment*.csv          # Enrichment data per height
|
|-- ========== COMMUNITY DETECTION ==========
|-- community_detection/
|   |-- README.txt                          # Documentation
|   |-- prepare_data.R                      # Edge preparation from iRegulon
|   |-- run_community.R                     # Leiden clustering + GO labeling
|   |-- tf_community_dotplot_res*.png       # Dotplots at various resolutions
|   |-- tf_community_enrichment.csv         # TF enrichment across communities
|
|-- ========== SUBNETWORK DATA ==========
|-- parent_subnetworks/                     # Per-GO-category network files
|   |-- cellular_process/                   # edges, metadata, regulators
|   |-- immune_system_process/
|   |-- response_to_external_stimulus/
|-- parent_category_genes/                  # Gene lists per parent category
|
|-- ========== EXPERIMENTAL DATA ==========
|-- neurite_outgrowth_plate_analysis/
|   |-- JP1.xlsx                            # Raw plate reader data
|   |-- makegraph.R                         # Two-way ANOVA + bar plots
|   |-- two_way_anova_plot.png              # Output visualization
|
|-- ========== REFERENCE MATERIAL ==========
|-- Background info/
|   |-- nextsteps.txt                       # Research notes and to-do items
|   |-- Biorxiv paper.pdf                   # Mohammadnia et al., 2024
|   |-- GIANFELICE_Christine_Neuroscience (1).pdf  # Thesis, 2025
|   |-- Neurite outgrowth assay immunostaining protocol.docx
|
|-- to_delete/                              # Staging area for files to remove
```

---

## 10. Dependencies

### R Packages

**CRAN:**
- `dplyr`, `tidyr`, `stringr` – Data manipulation
- `ggplot2`, `ggpubr` – Visualization
- `VennDiagram` – Venn diagram generation
- `pheatmap` – Heatmaps
- `readxl` – Excel file reading
- `tibble` – Tidy data frames

**Bioconductor:**
- `DESeq2` – Differential gene expression analysis
- `GEOquery` – Downloading datasets from NCBI GEO
- `EnsDb.Mmusculus.v79` – Mouse Ensembl gene annotation
- `biomaRt` – Ensembl ID ↔ gene symbol mapping
- `clusterProfiler`, `org.Mm.eg.db` – GO enrichment analysis
- `fgsea` – Fast pre-ranked gene set enrichment analysis
- `rrvgo` – ReviGO GO term semantic similarity reduction
- `RCy3` – R interface to Cytoscape
- `igraph` – Network analysis and community detection
- `GSVA` – Gene set variation analysis (ssGSEA)

**External Tools:**
- **Cytoscape** (with iRegulon and CentiScape plugins) – Network visualization and TF prediction
- **g:Profiler** – GO enrichment (run via web tool, results imported as CSV)
- **STRING DB** – Protein-protein interaction data (integrated via Cytoscape)

---

## 11. Running the Pipeline

### Full Pipeline
```r
source("0.Pipeline_execution.R")
```
This runs all steps end-to-end. Requires iRegulon to have been run externally in Cytoscape first, with results saved as `iregulon_results.csv`.

### Combinatorial Sweep
```r
source("run_all.R")
```
Runs all 24 combinations of gene mode × GO reduction method × similarity threshold. Saves results to `subsetted options/` and `subsetted options_unique/`.

### Individual Steps
Each script can be sourced independently after the prerequisite data files exist. Key entry points:
- `iregulon_reformat.R` – After iRegulon results are available
- `5.gene_intersection.R` – After DESeq2 results for both cell types
- `6.GO_analysis.R` – After GO_results.csv is available from g:Profiler
- `tf_clustering/run_clustering.R` – After reformatted_edges.csv exists
- `community_detection/run_community.R` – After reformatted_edges.csv exists

---

## 12. Current Status and Next Steps

### Completed
- Full pipeline runs end-to-end
- 215 shared upregulated genes identified between MNs and RGCs
- 50 transcription factors predicted as regulators
- 3,739 TF → target gene edges in the regulatory network
- GO enrichment with 3 major parent categories identified
- Subregulatory networks built for each GO category
- TF clustering (Jaccard) and community detection (Leiden) completed
- Combinatorial parameter sweep across 24 combinations
- Comparison with published 97-gene signature
- Neurite outgrowth plate analysis with two-way ANOVA
- Experimental validation of TFIs (DB2313, Fludarabine, SC75741) in cortical neurons

### Next Steps (from research notes)
1. **Rank TFs and plot centrality analysis** using average degree (your degree / total connections)
2. **Refine GO analysis** to better characterize which biological processes are affected
3. **Build sub-GRCs** (sub gene regulatory networks) for interesting processes — either by extracting genes and making new networks, or using the current network and removing non-included genes
4. **Group bubble plots** by x-axis size or clustering
5. **Change bubble plot to % out-degree** to remove bias from networks with many edges
6. **Run enrichment analysis** (fgsea on GO gene sets) with color encoding for enrichment scores
7. **Make unique gene lists per pathway** — genes that are only in that set and not in others
8. **Explore ML approaches** (e.g., random forest) to distinguish EAE vs. control, noting that bulk RNA-seq has limited sample sizes
9. **Consider VAE approaches** for dimensionality reduction, with caution about overgeneralization
10. **Final output**: Organized categories/summary of findings

---

*Last updated: July 2026*
*Fournier Lab, Montreal Neurological Institute, McGill University*
