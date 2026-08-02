COMMUNITY DETECTION (Leiden clustering)
========================================

Purpose:
  Cluster the TF-target network into gene communities using Leiden
  community detection on the undirected graph. Each community is a
  group of genes (TFs + targets) that are more connected to each
  other than to the rest of the network.

Outputs:
  - tf_community_dotplot.png: dot plot of TF activity across communities.
    Dot size = % of each TF's targets in that community.
    Dot color = enrichment score (NES, avg of MN+RGC fgsea).
  - tf_community_enrichment.csv: the underlying data for the dot plot.

How it works:
  1. Build undirected graph from TF-target edges
  2. Leiden clustering at user-chosen resolution (sweep shown first)
  3. Label each community via enrichGO() on its target genes
  4. Score each TF against each community via fgsea on MN+RGC expression
  5. Compute specificity (% of TF's total targets per community)

Files:
  - run_community.R          : the analysis script
  - prepare_data.R           : regenerates edges from iRegulon (run first if needed)
  - reformatted_edges.csv    : TF -> target gene edges
  - Expression_data.txt      : expression data for fgsea

Run: Rscript run_community.R
