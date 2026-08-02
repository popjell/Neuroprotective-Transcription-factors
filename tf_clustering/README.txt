TF JACCARD CLUSTERING
======================

Purpose:
  Group TFs by how similar their target gene sets are (Jaccard index).
  TFs that regulate the same genes cluster together. Each cluster is
  then labeled by the biological process its target genes are enriched in.

Outputs:
  - dendrogram.pdf: hierarchical clustering dendrogram of TFs.
  - tf_cluster_dotplot.png: dot plot of TF similarity to each cluster.
    Dot size = Jaccard similarity of TF targets to the cluster centroid.
    Dot color = enrichment score (NES, avg of MN+RGC fgsea).
  - tf_cluster_enrichment.csv: the underlying data for the dot plot.

How it works:
  1. For each pair of TFs, compute Jaccard(targets(TF1), targets(TF2))
  2. Hierarchical clustering (ward.D2), cut at k=5
  3. Find the centroid TF of each cluster (most representative member)
  4. Label each cluster via enrichGO() on its target genes
  5. Score each TF against each cluster via fgsea on MN+RGC expression
  6. Compute Jaccard of each TF to each cluster's centroid for dot sizing

Files:
  - run_clustering.R         : the analysis script
  - reformatted_edges.csv    : TF -> target gene edges
  - Expression_data.txt      : expression data for fgsea

Run: Rscript run_clustering.R
