import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/celltypist-mpl")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/celltypist-numba")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

# Load the sparse counts without making a dense copy.
counts = mmread("data/celltypist_counts.mtx").tocsr().T
genes = pd.read_csv("data/celltypist_genes.txt", header=None)[0].astype(str).to_numpy()
cells = pd.read_csv("data/celltypist_cells.txt", header=None)[0].astype(str).to_numpy()
cluster_table = pd.read_csv("results/cell_metadata.csv").set_index("cell")
keep = np.isin(cells, cluster_table.index)
counts = counts[keep]
cells = cells[keep]
clusters = cluster_table.loc[cells, "cluster"]

# CellTypist expects log-normalized counts per 10,000 transcripts.
adata = ad.AnnData(counts)
adata.obs_names = cells
adata.var_names = genes
adata.var_names_make_unique()
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

# Annotate cells and use the full-atlas Seurat clusters for majority voting.
prediction = celltypist.annotate(
    adata,
    model="Adult_Human_PrefrontalCortex.pkl",
    majority_voting=True,
    over_clustering=clusters.astype(str).to_numpy(),
)
labels = prediction.predicted_labels
labels.insert(0, "cell", labels.index)
labels.to_csv("results/celltypist_annotations.csv", index=False)
