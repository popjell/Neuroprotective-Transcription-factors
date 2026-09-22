# QC, cluster, and plot every Jakel nucleus.
library(Matrix)
library(Seurat)
library(ggplot2)

x <- readRDS("data/all_counts.rds")
counts <- x$counts
meta <- x$meta
rm(x)

meta$nCount_RNA <- colSums(counts)
meta$nFeature_RNA <- colSums(counts > 0)

# The GEO matrix is already filtered; remove only cells with fewer than 200 genes.
keep <- meta$nFeature_RNA >= 200
counts <- counts[, keep]
meta <- meta[keep, ]

brain <- CreateSeuratObject(counts, meta.data = meta, min.cells = 3)
brain <- NormalizeData(brain, verbose = FALSE)
brain <- FindVariableFeatures(brain, nfeatures = 2000, verbose = FALSE)
brain <- ScaleData(brain, features = VariableFeatures(brain), verbose = FALSE)
brain <- RunPCA(brain, npcs = 30, verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:20, seed.use = 1)
brain <- FindNeighbors(brain, dims = 1:20, verbose = FALSE)
brain <- FindClusters(brain, resolution = 0.4, random.seed = 1, verbose = FALSE)

umap <- data.frame(
  cell = colnames(brain),
  UMAP_1 = Embeddings(brain, "umap")[, 1],
  UMAP_2 = Embeddings(brain, "umap")[, 2],
  cluster = as.character(Idents(brain)),
  celltype = brain$Celltypes,
  patient = brain$Sample,
  condition = brain$Condition
)
write.csv(umap, "results/cell_metadata.csv", row.names = FALSE)

p <- DimPlot(brain, group.by = "Celltypes", label = TRUE, repel = TRUE,
             pt.size = 0.2) + NoLegend() + ggtitle("Jakel: all cell types")
ggsave("plots/umap.png", p, width = 10, height = 7, dpi = 200)

p <- DimPlot(brain, label = TRUE, repel = TRUE, pt.size = 0.2) +
  ggtitle("Jakel: Seurat clusters")
ggsave("plots/umap_clusters.png", p, width = 9, height = 7, dpi = 200)

p <- VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA"),
             group.by = "Condition", pt.size = 0, ncol = 2)
ggsave("plots/qc.png", p, width = 9, height = 4, dpi = 200)

# Recluster the neurons instead of using the authors' Neuron1-Neuron5 groups.
neuron_cells <- colnames(brain)[grepl("^Neuron", brain$Celltypes) &
                                 brain$nFeature_RNA >= 500 & brain$nCount_RNA >= 1000]
neurons <- subset(brain, cells = neuron_cells)
rm(brain)
neurons <- NormalizeData(neurons, verbose = FALSE)
neurons <- FindVariableFeatures(neurons, nfeatures = 2000, verbose = FALSE)
neurons <- ScaleData(neurons, features = VariableFeatures(neurons), verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 30, verbose = FALSE)
neurons <- RunUMAP(neurons, dims = 1:20, seed.use = 1)
neurons <- FindNeighbors(neurons, dims = 1:20, verbose = FALSE)
neurons <- FindClusters(neurons, resolution = 0.4, random.seed = 1, verbose = FALSE)

neuron_meta <- data.frame(
  cell = colnames(neurons),
  UMAP_1 = Embeddings(neurons, "umap")[, 1],
  UMAP_2 = Embeddings(neurons, "umap")[, 2],
  neuron_cluster = paste0("N", Idents(neurons)),
  patient = neurons$Sample,
  condition = neurons$Condition
)
write.csv(neuron_meta, "results/neuron_clusters.csv", row.names = FALSE)

p <- DimPlot(neurons, label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle("Jakel: reclustered neurons")
ggsave("plots/umap_neurons.png", p, width = 8, height = 6, dpi = 200)
