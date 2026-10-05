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

# Score simple marker sets and assign each cluster its best cell type.
markers <- list(
  Neuron = c("RBFOX3", "SNAP25", "SYT1", "TUBB3"),
  Astrocyte = c("AQP4", "GFAP", "SLC1A3", "ALDH1L1"),
  Oligodendrocyte = c("MOBP", "MBP", "PLP1", "MAG"),
  OPC = c("PDGFRA", "CSPG4", "VCAN"),
  Microglia = c("CX3CR1", "P2RY12", "TMEM119", "AIF1"),
  Endothelial = c("CLDN5", "EMCN", "KDR", "FLT1"),
  Pericyte = c("PDGFRB", "RGS5", "CSPG4", "MCAM")
)
markers <- lapply(markers, intersect, rownames(brain))
markers <- markers[lengths(markers) > 0]
brain <- AddModuleScore(brain, features = markers, name = "type", search = FALSE)
score_cols <- paste0("type", seq_along(markers))
cluster_scores <- aggregate(brain@meta.data[, score_cols],
                            list(cluster = Idents(brain)), mean)
cluster_label <- apply(cluster_scores[, -1, drop = FALSE], 1,
                       function(x) names(markers)[which.max(x)])
names(cluster_label) <- as.character(cluster_scores$cluster)
brain$self_annotation <- unname(cluster_label[as.character(Idents(brain))])

umap <- data.frame(
  cell = colnames(brain),
  UMAP_1 = Embeddings(brain, "umap")[, 1],
  UMAP_2 = Embeddings(brain, "umap")[, 2],
  cluster = as.character(Idents(brain)),
  self_annotation = brain$self_annotation,
  patient = brain$Sample,
  condition = brain$Condition
)
write.csv(umap, "results/cell_metadata.csv", row.names = FALSE)
write.csv(umap[, c("cell", "cluster", "patient", "condition")],
          "results/neuron_clusters.csv", row.names = FALSE)

p <- DimPlot(brain, label = TRUE, repel = TRUE, pt.size = 0.2) +
  ggtitle("Jakel: Seurat clusters")
ggsave("plots/umap_clusters.png", p, width = 9, height = 7, dpi = 200)

p <- DimPlot(brain, group.by = "self_annotation", label = TRUE, repel = TRUE,
             pt.size = 0.2) + ggtitle("Jakel: marker-based annotation")
ggsave("plots/umap_self_annotation.png", p, width = 10, height = 7, dpi = 200)

p <- VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA"),
             group.by = "Condition", pt.size = 0, ncol = 2)
ggsave("plots/qc.png", p, width = 9, height = 4, dpi = 200)
