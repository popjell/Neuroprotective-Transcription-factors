# QC, cluster, and plot the Schirmer neurons.
library(Matrix)
library(Seurat)
library(ggplot2)

x <- readRDS("data/neuron_counts.rds")
counts <- x$counts
meta <- x$meta
rm(x)
gc()

# Keep cells with at least 500 genes and 1,000 UMIs.
keep <- meta$genes >= 500 & meta$UMIs >= 1000
counts <- counts[, keep]
meta <- meta[keep, ]

neurons <- CreateSeuratObject(counts, meta.data = meta, min.cells = 3)
neurons <- NormalizeData(neurons, verbose = FALSE)
neurons <- FindVariableFeatures(neurons, nfeatures = 2000, verbose = FALSE)
neurons <- ScaleData(neurons, features = VariableFeatures(neurons), verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 30, verbose = FALSE)
neurons <- RunUMAP(neurons, dims = 1:20, seed.use = 1)
neurons <- FindNeighbors(neurons, dims = 1:20, verbose = FALSE)
neurons <- FindClusters(neurons, resolution = 0.4, random.seed = 1, verbose = FALSE)

# These are the same biological subtype groups used for pseudobulk.
neurons$pseudobulk_subtype <- as.character(neurons$cell_type)
neurons$pseudobulk_subtype[grepl("^EN-L2-3", neurons$pseudobulk_subtype)] <- "EN-L2-3"
neurons$pseudobulk_subtype[neurons$pseudobulk_subtype == "IN-PV"] <- "IN-PVALB"

umap <- data.frame(
  cell = colnames(neurons),
  UMAP_1 = Embeddings(neurons, "umap")[, 1],
  UMAP_2 = Embeddings(neurons, "umap")[, 2],
  cluster = as.character(Idents(neurons)),
  subtype = neurons$cell_type,
  patient = neurons$patient,
  condition = neurons$diagnosis
)
write.csv(umap, "results/cell_metadata.csv", row.names = FALSE)

p <- DimPlot(neurons, group.by = "pseudobulk_subtype", label = TRUE, repel = TRUE) +
  ggtitle("Schirmer neuronal subtypes used for pseudobulk")
ggsave("plots/umap.png", p, width = 9, height = 6, dpi = 200)

p <- DimPlot(neurons, label = TRUE, repel = TRUE) +
  ggtitle("Schirmer Seurat clusters")
ggsave("plots/umap_clusters.png", p, width = 9, height = 6, dpi = 200)

p <- VlnPlot(neurons, features = c("nFeature_RNA", "nCount_RNA"),
             group.by = "diagnosis", pt.size = 0, ncol = 2)
ggsave("plots/qc.png", p, width = 9, height = 4, dpi = 200)
