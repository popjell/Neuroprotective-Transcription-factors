# QC, cluster, and plot the Fagiani neurons.
library(Matrix)
library(Seurat)
library(ggplot2)

counts <- readRDS("data/neuron_counts.rds")
meta <- readRDS("data/neuron_metadata.rds")
rownames(meta) <- meta$barcodes

# Keep cells with at least 500 genes, 1,000 UMIs, and at most 5% mitochondrial RNA.
keep <- meta$nFeature_RNA >= 500 & meta$nCount_RNA >= 1000 & meta$percent.mt <= 5
meta <- meta[keep, ]
counts <- counts[, rownames(meta)]

# Make one full neuronal UMAP. Do not split excitatory and inhibitory cells.
meta$subtype <- meta$expertAnno.l2.old
meta$cluster <- meta$subtype
write.csv(meta, "results/cell_metadata.csv", row.names = TRUE)

neurons <- CreateSeuratObject(counts[, rownames(meta)], meta.data = meta,
                              min.cells = 3)
neurons <- NormalizeData(neurons, verbose = FALSE)
neurons <- FindVariableFeatures(neurons, nfeatures = 2000, verbose = FALSE)
neurons <- ScaleData(neurons, features = VariableFeatures(neurons), verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 20, verbose = FALSE)
neurons <- RunUMAP(neurons, dims = 1:20, seed.use = 1, verbose = FALSE)

p <- DimPlot(neurons, group.by = "subtype", label = TRUE, repel = TRUE,
             pt.size = 0.1) + ggtitle("Fagiani: full neuronal UMAP")
ggsave("plots/umap.png", p, width = 9, height = 7, dpi = 200)

qc <- data.frame(meta)
p <- ggplot(qc, aes(disease, nFeature_RNA, fill = disease)) +
  geom_violin(scale = "width") + theme_minimal() + guides(fill = "none")
ggsave("plots/qc.png", p, width = 6, height = 4, dpi = 200)
