# Load the cached raw RNA counts and metadata.
library(Matrix)
library(Seurat)
library(ggplot2)

x <- readRDS("data/new_regions.rds")
counts <- x$counts
meta <- x$meta

# Check the count source before interpreting the plots or running pseudobulk.
meta$recovered_UMIs <- colSums(counts)
count_source <- if (all(meta$recovered_UMIs == meta$nCount_RNA)) "raw RNA" else
  if (all(meta$recovered_UMIs == meta$nCount_SCT)) "SCT corrected" else "unknown"
message("Count source: ", count_source)
if (count_source != "raw RNA") stop("Run 01_extract_raw_counts.R first")
write.csv(meta[, c("cellId", "sample_id", "condition", "cns_region",
                   "regional_subtype", "cell_subtypes_intg", "nFeature_RNA", "nCount_RNA",
                   "mito_percent", "recovered_UMIs")],
          "results/qc_cells.csv", row.names = FALSE)

# Plot QC for the new tissue samples.
p <- ggplot(meta, aes(condition, nFeature_RNA, fill = condition)) +
  geom_violin(scale = "width") + facet_wrap(~cns_region) +
  theme_minimal() + labs(x = NULL, y = "Genes per nucleus")
ggsave("plots/qc.png", p, width = 8, height = 4, dpi = 200)

# Show the full published atlas so regions and reused cortex samples are visible.
all_meta <- read.delim("data/meta.tsv")
coords <- read.delim(gzfile("data/UMAP.coords.tsv.gz"), header = FALSE,
                     col.names = c("cellId", "UMAP_1", "UMAP_2"))
all_meta$UMAP_1 <- coords$UMAP_1[match(all_meta$cellId, coords$cellId)]
all_meta$UMAP_2 <- coords$UMAP_2[match(all_meta$cellId, coords$cellId)]
p <- ggplot(all_meta, aes(UMAP_1, UMAP_2, color = cns_region)) +
  geom_point(size = 0.15, alpha = 0.65) + theme_void() +
  labs(title = "Trobisch full atlas (published UMAP)", color = "Region")
ggsave("plots/full_atlas_umap.png", p, width = 9, height = 6, dpi = 200)

# Cluster the new cerebellum and spinal cord cells with Seurat.
brain <- CreateSeuratObject(counts, meta.data = meta, min.cells = 3)
brain <- NormalizeData(brain, verbose = FALSE)
brain <- FindVariableFeatures(brain, nfeatures = 2000, verbose = FALSE)
brain <- ScaleData(brain, features = VariableFeatures(brain), verbose = FALSE)
brain <- RunPCA(brain, npcs = 20, verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:20, seed.use = 1, verbose = FALSE)
brain <- FindNeighbors(brain, dims = 1:20, verbose = FALSE)
brain <- FindClusters(brain, resolution = 0.4, random.seed = 1, verbose = FALSE)

cell_info <- data.frame(cellId = colnames(brain),
                        UMAP_1 = Embeddings(brain, "umap")[, 1],
                        UMAP_2 = Embeddings(brain, "umap")[, 2],
                        cluster = as.character(Idents(brain)),
                        subtype = brain$regional_subtype,
                        region = brain$cns_region,
                        condition = brain$condition,
                        donor = brain$sample_id)
write.csv(cell_info, "results/cell_metadata.csv", row.names = FALSE)

p <- DimPlot(brain, group.by = "cns_region", pt.size = 0.2) +
  ggtitle("New cerebellum and spinal cord samples")
ggsave("plots/umap_regions.png", p, width = 9, height = 6, dpi = 200)
p <- DimPlot(brain, pt.size = 0.2, label = TRUE, repel = TRUE) +
  ggtitle("Seurat clusters")
ggsave("plots/umap_clusters.png", p, width = 9, height = 6, dpi = 200)

# Plot the neuronal labels on the same UMAP.
neuron_labels <- c("GRAN-EN", "FRMPD4-IN", "KCNC2-IN", "PC",
                   "VENTRAL-NEU", "DORSAL-NEU")
cell_info$neuron_subtype <- ifelse(cell_info$subtype %in% neuron_labels,
                                   cell_info$subtype, "Other")
p <- ggplot() +
  geom_point(data = cell_info[cell_info$neuron_subtype == "Other", ],
             aes(UMAP_1, UMAP_2), color = "grey85", size = 0.2) +
  geom_point(data = cell_info[cell_info$neuron_subtype != "Other", ],
             aes(UMAP_1, UMAP_2, color = neuron_subtype), size = 0.2) +
  facet_wrap(~region) + theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Neuronal populations", color = NULL)
ggsave("plots/umap_neurons.png", p, width = 11, height = 6, dpi = 200)

# Check marker expression for the main neuronal populations.
markers <- intersect(c("GABRA6", "SLC17A7", "CALB1", "PCP2", "GAD1",
                       "SNAP25", "MBP", "GFAP"), rownames(brain))
p <- DotPlot(brain, features = markers, group.by = "regional_subtype") +
  RotatedAxis() + theme_minimal() + labs(x = NULL, y = NULL)
ggsave("plots/annotation_markers.png", p, width = 9, height = 9, dpi = 200)
