library(DESeq2)
library(biomaRt)
library(tidyverse)

#make DDS object
dds_RGC <- DESeqDataSetFromMatrix(countData = count_matrix_RGC,
                 colData = coldata_RGC,
                 design = ~ condition)
#set condition factors for DESeq2
dds_RGC$condition <- factor(dds_RGC$condition, levels = c("Healthy", "EAE"))

# 1. Estimate size factors (required for depth normalization)
dds_RGC <- estimateSizeFactors(dds_RGC)
normalized_count <- counts(dds_RGC, normalized = TRUE)


#Run DESeq2
dds_RGC_analysis <- DESeq(dds_RGC)
RGC_results_raw <- results(dds_RGC_analysis, contrast = c("condition", "EAE", "Healthy"), )
RGC_results.df <- as.data.frame(RGC_results_raw)
RGC_results.df$symbol <- rownames(RGC_results.df)


mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_symbols <- rownames(RGC_results.df)

mapping_table <- getBM(
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  filters = "mgi_symbol",
  values = gene_symbols,
  mart = mart
)

#join mapping and name "ensemble_gene_id" row as "ensembl_id"
RGC_results.df <- RGC_results.df %>% 
  left_join(mapping_table, by = c("symbol" = "mgi_symbol")) %>% 
  dplyr::rename(ensembl_id = ensembl_gene_id)
