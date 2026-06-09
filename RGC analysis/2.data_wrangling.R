library(DESeq2)
library(EnsDb.Mmusculus.v79)
library(tidyverse)

#make DDS object
dds_RGC <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                 colData = coldata_filtered,
                 design = ~ condition)
#set condition factors for DESeq2
dds_RGC$condition <- factor(dds_RGC$condition, levels = c("Healthy", "EAE"))

#Run DESeq2
dds_RGC_analysis <- DESeq(dds_RGC)
RGC_results_raw <- results(dds_RGC_analysis, contrast = c("condition", "EAE", "Healthy"))
RGC_results.df <- as.data.frame(RGC_results_raw)
RGC_results.df$ensembl_id <- rownames(RGC_results.df)

#map ids
RGC_results.df$symbol <- mapIds(
  EnsDb.Mmusculus.v79,
  keys = RGC_results.df$ensembl_id,
  column = "GENENAME",
  keytype = "GENEID",
  multiVals = "first"
)





