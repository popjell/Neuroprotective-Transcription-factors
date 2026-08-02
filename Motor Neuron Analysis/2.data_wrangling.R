library(DESeq2)
library(EnsDb.Mmusculus.v79)
library(tidyverse)

#Create a DESeq2 object
dds_MN <- DESeqDataSetFromMatrix(countData = count_matrix_MN,
                                  colData = coldata_MN,
                                  design = ~ condition)

#set condition factors for DESeq2
dds_MN$condition <- factor(dds_MN$condition, levels = c("Healthy", "EAE"))

#Run DESeq2
dds_MN_analysis <- DESeq(dds_MN)
MN_results_raw <- results(dds_MN_analysis, contrast = c("condition", "EAE", "Healthy"))
MN_results.df <- as.data.frame(MN_results_raw)
MN_results.df$ensembl_id <- rownames(MN_results.df)

#map ids
MN_results.df$symbol <- mapIds(
  EnsDb.Mmusculus.v79,
  keys = MN_results.df$ensembl_id,
  column = "GENENAME",
  keytype = "GENEID",
  multiVals = "first"
)





