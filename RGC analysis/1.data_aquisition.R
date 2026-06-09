library(GEOquery)
library(tidyverse)
library(R.utils)
library(DESeq2)
library(readxl)

library(EnsDb.Mmusculus.v79)

gsenum_RGC = 'GSE247173'
gse_RGC <- getGEOSuppFiles(gsenum_RGC)


RGC_excel_data <- read_xlsx("GSE247173/GSE247173_RNA_Metadata_EAE_RGC.xlsx")


count_matrix_RGC <- RGC_excel_data |> 
 column_to_rownames(var = colnames(RGC_excel_data)[1]) |> 
 as.matrix()
storage.mode(count_matrix_RGC) <- "integer"



#map to ensemble ids
edb <- EnsDb.Mmusculus.v79

current_ids <- rownames(count_matrix_RGC)

mapping_table <- genes(edb, 
                       filter = GeneNameFilter(current_ids), 
                       columns = c("gene_id", "gene_name"), 
                       return.type = "data.frame")

mapping_table <- mapping_table[!duplicated(mapping_table$gene_name), ]
valid_genes <- intersect(current_ids, mapping_table$gene_name)

count_matrix_RGC <- count_matrix_RGC[valid_genes, ]
mapping_ordered <- mapping_table[match(valid_genes, mapping_table$gene_name), ]

rownames(count_matrix_RGC) <- mapping_ordered$gene_id


#make coldata
coldata_RGC <- data.frame(
 sample_id = colnames(count_matrix_RGC),
 row.names = colnames(count_matrix_RGC)
) |> 
 mutate(
  condition = case_when(
   str_detect(sample_id, "naive") ~ "Healthy",
   str_detect(sample_id, "pre")  ~ "EAE",
   str_detect(sample_id, "onset") ~ "EAE",
   str_detect(sample_id, "peak") ~ "EAE",
   TRUE              ~ NA_character_
  ),
  condition = factor(condition, levels = c("Healthy", "EAE"))
 )



outliers_to_keep <- c("Sample_naive_n1", "Sample_naive_n2", "Sample_naive_n5", "Sample_naive_n4", "Sample_onset_n1", "Sample_onset_n2", "Sample_onset_n3", "Sample_onset_n4","Sample_peak_n2")
count_matrix_filtered <- count_matrix_RGC[, colnames(count_matrix_RGC) %in% outliers_to_keep]
coldata_filtered <- coldata_RGC[rownames(coldata_RGC) %in% outliers_to_keep, ]

