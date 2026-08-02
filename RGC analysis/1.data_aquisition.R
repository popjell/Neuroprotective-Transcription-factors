library(GEOquery)
library(tidyverse)
library(R.utils)
library(DESeq2)
library(readxl)

library(EnsDb.Mmusculus.v79)

gsenum_RGC = 'GSE247173'
gse_RGC <- getGEOSuppFiles(gsenum_RGC)


RGC_excel_data <- read_xlsx("GSE247173/GSE247173_RNA_Metadata_EAE_RGC.xlsx")

samples_to_keep <- c("Sample_naive_n1", "Sample_naive_n2", "Sample_naive_n5", "Sample_naive_n4", "Sample_onset_n1", "Sample_onset_n2", "Sample_onset_n3", "Sample_onset_n4","Sample_peak_n2")

#keep only first column and columns in samples_to_keep
RGC_excel_data <- RGC_excel_data %>%
 dplyr::select(1, all_of(samples_to_keep))

# use a brace-expression so `.` and rownames(.) are available
count_matrix_RGC <- RGC_excel_data |>
  column_to_rownames(var = colnames(RGC_excel_data)[1])
# remove htseq-count summary rows (those that start with "__")
count_matrix_RGC <- count_matrix_RGC[!stringr::str_starts(rownames(count_matrix_RGC), "__"), , drop = FALSE]

count_matrix_RGC <- as.matrix(count_matrix_RGC)
storage.mode(count_matrix_RGC) <- "integer"


#make coldata
coldata_RGC <- data.frame(
 sample_id = colnames(count_matrix_RGC),
 row.names = colnames(count_matrix_RGC)
) |> 
 mutate(
  condition = case_when(
   str_detect(sample_id, "naive") ~ "Healthy",
   str_detect(sample_id, "onset") ~ "EAE",
   str_detect(sample_id, "peak") ~ "EAE",
   TRUE              ~ NA_character_
  ),
  condition = factor(condition, levels = c("Healthy", "EAE"))
 )


