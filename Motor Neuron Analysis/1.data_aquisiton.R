#Aquiring data from GEO

#packages
library(GEOquery)
library(tidyverse)
library(R.utils)
library(DESeq2)

#aquire data for Motor Neurons in EAE
#region
gsenum_MN = 'GSE104897'
gse_MN <- getGEOSuppFiles(gsenum_MN)

tar_file <- "GSE104897/GSE104897_RAW.tar"
target_dir <- "GSE104897/extracted_counts"

dir.create(target_dir, showWarnings = FALSE)
untar(tar_file, exdir = target_dir)

#aqurie raw read files
raw_files_MN <- list.files(target_dir, pattern = "raw.*\\.txt\\.gz$", full.names = TRUE)

#read into a count matrix
read_counts <- function(f) {
  sample_name <- basename(f) |> str_remove("_raw.*\\.txt\\.gz")
  
  df <- read.table(gzfile(f), header = FALSE, sep = "\t",
                   col.names = c("gene_id", "count")) |>
    dplyr::filter(!str_starts(gene_id, "__")) |>      # remove HTSeq summary rows
    column_to_rownames("gene_id")              
  
  colnames(df) <- sample_name
  df
}

# Rebuild the matrix
count_matrix_MN <- lapply(raw_files_MN, read_counts) |>
  bind_cols() |>
  mutate(across(everything(), as.integer)) |>
  as.matrix()

head(count_matrix_MN)

# Fetch the metadata
gse_meta_MN <- getGEO(gsenum_MN, GSEMatrix = TRUE)[[1]]
pheno_MN <- pData(gse_meta_MN)

# Clean up the characteristics columns
pheno_MN$cell_type <- pheno_MN$`cell type:ch1` |>
  str_replace("All spinal cord cells", "SC") |>
  str_replace("Motor neurons", "MN")

pheno_MN$condition <- pheno_MN$`disease state:ch1` |>
  str_replace("Healthy, not immunized", "Healthy") |>
  str_replace("EAE day 15", "EAE")

# Build a replicate number per group
pheno_MN <- pheno_MN |>
  group_by(cell_type, condition) |>
  mutate(rep = row_number()) |>
  ungroup()

# Construct sample names
pheno_MN$sample_name <- paste(pheno_MN$cell_type, pheno_MN$condition, paste("rep", pheno_MN$rep, sep = ""), sep = "-")

name_map_MN <- setNames(pheno_MN$sample_name, pheno_MN$geo_accession)
colnames(count_matrix_MN) <- name_map_MN[colnames(count_matrix_MN)]

#design matrix for DESeq2
coldata_MN <- data.frame(
  cell_type = factor(pheno_MN$cell_type),
  condition = factor(pheno_MN$condition, levels = c("Healthy", "EAE")),
  row.names = name_map_MN[pheno_MN$geo_accession]
)

#verify colum/row matching
all(rownames(coldata_MN) == colnames(count_matrix_MN))

#split data by cell type
count_matrix_MN <- count_matrix_MN[, coldata_MN$cell_type == "MN"]
coldata_MN <- coldata_MN[coldata_MN$cell_type == "MN", ]


#endregion