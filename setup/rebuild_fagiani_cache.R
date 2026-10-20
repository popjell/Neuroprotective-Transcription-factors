atomic_save <- function(value, path) {
  temporary <- paste0(path, '.part')
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop('Could not finalize cache: ', path)
}
# Recover counts only if the deposited log-normalization reverses exactly.
# The source is normalized expression, not an independently deposited count matrix.
library(Matrix)
base <- 'scRNAseq/Fagiani_2025/data'
if (!file.exists(file.path(base,'neuron_counts.rds')) || '--force' %in% commandArgs(TRUE)) {
  message('Reading the deposited normalized matrix...')
  x <- readRDS(file.path(base,'GSE277430_normCount_GEOsubmission.rds'))
  message('Loaded matrix: ',nrow(x),' genes, ',ncol(x),' cells')
  meta <- read.delim(gzfile(file.path(base,'GSE277430_metadata_GEOsubmission.tsv.gz')))
  stopifnot(!anyDuplicated(meta$barcodes), all(colnames(x) %in% meta$barcodes))
  meta <- meta[match(colnames(x),meta$barcodes),]
  keep <- grepl('^NEU_',meta$expertAnno.l2.old) & meta$nFeature_RNA >= 500 &
          meta$nCount_RNA >= 1000 & meta$percent.mt <= 5
  prior_cells <- read.csv('scRNAseq/Fagiani_2025/results/cell_metadata.csv')$barcodes
  stopifnot(all(prior_cells %in% meta$barcodes))
  # Preserve the prior selection across the authors' annotation versions.
  keep <- meta$barcodes %in% prior_cells
  counts <- x[,keep,drop=FALSE]
  meta <- meta[keep,]
  rownames(meta) <- meta$barcodes
  rm(x); gc()
  message('Reversing normalization for ',ncol(counts),' prior-analysis cells...')
  values <- expm1(counts@x)
  counts@x <- numeric(0) # Detach the slot while transforming the unshared vector.
  gc()
  max_residual <- 0
  for (j in seq_len(ncol(counts))) {
    start <- counts@p[j]+1L; end <- counts@p[j+1L]
    if (end<start) next
    index <- start:end
    recovered <- values[index]*meta$nCount_RNA[j]/10000
    max_residual <- max(max_residual, max(abs(recovered-round(recovered))))
    if (max_residual>1e-5) stop('Normalization cannot be reversed exactly; raw counts required')
    values[index] <- round(recovered)
  }
  counts@x <- values
  rm(values); gc()
  counts <- drop0(counts)
  stopifnot(all(colSums(counts)==meta$nCount_RNA),
            all(colSums(counts>0)==meta$nFeature_RNA),
            identical(colnames(counts),meta$barcodes))
  atomic_save(counts,file.path(base,'neuron_counts.rds'))
  atomic_save(meta,file.path(base,'neuron_metadata.rds'))
  writeLines(c('Source: GSE277430_normCount_GEOsubmission.rds.gz',
     'Reconstruction: round(expm1(normalized)*nCount_RNA/10000)',
     paste('Maximum integer residual:',format(max_residual,digits=16)),
     paste('Neurons:',ncol(counts)),
     'Every reconstructed total and detected-gene count equals published metadata.'),
     file.path(base,'count_reconstruction_provenance.txt'))
  message('Fagiani: ',ncol(counts),' neurons; exact integer reconstruction and all published totals verified')
}
