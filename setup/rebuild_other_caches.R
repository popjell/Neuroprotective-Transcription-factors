atomic_save <- function(value, path) {
  temporary <- paste0(path, '.part')
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop('Could not finalize cache: ', path)
}
# Run from the repository root after restore_other_data.py.
library(Matrix)
base <- 'scRNAseq/Schirmer_2019/data'
if (!file.exists(file.path(base, 'neuron_counts.rds'))) {
  rawdir <- file.path(base, 'raw')
  dir.create(rawdir, showWarnings=FALSE)
  unzip(file.path(base, 'rawMatrix.zip'), exdir=rawdir)
  genes <- read.delim(file.path(rawdir, 'genes.tsv'), header=FALSE,
                      col.names=c('ensembl_gene_id','gene_symbol'))
  cells <- readLines(file.path(rawdir,'barcodes.tsv'))
  meta <- read.delim(file.path(base,'meta.tsv'), check.names=FALSE)
  stopifnot(!anyDuplicated(meta$cell), all(cells %in% meta$cell))
  meta <- meta[match(cells,meta$cell),]
  raw <- as(readMM(file.path(rawdir,'matrix.mtx')),'CsparseMatrix')
  stopifnot(ncol(raw)==length(cells), nrow(raw)==nrow(genes), all(raw@x==round(raw@x)))
  rownames(raw) <- genes$ensembl_gene_id
  colnames(raw) <- cells
  keep <- grepl('^EN-|^IN-', meta$cell_type)
  counts <- raw[,keep,drop=FALSE]
  meta <- meta[keep,]
  meta$patient <- meta$sample
  rownames(meta) <- meta$cell
  message('Schirmer source total agreement: ', sum(colSums(counts)==meta$UMIs), '/', ncol(counts))
  write.csv(data.frame(cell=colnames(counts), source_total=colSums(counts), published_UMIs=meta$UMIs),file.path(base,'restored_count_totals.csv'),row.names=FALSE)
  stopifnot(all(colSums(counts)>0))
  atomic_save(list(counts=counts, meta=meta, genes=genes),file.path(base,'neuron_counts.rds'))
  message('Schirmer: ',ncol(counts),' neurons; integer raw counts and cell IDs verified')
  rm(raw,counts,meta); gc()
}
base <- 'scRNAseq/Jakel_2019/data'
if (!file.exists(file.path(base,'all_counts.rds'))) {
  counts <- as(readMM(file.path(base,'restored_counts.mtx')),'CsparseMatrix')
  rownames(counts) <- readLines(file.path(base,'restored_genes.txt'))
  colnames(counts) <- readLines(file.path(base,'restored_cells.txt'))
  meta <- read.delim(gzfile(file.path(base,'GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz')), row.names=1)
  stopifnot(all(colnames(counts) %in% rownames(meta)),all(counts@x==round(counts@x)))
  meta <- meta[colnames(counts),,drop=FALSE]
  message('Jakel detected-gene agreement: ',sum(colSums(counts>0)==meta$genes),'/',ncol(counts))
  write.csv(data.frame(cell=colnames(counts),source_detected=colSums(counts>0),published_detected=meta$genes),file.path(base,'restored_detected_genes.csv'),row.names=FALSE)
  stopifnot(all(colSums(counts)>0))
  atomic_save(list(counts=counts,meta=meta),file.path(base,'all_counts.rds'))
  message('Jakel: ',ncol(counts),' nuclei; integer raw counts and cell IDs verified')
}

