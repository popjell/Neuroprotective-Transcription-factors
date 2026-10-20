library(Matrix)
rows <- list()
for (study in c('Trobisch_2022','Schirmer_2019','Jakel_2019','Fagiani_2025')) {
  folder <- file.path('scRNAseq',study,'data')
  caches <- switch(study,
    Trobisch_2022=c('new_regions.rds','cortex_neurons.rds','cortex_senders.rds'),
    Schirmer_2019='neuron_counts.rds',Jakel_2019='all_counts.rds',Fagiani_2025='neuron_counts.rds')
  for (cache in caches) {
    x <- readRDS(file.path(folder,cache))
    counts <- if (is(x,'sparseMatrix')) x else x$counts
    meta <- if (study=='Fagiani_2025') readRDS(file.path(folder,'neuron_metadata.rds')) else x$meta
    stopifnot(ncol(counts)==nrow(meta),!anyDuplicated(colnames(counts)),
              all(counts@x>=0),all(counts@x==round(counts@x)))
    if (study=='Fagiani_2025') {
      stopifnot(identical(colnames(counts),meta$barcodes),
                all(colSums(counts)==meta$nCount_RNA),
                all(colSums(counts>0)==meta$nFeature_RNA))
      keep <- meta$nFeature_RNA>=500 & meta$nCount_RNA>=1000 & meta$percent.mt<=5
      if (!all(keep)) {
        counts <- counts[,keep,drop=FALSE]; meta <- meta[keep,]
        saveRDS(counts,file.path(folder,'neuron_counts.rds.part'))
        stopifnot(file.rename(file.path(folder,'neuron_counts.rds.part'),file.path(folder,'neuron_counts.rds')))
        saveRDS(meta,file.path(folder,'neuron_metadata.rds.part'))
        stopifnot(file.rename(file.path(folder,'neuron_metadata.rds.part'),file.path(folder,'neuron_metadata.rds')))
      }
      old <- read.csv(file.path('scRNAseq',study,'results/cell_metadata.csv'))
      stopifnot(all(old$barcodes %in% meta$barcodes))
      if (!setequal(meta$barcodes,old$barcodes)) {
        selection <- match(old$barcodes,meta$barcodes)
        counts <- counts[,selection,drop=FALSE]; meta <- meta[selection,]
        saveRDS(counts,file.path(folder,'neuron_counts.rds.part'))
        stopifnot(file.rename(file.path(folder,'neuron_counts.rds.part'),file.path(folder,'neuron_counts.rds')))
        saveRDS(meta,file.path(folder,'neuron_metadata.rds.part'))
        stopifnot(file.rename(file.path(folder,'neuron_metadata.rds.part'),file.path(folder,'neuron_metadata.rds')))
        message('Fagiani cache aligned to the committed prior-analysis cell selection: ',nrow(meta))
      }
    }
    if (study=='Trobisch_2022') stopifnot(all(colSums(counts)==meta$nCount_RNA))
    rows[[paste(study,cache)]] <- data.frame(study=study,cache=cache,genes=nrow(counts),cells=ncol(counts))
    rm(x,counts,meta); gc()
  }
}
summary <- do.call(rbind,rows)
print(summary,row.names=FALSE)
write.csv(summary,'setup/logs/restored-data-inventory.csv',row.names=FALSE)
cat('RESTORED DATA VERIFICATION PASSED\n')
