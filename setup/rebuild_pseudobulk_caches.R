# Reuse each analysis's own count aggregation, stopping before statistical tests.
project_root <- getwd()
scripts <- c(Trobisch_2022='09_integrated_neurons.R',
             Trobisch_2022='12_cortex_integrated_neurons.R',
             Schirmer_2019='02_pseudobulk_deseq2.R',
             Jakel_2019='02_pseudobulk_deseq2.R',
             Fagiani_2025='02_pseudobulk_deseq2.R')
for (i in seq_along(scripts)) {
  setwd(file.path(project_root,'scRNAseq',names(scripts)[i]))
  expressions <- parse(scripts[i])
  save_index <- which(vapply(expressions,function(e) {
    is.call(e) && identical(e[[1]],as.name('saveRDS'))
  },logical(1)))[1]
  stopifnot(!is.na(save_index))
  destination <- as.character(expressions[[save_index]][[3]])
  if (file.exists(destination)) next
  env <- new.env(parent=globalenv())
  for (j in seq_len(save_index)) {
    expr <- expressions[[j]]
    if (is.call(expr) && identical(expr[[1]],as.name('write.csv'))) next
    eval(expr,envir=env)
  }
  message('Restored ',names(scripts)[i],': ',destination)
  rm(env); gc()
}
setwd(project_root)
