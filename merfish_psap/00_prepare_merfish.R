# Read the 17 GSE284005 MERSCOPE samples into one cached object (counts, cell metadata, centroids). Needs data/raw/ and the series matrix.
source("common.R")

out <- "data/merfish_input.rds"
no_overwrite(out)

# sample metadata from the series matrix
sm <- readLines(gzfile("data/GSE284005_series_matrix.txt.gz"))
field <- function(key) lapply(sm[startsWith(sm, key)], function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
ch <- field("!Sample_characteristics_ch1")
samples <- data.frame(sample = field("!Sample_title")[[1]], gsm = field("!Sample_geo_accession")[[1]],
                      tissue = sub("tissue: ", "", ch[[1]]), sex = sub("gender: ", "", ch[[2]]),
                      ms_duration_years = as.numeric(sub("ms duration: ", "", ch[[3]])), age = as.numeric(sub("age: ", "", ch[[4]])))
# the 3 samples with duration 0 are the ones the authors label healthyWM; treated as controls (inferred, not stated in GEO)
samples$condition <- ifelse(samples$ms_duration_years == 0, "Control", "MS")
# no donor ID is given; every sample differs in sex, age and duration, so one donor per sample
samples$donor <- samples$sample

# read per-sample counts, labels and centroids
counts <- list()
meta <- list()
for (i in seq_len(nrow(samples))) {
  pre <- file.path("data/raw", paste0(samples$gsm[i], "_", samples$sample[i]))
  cnt <- fread(paste0(pre, "_count.tsv.gz"), data.table = FALSE)
  ct <- fread(paste0(pre, "_celltypes.tsv.gz"), data.table = FALSE)
  xy <- fread(paste0(pre, "_coordinates.tsv.gz"), data.table = FALSE)
  m <- as(as.matrix(cnt[, -1]), "CsparseMatrix")
  rownames(m) <- cnt[[1]]
  stopifnot(identical(colnames(m), ct$cells), identical(ct$cells, xy$cells))
  counts[[i]] <- m
  meta[[i]] <- data.frame(cell = ct$cells, sample = samples$sample[i], major = ct$majorCluster_final, region = ct$Region_banksy,
                          sub = ct$clean_sub, x = xy$X, y = xy$Y)
}
stopifnot(all(vapply(counts, function(m) identical(rownames(m), rownames(counts[[1]])), TRUE)))
counts <- do.call(cbind, counts)
meta <- merge(do.call(rbind, meta), samples, by = "sample", sort = FALSE)
meta <- meta[match(colnames(counts), meta$cell), c("cell", "sample", "donor", "condition", "tissue", "sex", "age",
                                                    "ms_duration_years", "major", "region", "sub", "x", "y")]
stopifnot(identical(meta$cell, colnames(counts)))
saveRDS(list(counts = counts, meta = meta, samples = samples), out)
cat(nrow(samples), "samples,", nrow(counts), "genes,", ncol(counts), "cells\n")
