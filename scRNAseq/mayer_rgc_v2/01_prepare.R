# Step 1: check labels, map the author groups, reproduce the baseline pooled PSAP, and cache the pseudobulks (retina, retina x 22-label, retina x 27-label, retina x group).
source("common_v2.R")
suppressPackageStartupMessages(library(DESeq2))

no_overwrite(file.path(out_root, "pseudobulk.rds"))
cnt <- load_counts_ro()
x <- cnt$x; genes <- cnt$genes
meta <- read.csv(file.path(in_dir, "RGC_metadata_per_cell.csv.gz"), stringsAsFactors = FALSE)
stopifnot(identical(cnt$cells, meta$cell_barcode), anyDuplicated(genes) == 0)

# 1. metadata columns and label sets
stopifnot(all(c("Sample", "Condition", "celltypes", "cellsubtypes") %in% names(meta)))
meta$retina <- paste0("S", meta$Sample)
meta$condition <- factor(ifelse(meta$Condition == "Ctrl", "Control", "MS"), levels = c("Control", "MS"))
stopifnot(all(tapply(as.character(meta$condition), meta$retina, function(v) length(unique(v))) == 1))
stopifnot(length(unique(meta$celltypes)) == 22, length(unique(meta$cellsubtypes)) == 27)
writeLines(c("celltypes (22): " , sort(unique(meta$celltypes)), "", "cellsubtypes (27):", sort(unique(meta$cellsubtypes)),
             "", "Sample ids:", sort(unique(meta$Sample)), "", "Condition values:", unique(meta$Condition)), file.path(out_root, "label_sets.txt"))
log_line("labels: celltypes 22, cellsubtypes 27, retinas ", length(unique(meta$retina)), " (", sum(!duplicated(meta$retina) & meta$condition == "Control"), " Control, ",
         sum(!duplicated(meta$retina) & meta$condition == "MS"), " MS)")

# 2. author groups in the 27-label column
want <- data.frame(author_name = c("MG-ON4", "MG-OFF2", "ipRGC", "PG-ON", "RGC26", "RGC13"),
                   group = c(rep("susceptible", 3), rep("resilient", 3)), stringsAsFactors = FALSE)
lab27 <- sort(unique(meta$cellsubtypes))
want$matched_label <- lab27[match(norm_label(want$author_name), norm_label(lab27))]
print(want)
stopifnot(!anyNA(want$matched_label))
stopifnot(setequal(want$matched_label[want$group == "resilient"], resilient_groups), setequal(want$matched_label[want$group == "susceptible"], susceptible_groups))
write.csv(want, file.path(out_root, "author_group_mapping.csv"), row.names = FALSE)
meta$group <- ifelse(meta$cellsubtypes %in% resilient_groups, "resilient", ifelse(meta$cellsubtypes %in% susceptible_groups, "susceptible", NA))
log_line("author groups matched. cells resilient ", sum(meta$group == "resilient", na.rm = TRUE), ", susceptible ", sum(meta$group == "susceptible", na.rm = TRUE))

# 4. cells per retina per 27-label, retinas by condition, donor map
tab <- table(retina = meta$retina, subtype27 = meta$cellsubtypes)
write.csv(as.data.frame.matrix(tab), file.path(out_root, "cells_per_retina_subtype27.csv"))
ret <- unique(meta[, c("retina", "condition")]); ret$n_cells <- as.integer(table(meta$retina)[ret$retina])
write.csv(ret[order(ret$condition, ret$retina), ], file.path(out_root, "retinas.csv"), row.names = FALSE)
print(table(ret$condition))
if (file.exists(file.path(in_dir, "donor_map.csv")) || file.exists("_input/donor_map.csv")) stop("donor_map.csv exists: Analysis J should run; stop and handle") else log_line("no donor map; donor collapse skipped")

# 3. reproduce the baseline pooled DESeq2
pb_retina <- sum_counts(x, factor(meta$retina))
rownames(pb_retina) <- genes
stopifnot(all(pb_retina == round(pb_retina)), sum(pb_retina) == sum(x))
keep <- rowSums(pb_retina >= 10) >= 2
log_line("pooled pseudobulk: ", ncol(pb_retina), " retinas; genes passing count filter ", sum(keep))
cd <- data.frame(row.names = colnames(pb_retina), condition = ret$condition[match(colnames(pb_retina), ret$retina)])
dds <- DESeq(DESeqDataSetFromMatrix(pb_retina[keep, ], cd, ~condition), quiet = TRUE)
res <- results(dds, contrast = c("condition", "MS", "Control"))
base <- read.csv(file.path(base_res, "all_RGC", "deseq2.csv"), stringsAsFactors = FALSE)
b_psap <- base$log2FC[base$gene == "PSAP"]; n_psap <- res["PSAP", "log2FoldChange"]
log_line("PSAP log2FC: baseline ", round(b_psap, 6), " rebuilt ", round(n_psap, 6), "; padj rebuilt ", round(res["PSAP", "padj"], 4))
stopifnot(round(b_psap, 3) == round(n_psap, 3))
write.csv(data.frame(gene = rownames(res), log2FC = res$log2FoldChange, lfcSE = res$lfcSE, pvalue = res$pvalue, padj = res$padj), file.path(out_root, "pooled_deseq2_rebuilt.csv"), row.names = FALSE)

# unit tables and pseudobulks
unit_pb <- function(label_col, label_name) {
  u <- paste(meta$retina, meta[[label_col]], sep = "|")
  pb <- sum_counts(x, factor(ifelse(is.na(meta[[label_col]]), NA, u)))
  rownames(pb) <- genes
  parts <- do.call(rbind, strsplit(colnames(pb), "|", fixed = TRUE))
  n <- as.integer(table(u)[colnames(pb)])
  info <- data.frame(unit = colnames(pb), retina = parts[, 1], label = parts[, 2], n_cells = n,
                     condition = ret$condition[match(parts[, 1], ret$retina)], stringsAsFactors = FALSE)
  log_line(label_name, " units: ", nrow(info), "; with >= ", min_cells, " cells: ", sum(info$n_cells >= min_cells))
  list(pb = pb, info = info)
}
pb22 <- unit_pb("celltypes", "retina x celltypes (22)")
pb27 <- unit_pb("cellsubtypes", "retina x cellsubtypes (27)")
pbg <- unit_pb("group", "retina x author group")
prop22 <- prop.table(table(meta$retina, meta$celltypes), 1)
saveRDS(list(pb_retina = pb_retina, pb22 = pb22, pb27 = pb27, pbg = pbg, retina = ret, genes = genes, prop22 = unclass(prop22), prop27 = unclass(prop.table(tab, 1)),
             meta_small = meta[, c("cell_barcode", "retina", "condition", "celltypes", "cellsubtypes", "group")]), file.path(out_root, "pseudobulk.rds"))
log_line("step 1 passed")
