# Secondary: within-donor GML vs NAGM in upper_EN, one GML and one NAGM pseudobulk per MS donor, DESeq2 ~donor + tissue (NAGM reference).
# Needs results/pseudobulk_sample_group.rds and results/samples.csv from 00-01, and FROZEN_PLAN.md. Local pathology, not MS vs control.
source("common.R")

out_dir <- res("within_donor_GML_NAGM")
dir.create(out_dir, showWarnings = FALSE)
out <- file.path(out_dir, "psap.csv")
no_overwrite(out)

# Samples with at least min_cells upper_EN cells, from MS donors holding both tissues
pb <- load_pb()
up <- pb$cols[pb$cols$group == "upper_EN" & pb$cols$condition == "MS" & pb$cols$n_cells >= min_cells, ]
both <- intersect(up$donor[up$tissue == "GML"], up$donor[up$tissue == "NAGM"])
stopifnot(length(both) == 10, length(both) >= min_donors)
smp <- up[up$donor %in% both, c("sample", "donor", "tissue", "n_cells", "id")]

# Sequencing pool of each sample
pool <- unique(read.csv(gzfile(file.path(macnair_dir, "ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz")))[, c("sample_id_anon", "seq_pool")])
stopifnot(!anyDuplicated(pool$sample_id_anon))
smp$seq_pool <- pool$seq_pool[match(smp$sample, pool$sample_id_anon)]
smp$pseudobulk <- paste(smp$donor, smp$tissue, sep = "|")

# One pseudobulk per donor x tissue
s <- unique(smp[, c("pseudobulk", "donor", "tissue")])
s$n_cells <- tapply(smp$n_cells, smp$pseudobulk, sum)[s$pseudobulk]
s$tissue <- factor(s$tissue, levels = c("NAGM", "GML"))
rownames(s) <- s$pseudobulk
mat <- pseudobulk(pb$counts[, smp$id], smp$pseudobulk)[, s$pseudobulk]

# DESeq2 ~donor + tissue, NAGM reference
fit <- fit_deseq(mat, s, ~donor + tissue, c("tissue", "GML", "NAGM"))
psap <- cbind(dataset = "Macnair GM upper_EN, GML vs NAGM within MS donors (local pathology, not MS vs control)", donors = length(both), psap_row(fit$res))
psap$direction <- ifelse(psap$log2FC > 0, "higher in GML", "lower in GML")
print(psap)
write.csv(psap, out, row.names = FALSE)
write.csv(fit$res, file.path(out_dir, "deseq2.csv"), row.names = FALSE)
write.csv(smp[order(smp$donor, smp$tissue, smp$sample), c("donor", "tissue", "sample", "seq_pool", "n_cells", "pseudobulk")], file.path(out_dir, "samples.csv"), row.names = FALSE)
cat("within-donor:", length(both), "donors,", nrow(smp), "samples summed into", nrow(s), "pseudobulks\n")
