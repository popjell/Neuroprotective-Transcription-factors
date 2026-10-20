# Stage 1 tables: exclusions, donor eligibility cross-tabs, GML/NAGM pairs and input_inventory.md. No gene-level results are read.
# Needs results/labels.csv, results/samples.csv and results/macnair_dims.txt from 00-01.
source("common.R")

out_md <- res("input_inventory.md")
no_overwrite(out_md)
dir.create(res("stage1_tables"), showWarnings = FALSE)
tab <- function(...) file.path(res("stage1_tables"), ...)

samples <- read.csv(res("samples.csv"))
lab <- read.csv(res("labels.csv"))
stopifnot(nrow(lab) == sum(samples$n_neurons), setequal(lab$sample, samples$sample))

# Overlap with Trobisch / Schirmer donors by ID (Macnair IDs are anonymised, so a match is only possible by coincidence)
schirmer <- read.csv(schirmer_donor_file)
macnair_ids <- unique(c(samples$donor, samples$sample))
hits <- macnair_ids[macnair_ids %in% c(shared_donors, schirmer$donor, schirmer$sample)]
excluded <- data.frame(donor = character(), reason = character())
if (length(hits)) excluded <- data.frame(donor = hits, reason = "ID matches a Trobisch/Schirmer donor or sample")
write.csv(excluded, res("excluded_donors.csv"), row.names = FALSE)
meta_names <- names(read.csv(gzfile(file.path(macnair_dir, "ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz")), nrows = 1))
source_cols <- grep("study|dataset|origin|schirmer|jakel|absinta", meta_names, ignore.case = TRUE, value = TRUE)

# Per-sample group counts
per_sample <- as.data.frame.matrix(table(lab$sample, factor(lab$coarse_group, c("upper_EN", "deep_EN", "IN", "excluded", "unassigned"))))
per_sample <- cbind(sample = rownames(per_sample), per_sample)
per_sample <- merge(samples[, c("sample", "donor", "condition", "tissue", "bank", "sex", "age", "pmi_minutes")], per_sample)
write.csv(per_sample, tab("per_sample_groups.csv"), row.names = FALSE)
write.csv(samples, tab("samples.csv"), row.names = FALSE)

# Donors with at least min_cells pooled upper_EN cells
donors <- aggregate(upper_EN ~ donor + condition + bank + sex + age + pmi_minutes, per_sample, sum)
stopifnot(!anyDuplicated(donors$donor))
eligible <- donors[donors$upper_EN >= min_cells, ]
write.csv(donors, tab("donors_upper_EN.csv"), row.names = FALSE)
by_bank <- as.data.frame.matrix(table(condition = eligible$condition, bank = eligible$bank))
by_sex <- as.data.frame.matrix(table(condition = eligible$condition, sex = eligible$sex))
write.csv(cbind(condition = rownames(by_bank), by_bank), tab("condition_by_bank.csv"), row.names = FALSE)
write.csv(cbind(condition = rownames(by_sex), by_sex), tab("condition_by_sex.csv"), row.names = FALSE)
summ <- function(x) c(n = length(x), median = median(x), min = min(x), max = max(x))
cov <- do.call(rbind, lapply(split(eligible, eligible$condition), function(e)
  data.frame(condition = e$condition[1], rbind(age = summ(e$age), pmi_minutes = summ(e$pmi_minutes)), variable = c("age", "pmi_minutes"))))
write.csv(cov, tab("age_pmi_by_condition.csv"), row.names = FALSE)


# MS donors with a GML and an NAGM sample that each have at least min_cells upper_EN cells
ms <- per_sample[per_sample$condition == "MS", ]
pairs <- do.call(rbind, lapply(split(ms, ms$donor), function(m) data.frame(donor = m$donor[1],
  gml_samples = sum(m$tissue == "GML" & m$upper_EN >= min_cells), nagm_samples = sum(m$tissue == "NAGM" & m$upper_EN >= min_cells))))
pairs$paired <- pairs$gml_samples > 0 & pairs$nagm_samples > 0
write.csv(pairs, tab("gml_nagm_pairs.csv"), row.names = FALSE)

# Inventory
md <- function(x) paste(capture.output(print(x, row.names = FALSE)), collapse = "\n")
dims <- readLines(res("macnair_dims.txt"))
writeLines(c("# Input inventory (Macnair)", "",
  "- source: Macnair et al., Zenodo 8338963; the 7 files below were moved (mv, not copied) into macnair_harmonized/data/; md5 checked against the verified json files",
  "- original location of the 7 moved files: /home/jathav/lab/.tmp_macnair_neurons_20261010/data/",
  paste0("  - ", c("ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz", "ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz",
                  "ms_lesions_snRNAseq_col_data_column_descriptions_2023-09-12.txt", "ms_lesions_snRNAseq_row_data_2023-09-12.txt.gz",
                  "ms_lesions_snRNAseq_scripts_2023-09-12.zip", "matrix_download_verified.json", "metadata_download_verified.json")),
  paste0("- matrix: ", normalizePath(file.path(macnair_dir, "ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz")), " (", round(file.size(file.path(macnair_dir, "ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz")) / 1e9, 2), " GB)"),
  paste0("- ", dims), paste0("- metadata columns: ", paste(meta_names, collapse = ", ")),
  "- tissue type: matter == GM and lesion_type GM = control GM, NAGM, GML (diagnosis CTR only in control GM; MS only in NAGM/GML)",
  "- neurons: Macnair type_broad Excitatory neurons + Inhibitory neurons; flagged exclude_pseudobulk cells are kept (counts in stage1_tables/samples.csv)",
  paste0("- source-dataset columns: ", if (length(source_cols)) paste(source_cols, collapse = ", ") else "none (sample_source is the brain bank, not a dataset)"),
  paste0("- ID matches with Trobisch/Schirmer lists: ", if (length(hits)) paste(hits, collapse = ", ") else "none (IDs are anonymised IndNN, so overlap cannot be ruled out by ID)"),
  "- label transfer as schirmer_harmonized/01_label_transfer.R, but with the 31 Rank 2 genes removed from reference and query (the Schirmer run did not remove them)", "",
  sprintf("Neurons %d; unassigned %d (%.1f%%)", nrow(lab), sum(lab$coarse_group == "unassigned"), 100 * mean(lab$coarse_group == "unassigned")), "",
  "## Per sample", "", "```", md(per_sample), "```", "",
  paste0("## Donors with at least ", min_cells, " pooled upper_EN cells: ", nrow(eligible), " of ", nrow(donors)), "",
  "condition x bank", "```", md(cbind(condition = rownames(by_bank), by_bank)), "```",
  "condition x sex", "```", md(cbind(condition = rownames(by_sex), by_sex)), "```",
  "age and PMI", "```", md(cov), "```", "",
  paste0("MS donors with GML and NAGM samples each >= ", min_cells, " upper_EN cells: ", sum(pairs$paired), " of ", nrow(pairs)), "```", md(pairs), "```"), out_md)
cat("eligible donors:", nrow(eligible), "; GML/NAGM donors:", sum(pairs$paired), "\n")
