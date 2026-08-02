# Prepare data for community detection analysis
# Reads iRegulon results and creates the edges file needed by run_community.R

library(dplyr)
library(tidyr)

if (!file.exists("iregulon_results.csv")) {
  stop("iregulon_results.csv not found in current directory.")
}
if (!file.exists("intersected_upregulated_IDs.txt")) {
  stop("intersected_upregulated_IDs.txt not found.")
}

df <- read.csv("iregulon_results.csv", skip = 11, stringsAsFactors = FALSE, check.names = FALSE)
all_upregulated_genes <- readLines("intersected_upregulated_IDs.txt")
expression_data <- read.csv("Expression_data.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

expression_data <- expression_data %>%
  mutate(
    is_sig_RGC = !is.na(pvalue_RGC) & !is.na(log2FoldChange_RGC) & pvalue_RGC < 0.05 & abs(log2FoldChange_RGC) > 0.58,
    is_sig_MN  = !is.na(pvalue_MN)  & !is.na(log2FoldChange_MN)  & pvalue_MN  < 0.05 & abs(log2FoldChange_MN)  > 0.58,
    Specificity = case_when(
      is_sig_RGC & is_sig_MN & (sign(log2FoldChange_RGC) != sign(log2FoldChange_MN)) ~ "Opposite",
      is_sig_RGC & is_sig_MN ~ "Both",
      is_sig_RGC ~ "RGC_Only",
      is_sig_MN  ~ "MN_Only"
    ),
    EAE_Combined_Score = case_when(
      Specificity == "Opposite" ~ 0,
      Specificity == "Both"     ~ (log2FoldChange_RGC + log2FoldChange_MN) / 2,
      Specificity == "RGC_Only" ~ log2FoldChange_RGC,
      Specificity == "MN_Only"  ~ log2FoldChange_MN
    )
  )

edges <- df %>%
  filter(!is.na(`Transcription factor`) & `Transcription factor` != "") %>%
  select(`Transcription factor`, `Target genes`, `Motif id`, `NES`) %>%
  separate_rows(`Transcription factor`, sep = ",\\s*") %>%
  separate_rows(`Target genes`, sep = ",\\s*") %>%
  mutate(
    `Regulator Gene` = trimws(`Transcription factor`),
    `Target Gene` = trimws(`Target genes`),
    interaction = paste0("regulates via ", trimws(`Motif id`)),
    Assembly = "mm9"
  ) %>%
  filter(`Regulator Gene` != "" & `Target Gene` != "") %>%
  select(`Regulator Gene`, `Target Gene`, interaction, `Motif id`, NES, Assembly) %>%
  distinct(`Regulator Gene`, `Target Gene`, .keep_all = TRUE)

all_tfs <- unique(edges$`Regulator Gene`)

tf_nodes <- data.frame(
  id = all_tfs,
  `Regulatory function` = "Regulator",
  check.names = FALSE
) %>%
  left_join(expression_data, by = c("id" = "symbol")) %>%
  filter((is_sig_RGC | is_sig_MN) & Specificity != "Opposite")

edges <- edges %>%
  filter(`Regulator Gene` %in% tf_nodes$id)

write.csv(edges, "reformatted_edges.csv", row.names = FALSE)
cat("Created reformatted_edges.csv with", nrow(edges), "edges\n")

exp_path <- "Expression_data.txt"
if (!file.exists(exp_path)) {
  # Copy from parent if not present
  if (file.exists("../Expression_data.txt")) {
    file.copy("../Expression_data.txt", exp_path)
  }
}
cat("Data preparation complete.\n")
