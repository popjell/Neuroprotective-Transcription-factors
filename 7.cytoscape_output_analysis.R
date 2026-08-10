#note: run centiscape analysis and add in stringdb interactions to the cytoscape output table and export it
library(dplyr)
library(tidyr)
library(readr)
if(!file.exists("iregulon+centiscape+stringdb_node_table.csv")){
  stop("iregulon+centiscape+stringdb_node_table.csv not found in the root directory. Please run the centiscape analysis and add stringdb interactions to the cytoscape output table before proceeding.")
}
print(colnames(read_csv("iregulon+centiscape+stringdb_node_table.csv")))
cytoscape_output_full <- read_csv("iregulon+centiscape+stringdb_node_table.csv") %>%
  dplyr::rename(
    symbol = `merge`
  ) %>%
  # Fix the NA values before doing the final column selection
  dplyr::mutate(
    InDegree = tidyr::replace_na(InDegree, 0),
    OutDegree = tidyr::replace_na(OutDegree, 0)
  ) %>%
  dplyr::select(
    #columns to remove: 
    -`@id`,
    -`display name`,
    -`name`,
    -`query term`,
    -`shared name`
  ) %>%
  dplyr::select(
    #regular ordering
    symbol,
    ensembl_id,
    `Regulatory function`,
    EAE_Combined_Score,
    Specificity,
    InDegree,
    OutDegree,
    
    starts_with("log2FoldChange"),
    starts_with("pvalue"),
    starts_with("padj"),
    starts_with("stat"),
    starts_with("lfcSE"),
    
    starts_with("baseMean"),
    
    starts_with("is_sig"),
    selected,
    
    starts_with("compartment::"),
    starts_with("tissue::"),
    
    starts_with("stringdb::"),
    
    everything()
  )


cytoscape_output <- cytoscape_output_full %>%
  dplyr::select(
    symbol,
    ensembl_id,
    `Regulatory function`,
    EAE_Combined_Score,
    Specificity,
    InDegree,
    OutDegree,
    
    starts_with("log2FoldChange"),
    starts_with("pvalue"),
    starts_with("padj"),
    starts_with("stat"),
    starts_with("lfcSE"),
    
    starts_with("baseMean"),
    
    starts_with("is_sig"),
  )

#View(cytoscape_output)
#ranking TF's by centrality analysis

#1. calculate avg degree (outdegree / total outdegee), and (indegree/ total indegree)
cytoscape_output <- cytoscape_output %>%
  dplyr::mutate(
    avg_outdegree = OutDegree / sum(OutDegree),
    avg_indegree = InDegree / sum(InDegree)
  )

tf_ranking <- cytoscape_output %>%
  dplyr::filter(`Regulatory function` == "Regulator") %>%
  dplyr::arrange(desc(avg_outdegree))

ggplot(tf_ranking, aes(x = reorder(symbol, -avg_outdegree), y = avg_outdegree))+
  geom_bar(stat = "identity", fill = "mediumpurple") +
  labs(title = "Average OutDegree of Transcription Factors",
       subtitle = "EAE upregulated genes in motor neurons and RGC's",
       x = "Transcription Factor",
        y = "Average OutDegree",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

top_targets <- cytoscape_output %>%
  dplyr::slice_max(order_by = avg_indegree, n = 20)

ggplot(top_targets, aes(x = reorder(symbol, -avg_indegree), y = avg_indegree)) +
  geom_bar(stat = "identity", fill = "coral") +
  labs(title = "Top 20 Most Heavily Co-Regulated Proteins",
       x = "Regulated Protein",
       y = "Average In-Degree") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





tf_data <- cytoscape_output_full %>% 
  dplyr::filter(`Regulatory function` == "Regulator")

ggplot(tf_data, aes(x = log2FoldChange_MN, y = log2FoldChange_RGC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = OutDegree), color = "#6a3d9a", alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = symbol), size = 3) +
  labs(title = "Transcription Factor Expression & Regulatory Impact",
       x = "Log2 Fold Change (Motor Neurons)",
       y = "Log2 Fold Change (RGCs)",
       size = "Out-Degree (Targets)") +
  theme_bw()

