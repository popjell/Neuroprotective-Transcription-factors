#GEO Plot
#region
#Export Ensemble IDs for gprofiler use
write(intersected_genes, "intersected_upregulated_IDs.txt")
print(paste("Number of overlapping genes:", length(intersected_genes)))

library(readxl)

set <- readxl::read_xlsx("List of 97 common genes.xlsx")

# Extract the gene names as a vector
set_genes <- set$Genes

# Now the intersection will work
overlap <- intersect(intersected_symbols, set_genes)
overlap_table <- tibble::tibble(Genes = overlap)

print(paste("Number of overlapping genes with the paper's list:", nrow(overlap_table)))



#Run gprofiler analysis in web

#Import the .csv from gprofiler into R
gprofiler_results <- read.csv("GO_results.csv")

gprofiler_selected <- gprofiler_results %>%
dplyr::filter(source == "GO:BP") %>%
arrange(desc(negative_log10_of_adjusted_p_value))%>%
mutate(
term_name = reorder(term_name, -negative_log10_of_adjusted_p_value)
) %>%
slice_head(n = 20) 



ggplot(gprofiler_selected, aes(x = term_name, y = negative_log10_of_adjusted_p_value)) +
  geom_point(aes(size = intersection_size), color = "#f8766d") +
  scale_size_continuous(range = c(3, 10), name = "Count") +
    labs(
    x = NULL,
    y = expression(-log[10] ~ "adjusted pvalue"),
    title = "Top 20 Enriched GO:BP Terms",
    subtitle = "EAE upregulated genes in motor neurons",
    caption = paste0("produced on ", Sys.time())
  ) +
  coord_flip() + 
  theme_bw() +
  theme(
    # 2. FIX TEXT ALIGNMENT: Horizontal text on the Y-axis, clean text on X
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


