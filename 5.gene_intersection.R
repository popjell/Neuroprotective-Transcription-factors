#region
my_gene_list_MN <- MN_upregulated_candidates$ensembl_id
my_gene_list_RGC <- RGC_upregulated_candidates$ensembl_id
intersected_genes <- intersect(my_gene_list_MN, my_gene_list_RGC)

#intersect genes from both datasets  to make new list, and make venn diagram
library(VennDiagram)


# 2. Generate the Venn Diagram
# This creates a high-res TIFF file in your working directory matching the paper's style

# clean inputs, remove NAs and duplicates
mn <- unique(na.omit(my_gene_list_MN))
rgc <- unique(na.omit(my_gene_list_RGC))

#supress log file output
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# recreate and draw the Venn (filename=NULL returns grobs)
venn <- VennDiagram::venn.diagram(
  x = list(MNs = mn, RGCs = rgc),
  category.names = c("MNs", "RGCs"),
  filename = NULL,
  imagetype = "png",
  lwd = 2, lty = "solid",
  fill = c("#7b91c4", "#f47c7c"),
  alpha = c(0.8, 0.8),
  cex = 1.5, fontface = "bold", fontfamily = "sans",
  scaled = FALSE,
  cat.cex = 1.5, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05), cat.fontfamily = "sans",
  main = "Significantly upregulated genes", main.cex = 1.5, main.fontface = "bold",
  main.fontfamily = "sans"
)

grid::grid.newpage()
grid::grid.draw(venn)
venn



#output for iregulon analysis
intersected_df <- RGC_results.df[RGC_results.df$ensembl_id %in% intersected_genes, ]
intersected_symbols <- intersected_df$symbol
intersected_symbols <- intersected_symbols[!is.na(intersected_symbols) & intersected_symbols != ""]
print(paste("Total valid symbols found:", length(intersected_symbols)))
write(intersected_symbols, "intersected_upregulated_genes.txt")
write(intersected_genes, "intersected_upregulated_ensembl_ids.txt")
#now output a list of all genes with their log2 fold change and p-value for iregulon analysis, from both MN and RGC
intersected_combined_df <- merge(MN_results.df, RGC_results.df, by = "ensembl_id", suffixes = c("_MN", "_RGC"))
#change so that symbol column is not dupicated
intersected_combined_df <- intersected_combined_df %>%
 dplyr::mutate(symbol = ifelse(!is.na(symbol_MN) & symbol_MN != "", symbol_MN, symbol_RGC))  %>%
#remove the old symbol columns
  dplyr::select(-symbol_MN, -symbol_RGC)

#write as .txt file with no commas
write.table(intersected_combined_df, "Expression_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Intersect with 97 common
library(readxl)
