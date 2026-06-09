#region
my_gene_list_MN <- MN_upregulated_candidates$ensembl_id
my_gene_list_RGC <- RGC_upregulated_candidates$ensembl_id
intersected_genes <- intersect(my_gene_list_MN, my_gene_list_RGC)

#intersect genes from both datasets  to make new list, and make venn diagram
library(VennDiagram)


# 2. Generate the Venn Diagram
# This creates a high-res TIFF file in your working directory matching the paper's style
venn <- venn.diagram(
  x = list(MNs = my_gene_list_MN, RGCs = my_gene_list_RGC),
  category.names = c("MNs", "RGCs"),
  filename = NULL,
  
  # Structural layout matching your image
  lwd = 2,
  lty = 'solid',
  fill = c("#7b91c4", "#f47c7c"), # Custom hex codes matching their blue/red palette
  alpha = c(0.8, 0.8),            # Transparency so the overlap looks blended
  
  # Text adjustments (Numbers inside the circles)
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  scaled = FALSE,                 
  
  # Category labels adjustments (MNs and RGCs text placement)
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),          
  cat.dist = c(0.05, 0.05),       
  cat.fontfamily = "sans",
  
  main = "Significantly upregulated genes",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)
grid.newpage()
grid.draw(venn)


#output for iregulon analysis
intersected_df <- RGC_results.df[RGC_results.df$ensembl_id %in% intersected_genes, ]
intersected_symbols <- intersected_df$symbol
intersected_symbols <- intersected_symbols[!is.na(intersected_symbols) & intersected_symbols != ""]
print(paste("Total valid symbols found:", length(intersected_symbols)))
write(intersected_symbols, "intersected_upregulated_genes.txt")

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

set <- readxl::read_xlsx("List of 97 common genes.xlsx")

# Extract the gene names as a vector
set_genes <- set$Genes

# Now the intersection will work
overlap <- intersect(intersected_symbols, set_genes)
overlap_table <- tibble::tibble(Genes = overlap)

print(paste("Number of overlapping genes with the paper's 97 gene list:", nrow(overlap_table)))

