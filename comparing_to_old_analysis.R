
set <- readxl::read_xlsx("List of 97 common genes.xlsx")

# Extract the gene names as a vector
set_genes <- set$Genes

# Now the intersection will work
overlap <- intersect(intersected_symbols, set_genes)
overlap_table <- tibble::tibble(Genes = overlap)

print(paste("Number of overlapping genes with the paper's 97 gene list:", nrow(overlap_table)))


excel <- read_csv("biorxiv_analysis.csv")


old_results <- excel %>%
  dplyr::filter(log2FoldChange > 0 & pvalue < 0.05)

  compare <- as_tibble(intersect(old_results$...1, RGC_upregulated_candidates$symbol))

  print(paste("Number of overlapping genes:", nrow(compare)))
  print(paste("% of my upregulated genes that are in the old's list:", round(nrow(compare) / nrow(RGC_upregulated_candidates) * 100, 2)))
  print(paste("% of the old's upregulated genes that are in my list:", round(nrow(compare) / nrow(old_results) * 100, 2)))






  #change below to use whole dataset rather than just the upregulated 
  comparison <- inner_join(
    RGC_results.df %>% dplyr::select(symbol, pvalue, log2FoldChange) %>% 
        dplyr::rename(pvalue_mine = pvalue, lfc_mine = log2FoldChange),
          excel %>% dplyr::select(...1, pvalue, log2FoldChange) %>%
              dplyr::rename(pvalue_old = pvalue, lfc_old = log2FoldChange, symbol = ...1),
                by = "symbol"
                )
                # Plot p-value correlation
                #add lines to show significance thresholds, and red to show where they wont be matching between the 2 sets
                ggplot(comparison, aes(x = -log10(pvalue_old), y = -log10(pvalue_mine))) +
                  geom_point(alpha = 0.3, label = comparison$symbol) +
                    geom_abline(slope = 1, intercept = 0, color = "red") +
                      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "blue") +
                        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
                          labs(x = "-log10(p-value) old paper", y = "-log10(p-value) mine") +
                            theme_bw() +
                              labs(title = "P-value comparison: mine vs old")

                                  # And LFC correlation
                                  ggplot(comparison, aes(x = lfc_old, y = lfc_mine)) +
                                    geom_point(alpha = 0.3) +
                                      geom_abline(slope = 1, intercept = 0, color = "red") +
                                        theme_bw() +
                                          labs(title = "Log2FC comparison: mine vs old")
