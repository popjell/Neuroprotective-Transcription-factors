expression_data <- readr::read_table("Expression_data.txt", col_names = TRUE)

expression_data


files <- c("GO_cell_death", "GO_metabolism", "LLM_cell_death", "LLM_metabolism")
genelists <- list()
for (f in files) {
  genes <- unique(readLines(file.path("filtered_genelists", paste0(f, "_genes.txt"))))
  genelists[[f]] <- expression_data[expression_data$symbol %in% genes, ]
  assign(paste0(f, "_genes"), genelists[[f]])
  write.csv(genelists[[f]], file.path("filtered_genelists", paste0(f, "_genes_expression_data.csv")), row.names = FALSE)
}


