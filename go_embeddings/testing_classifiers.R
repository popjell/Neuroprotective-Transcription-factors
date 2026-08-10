
#GO Based Annotations tests

go_annotated_llm <- read.csv("all_BP_terms_annotated_LLM.csv") %>%
  dplyr::select(
    term_name, metabolism_related, cell_killing_override
  ) %>%
  dplyr::mutate(
    metabolism = as.logical(metabolism_related),
    cell_killing = as.logical(cell_killing_override)
  ) %>%
  dplyr::select(
    term_name, metabolism, cell_killing
  )
py_run_file("go_embeddings/classifier.py")
go_annotated_semantic <- py$classified %>%
  dplyr::mutate(
    metabolism = ifelse((best_category == "metabolism" | best_category == "lipid_metabolism" | best_category == "mitochondria_energy" | best_category == "redox_metabolism"), TRUE, FALSE),
    cell_killing = ifelse((best_category == "cell_death" | best_category == "cell_killing"), TRUE, FALSE)
  ) %>%
  dplyr::select(
    term_name, metabolism, cell_killing
  )
compare_llm_semantic <- go_annotated_llm %>%
  dplyr::left_join(go_annotated_semantic, by = "term_name") %>%
  dplyr::mutate(
    cell_killing_conflict = ifelse(cell_killing.x != cell_killing.y, TRUE, FALSE),
    metabolism_conflict = ifelse(metabolism.x != metabolism.y, TRUE, FALSE)
  )
print(
  paste0(
    "Number of conflicts between LLM and semantic labeling in cell death: ",
    sum(compare_llm_semantic$cell_killing_conflict)
  )
)
print(
  paste0(
    "Number of conflicts between LLM and semantic labeling in metabolism: ",
    sum(compare_llm_semantic$metabolism_conflict)
  )
)
#semantic labeling is innefective
py_run_string("
import pandas as pd
input_df = r.all_BP
")
py_run_file("go_embeddings/GO_graph_traversal.py")
go_traversal <- py$traversal_df %>%
  dplyr::mutate(
    metabolism = ifelse(grepl("metabolism", category) & !grepl("exclude", category), TRUE, FALSE),
    cell_killing = ifelse(grepl("cell_death", category) & !grepl("exclude", category), TRUE, FALSE)
  ) %>%
  dplyr::select(
    term_name, metabolism, cell_killing, category
  )

compare_llm_traversal <- go_annotated_llm %>%
  dplyr::left_join(go_traversal, by = "term_name", suffix = c(".llm", ".traversal")) %>%
  dplyr::mutate(
    cell_killing_conflict = ifelse(cell_killing.llm != cell_killing.traversal, TRUE, FALSE),
    metabolism_conflict = ifelse(metabolism.llm != metabolism.traversal, TRUE, FALSE)
  )
print(
  paste0(
    "Number of conflicts between LLM and graph traversal labeling in cell death: ",
    sum(compare_llm_traversal$cell_killing_conflict)
  )
)
print(
  paste0(
    "Number of conflicts between LLM and graph traversal labeling in metabolism: ",
    sum(compare_llm_traversal$metabolism_conflict)
  )
)
