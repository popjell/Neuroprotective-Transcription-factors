

  data_dir <- "nichenet_data"
  organism <- "mouse"
  if(!dir.exists("nichenet_data")) {
  dir.create("nichenet_data", showWarnings = FALSE)
    options(timeout = 3600) # Sets timeout to 1 hour
    urls <- c(
      "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds",
      "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds",
      "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
    )

    for (u in urls) {
      dest <- file.path("nichenet_data", basename(u))
      if (!file.exists(dest)) {
        download.file(u, destfile = dest, mode = "wb", timeout = 3600)
      }
    }
  }

  if (organism == "human") {
    lr_network <- readRDS(file.path(data_dir, "lr_network_human_21122021.rds"))
    ligand_target_matrix <- readRDS(file.path(data_dir, "ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks <- readRDS(file.path(data_dir, "weighted_networks_nsga2r_final.rds"))
  } else if (organism == "mouse") {
    lr_network <- readRDS(file.path(data_dir, "lr_network_mouse_21122021.rds"))
    ligand_target_matrix <- readRDS(file.path(data_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"))
    weighted_networks <- readRDS(file.path(data_dir, "weighted_networks_nsga2r_final_mouse.rds"))
  }