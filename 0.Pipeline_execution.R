
#utility function to run whole pipeline


# 1. Run the package installation script first
install_script <- "install_packages.R"

if (file.exists(install_script)) {
  message("=== Step 0: Running Package Installation ===")
  tryCatch({
    source(install_script, local = FALSE, print.eval = TRUE)
  }, error = function(e) {
    beep(10)
    stop("CRITICAL: Package installation failed or threw an error. Halting pipeline execution.\nError details: ", conditionMessage(e))
  })
} else {
  beep(10)
  stop("CRITICAL: '", install_script, "' not found in the root directory. Halting execution.")
}

library(beepr)
if (!requireNamespace("gtools", quietly = TRUE)) {
  message("gtools package not found. Installing it now...")
  install.packages("gtools", repos = "https://cloud.r-project.org")
}

iregulon_refomat <- "iregulon_reformat.R"
if (file.exists(iregulon_refomat)) {
  message("=== Step 1: Running iRegulon Reformatting ===")
  tryCatch({
    source(iregulon_refomat, local = FALSE, print.eval = TRUE)
  }, error = function(e) {
    stop("CRITICAL: iRegulon reformatting failed or threw an error. Halting pipeline execution.\nError details: ", conditionMessage(e))
  })
} else {
  stop("CRITICAL: '", iregulon_refomat, "' not found in the root directory. Halting execution.")
}
# Define the two parallel analysis directories
sub_dirs <- c("Motor Neuron Analysis", "RGC analysis")

# 2. Run the nested scripts (1 to 4) in their respective directories
for (dir in sub_dirs) {
  message("\n=== Running pipeline for: ", dir, " ===")
  
  if (!dir.exists(dir)) {
    warning("Directory '", dir, "' not found. Skipping.")
    next
  }
  
  sub_scripts <- list.files(path = dir, pattern = "^[1-4]\\..*\\.R$", full.names = TRUE)
  sub_scripts <- gtools::mixedsort(sub_scripts)
  
  for (script in sub_scripts) {
    message("Executing: ", script)
    
    if (requireNamespace("withr", quietly = TRUE)) {
      withr::with_dir(dir, source(basename(script), local = FALSE, print.eval = TRUE))
    } else {
      old_wd <- getwd()
      setwd(dir)
      tryCatch(source(basename(script), local = FALSE, print.eval = TRUE), finally = setwd(old_wd))
    }
  }
}

# 3. Always run scripts 5 (gene intersection) and 6 (GO analysis)
message("\n=== Running scripts 5-6 (gene intersection + GO analysis) ===")

for (script in c("5.gene_intersection.R", "6.GO_analysis.R")) {
  if (file.exists(script)) {
    message("Executing: ", script)
    source(script, local = FALSE, print.eval = TRUE)
  } else {
    warning("Script not found: ", script)
  }
}

# 4. Run scripts 7-9 (Cytoscape output analysis, subnetwork generation, visualization)

if (file.exists("7.cytoscape_output_analysis.R")) {
  message("Executing: 7.cytoscape_output_analysis.R")
  source("7.cytoscape_output_analysis.R", local = FALSE, print.eval = TRUE)
}

for (script in c("8.generate_subregulatory_networks.R", "9.visualize_subnetworks.R", "10.tf_cluster_analysis.R", "11.cluster_deep_dive.R")) {
  if (file.exists(script)) {
    message("Executing: ", script)
    source(script, local = FALSE, print.eval = TRUE)
  } else {
    warning("Script not found: ", script)
  }
}

beep(1)
message("\nFull pipeline execution complete.")