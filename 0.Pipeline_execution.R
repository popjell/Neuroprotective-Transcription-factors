# Setup dependencies early
library(beepr)
library(tidyverse)
library(stringr)

if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools", repos = "https://cloud.r-project.org")
}

# Helper to run scripts safely with beep error alert
run_script_safely <- function(script_path, work_dir = NULL) {
  message("\n=== Executing: ", script_path, " ===")
  
  run_code <- function() {
    target_file <- if (!is.null(work_dir)) basename(script_path) else script_path
    source(target_file, local = FALSE, print.eval = TRUE)
  }
  
  tryCatch({
    if (!is.null(work_dir)) {
      if (requireNamespace("withr", quietly = TRUE)) {
        withr::with_dir(work_dir, run_code())
      } else {
        old_wd <- getwd()
        on.exit(setwd(old_wd))
        setwd(work_dir)
        run_code()
      }
    } else {
      run_code()
    }
  }, error = function(e) {
    beep(10)
    stop("CRITICAL: Script '", script_path, "' failed.\nError details: ", conditionMessage(e), call. = FALSE)
  })
}

# 1. Package Installation
install_script <- "install_packages.R"
if (file.exists(install_script)) {
  run_script_safely(install_script)
} else {
  beep(10)
  stop("CRITICAL: '", install_script, "' not found in root.", call. = FALSE)
}

# 2. iRegulon Reformatting
iregulon_refomat <- "iregulon_reformat.R"
if (file.exists(iregulon_refomat)) {
  run_script_safely(iregulon_refomat)
} else {
  beep(10)
  stop("CRITICAL: '", iregulon_refomat, "' not found in root.", call. = FALSE)
}

# 3. Sub-directory execution (Scripts 1-4)
sub_dirs <- c("Motor Neuron Analysis", "RGC analysis")

for (dir in sub_dirs) {
  if (!dir.exists(dir)) {
    warning("Directory '", dir, "' not found. Skipping.")
    next
  }
  
  sub_scripts <- list.files(path = dir, pattern = "^[1-4]\\..*\\.[Rr](md)?$", full.names = TRUE)
  sub_scripts <- gtools::mixedsort(sub_scripts)
  
  for (script in sub_scripts) {
    run_script_safely(script, work_dir = dir)
  }
}

# 4. Root scripts (> 4)
root_scripts <- list.files(".", pattern = "\\.[Rr]$|\\.[Rr]md$", full.names = TRUE) %>% 
  .[suppressWarnings(as.numeric(str_extract(basename(.), "^\\d+"))) > 4 %in% TRUE] %>% 
  str_sort(numeric = TRUE)

for (script in root_scripts) {
  run_script_safely(script)
}

beep(1)
message("\nFull pipeline execution complete.")