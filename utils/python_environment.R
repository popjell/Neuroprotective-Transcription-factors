# Run from the repository root, as expected by the analysis scripts.
project_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
if (!nzchar(project_python)) {
  project_python <- file.path(".venv", if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
}
if (!file.exists(project_python)) {
  stop("Project Python environment missing: ", project_python,
       ". Run setup/install_wsl.sh from the repository root.")
}
reticulate::use_python(normalizePath(project_python, mustWork = TRUE), required = TRUE)
