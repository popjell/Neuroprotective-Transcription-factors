# Project-local packages and Python, when R starts at the repository root.
local({
  project_root <- normalizePath(getwd())
  project_library <- file.path(project_root, ".Rlibs")
  if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
  project_python <- file.path(project_root, ".venv", "bin", "python")
  if (file.exists(project_python) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = project_python)
  }
})
