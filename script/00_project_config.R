locate_project_root <- function() {
  env_root <- Sys.getenv("SMX_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(file.path(getwd(), "script", "github", "R"), winslash = "/", mustWork = FALSE)
  }

  normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE)
}

root_dir <- locate_project_root()

project_path <- function(...) {
  file.path(root_dir, ...)
}

dir_map <- list(
  ref = project_path("ref"),
  log = project_path("Log"),
  script = project_path("script"),
  data_raw = project_path("data", "raw"),
  data_processed = project_path("data", "processed"),
  data_metadata = project_path("data", "metadata"),
  res_tables = project_path("res", "tables"),
  res_models = project_path("res", "models"),
  res_qc = project_path("res", "qc"),
  manuscript_notes = project_path("manuscript", "notes"),
  manuscript_outline = project_path("manuscript", "outline"),
  figure_source = project_path("figure", "source"),
  figure_export = project_path("figure", "export")
)

dataset_manifest <- project_path("data", "metadata", "dataset_manifest.tsv")
