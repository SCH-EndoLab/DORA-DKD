script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(file.path(getwd(), "script", "github", "R"), winslash = "/", mustWork = FALSE)
  }
})
source(file.path(script_dir, "00_project_config.R"))

manifest <- read.delim(
  dataset_manifest,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

cat("Dataset manifest summary\n")
cat("========================\n")
cat("Rows:", nrow(manifest), "\n")
cat("Layers:\n")
print(table(manifest$layer))
cat("\nStatuses:\n")
print(table(manifest$status))
