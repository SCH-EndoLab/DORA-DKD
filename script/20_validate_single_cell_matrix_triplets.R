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

read_text_lines <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(readLines(gzfile(path), warn = FALSE))
  }
  readLines(path, warn = FALSE)
}

read_mtx_dims <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) == 0 || !grepl("^%%MatrixMarket", header[1])) {
    stop("Invalid Matrix Market header")
  }

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) {
      stop("Matrix Market size line not found")
    }
    if (!startsWith(line, "%")) {
      parts <- strsplit(trimws(line), "[[:space:]]+")[[1]]
      return(as.integer(parts[1:2]))
    }
  }
}

manifest <- read.delim(
  project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest_ready <- manifest[manifest$files_ready, ]
manifest_ready <- manifest_ready[order(manifest_ready$priority_rank, manifest_ready$dataset_id, manifest_ready$sample_id, manifest_ready$matrix_id), ]

validation_subset <- do.call(
  rbind,
  lapply(split(manifest_ready, manifest_ready$dataset_id), function(df) {
    head(df, 3)
  })
)

validate_one <- function(row_df) {
  matrix_dims <- tryCatch(
    read_mtx_dims(row_df$matrix_path),
    error = function(e) c(NA_integer_, NA_integer_)
  )
  feature_n <- tryCatch(length(read_text_lines(row_df$features_path)), error = function(e) NA_integer_)
  barcode_n <- tryCatch(length(read_text_lines(row_df$barcodes_path)), error = function(e) NA_integer_)

  data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    matrix_id = row_df$matrix_id,
    matrix_nrow = matrix_dims[1],
    matrix_ncol = matrix_dims[2],
    feature_lines = feature_n,
    barcode_lines = barcode_n,
    row_match = identical(as.integer(matrix_dims[1]), as.integer(feature_n)),
    col_match = identical(as.integer(matrix_dims[2]), as.integer(barcode_n)),
    stringsAsFactors = FALSE
  )
}

validation_list <- lapply(seq_len(nrow(validation_subset)), function(i) validate_one(validation_subset[i, ]))
validation_df <- do.call(rbind, validation_list)

dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)
write.table(
  validation_df,
  file = project_path("res", "qc", "single_cell", "single_cell_matrix_validation.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell matrix validation written to res/qc/single_cell/single_cell_matrix_validation.tsv\n")
