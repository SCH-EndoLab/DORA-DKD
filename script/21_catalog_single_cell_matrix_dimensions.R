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
    return(length(readLines(gzfile(path), warn = FALSE)))
  }
  length(readLines(path, warn = FALSE))
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
      return(as.integer(parts[1:3]))
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

catalog_list <- lapply(seq_len(nrow(manifest)), function(i) {
  row_df <- manifest[i, ]
  dims <- read_mtx_dims(row_df$matrix_path)
  feature_n <- read_text_lines(row_df$features_path)
  barcode_n <- read_text_lines(row_df$barcodes_path)

  data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    n_genes = dims[1],
    n_barcodes = dims[2],
    n_nonzero = dims[3],
    feature_lines = feature_n,
    barcode_lines = barcode_n,
    row_match = dims[1] == feature_n,
    col_match = dims[2] == barcode_n,
    stringsAsFactors = FALSE
  )
})

catalog <- do.call(rbind, catalog_list)

dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)
write.table(
  catalog,
  file = project_path("res", "qc", "single_cell", "single_cell_matrix_dimension_catalog.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_df <- aggregate(
  list(
    n_matrices = rep(1L, nrow(catalog)),
    median_genes = catalog$n_genes,
    median_barcodes = catalog$n_barcodes,
    median_nonzero = catalog$n_nonzero
  ),
  by = list(
    dataset_id = catalog$dataset_id,
    matrix_format = catalog$matrix_format,
    analysis_role = catalog$analysis_role
  ),
  FUN = function(x) {
    if (all(x %in% c(0L, 1L))) {
      sum(x)
    } else {
      median(x)
    }
  }
)

write.table(
  summary_df,
  file = project_path("res", "qc", "single_cell", "single_cell_matrix_dimension_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell matrix dimension catalog written to res/qc/single_cell/single_cell_matrix_dimension_catalog.tsv\n")
cat("Single-cell matrix dimension summary written to res/qc/single_cell/single_cell_matrix_dimension_summary.tsv\n")
