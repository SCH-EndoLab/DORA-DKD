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

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

build_record <- function(dataset_id, sample_id, matrix_id, matrix_path, features_path, barcodes_path, matrix_format, source_level) {
  data.frame(
    dataset_id = dataset_id,
    sample_id = sample_id,
    matrix_id = matrix_id,
    matrix_path = matrix_path,
    features_path = features_path,
    barcodes_path = barcodes_path,
    matrix_format = matrix_format,
    source_level = source_level,
    stringsAsFactors = FALSE
  )
}

collect_gse209781 <- function() {
  matrix_paths <- Sys.glob(project_path("data", "processed", "single_cell_archives", "GSE209781", "*", "*", "*", "matrix.mtx.gz"))
  if (length(matrix_paths) == 0) {
    return(data.frame())
  }
  out <- lapply(matrix_paths, function(matrix_path) {
    sample_id <- strsplit(matrix_path, .Platform$file.sep, fixed = TRUE)[[1]]
    sample_id <- sample_id[length(sample_id) - 3]
    matrix_dir <- dirname(matrix_path)
    build_record(
      dataset_id = "GSE209781",
      sample_id = sample_id,
      matrix_id = basename(matrix_dir),
      matrix_path = matrix_path,
      features_path = file.path(matrix_dir, "features.tsv.gz"),
      barcodes_path = file.path(matrix_dir, "barcodes.tsv.gz"),
      matrix_format = "10x_mtx",
      source_level = "unpacked_tar"
    )
  })
  do.call(rbind, out)
}

collect_gse279086 <- function() {
  processed_invalid_samples <- c("GSM8561119", "GSM8561132")
  matrix_paths <- Sys.glob(project_path("data", "raw", "single_cell", "GSE279086", "*", "*_matrix_processed.mtx.gz"))
  if (length(matrix_paths) == 0) {
    return(data.frame())
  }
  out <- lapply(matrix_paths, function(matrix_path) {
    sample_dir <- dirname(matrix_path)
    sample_id <- basename(sample_dir)
    prefix <- sub("_matrix_processed\\.mtx\\.gz$", "", basename(matrix_path))
    processed_features <- file.path(sample_dir, paste0(prefix, "_features_processed.tsv.gz"))
    processed_barcodes <- file.path(sample_dir, paste0(prefix, "_barcodes_processed.tsv.gz"))
    raw_matrix <- file.path(sample_dir, paste0(prefix, "_matrix.mtx.gz"))
    raw_features <- file.path(sample_dir, paste0(prefix, "_features.tsv.gz"))
    raw_barcodes <- file.path(sample_dir, paste0(prefix, "_barcodes.tsv.gz"))

    processed_ok <- file.exists(matrix_path) &&
      file.exists(processed_features) &&
      file.exists(processed_barcodes) &&
      file.info(matrix_path)$size > 0 &&
      file.info(processed_features)$size > 0 &&
      file.info(processed_barcodes)$size > 0

    raw_ok <- file.exists(raw_matrix) &&
      file.exists(raw_features) &&
      file.exists(raw_barcodes) &&
      file.info(raw_matrix)$size > 0 &&
      file.info(raw_features)$size > 0 &&
      file.info(raw_barcodes)$size > 0

    if (processed_ok && !(sample_id %in% processed_invalid_samples)) {
      return(build_record(
        dataset_id = "GSE279086",
        sample_id = sample_id,
        matrix_id = prefix,
        matrix_path = matrix_path,
        features_path = processed_features,
        barcodes_path = processed_barcodes,
        matrix_format = "10x_mtx_processed",
        source_level = "downloaded_processed_triplet"
      ))
    }

    if (raw_ok && ((sample_id %in% processed_invalid_samples) || !processed_ok)) {
      return(build_record(
        dataset_id = "GSE279086",
        sample_id = sample_id,
        matrix_id = paste0(prefix, "_raw_fallback"),
        matrix_path = raw_matrix,
        features_path = raw_features,
        barcodes_path = raw_barcodes,
        matrix_format = "10x_mtx_raw_fallback",
        source_level = "downloaded_raw_triplet_fallback"
      ))
    }

    build_record(
      dataset_id = "GSE279086",
      sample_id = sample_id,
      matrix_id = prefix,
      matrix_path = matrix_path,
      features_path = processed_features,
      barcodes_path = processed_barcodes,
      matrix_format = "10x_mtx_processed",
      source_level = "downloaded_processed_triplet_incomplete"
    )
  })
  do.call(rbind, out)
}

collect_gse266146 <- function() {
  matrix_paths <- Sys.glob(project_path("data", "processed", "single_cell_matrices", "GSE266146", "*", "*", "mex", "matrix.mtx.gz"))
  if (length(matrix_paths) == 0) {
    return(data.frame())
  }
  out <- lapply(matrix_paths, function(matrix_path) {
    path_parts <- strsplit(matrix_path, .Platform$file.sep, fixed = TRUE)[[1]]
    sample_id <- path_parts[length(path_parts) - 3]
    matrix_id <- path_parts[length(path_parts) - 2]
    matrix_dir <- dirname(matrix_path)
    build_record(
      dataset_id = "GSE266146",
      sample_id = sample_id,
      matrix_id = matrix_id,
      matrix_path = matrix_path,
      features_path = file.path(matrix_dir, "features.tsv.gz"),
      barcodes_path = file.path(matrix_dir, "barcodes.tsv.gz"),
      matrix_format = "bd_rsec_mex",
      source_level = "expanded_nested_zip"
    )
  })
  do.call(rbind, out)
}

study_manifest <- read_tsv(project_path("data", "metadata", "single_cell_study_manifest.tsv"))

manifest_list <- list(
  collect_gse209781(),
  collect_gse279086(),
  collect_gse266146()
)
manifest_list <- manifest_list[vapply(manifest_list, nrow, integer(1)) > 0]
if (length(manifest_list) == 0) {
  stop("No single-cell matrix triplets detected. Run download and expansion steps first.")
}
matrix_manifest <- do.call(rbind, manifest_list)

matrix_manifest$files_ready <- file.exists(matrix_manifest$matrix_path) &
  file.exists(matrix_manifest$features_path) &
  file.exists(matrix_manifest$barcodes_path)

matrix_manifest <- merge(
  matrix_manifest,
  study_manifest[, c("dataset_id", "sample_id", "sample_title", "group_std", "diabetes_context", "analysis_role", "priority_rank")],
  by = c("dataset_id", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)

matrix_manifest <- matrix_manifest[order(matrix_manifest$priority_rank, matrix_manifest$dataset_id, matrix_manifest$sample_id, matrix_manifest$matrix_id), ]
row.names(matrix_manifest) <- NULL

ensure_dir(project_path("data", "metadata"))
ensure_dir(project_path("res", "qc", "single_cell"))

write.table(
  matrix_manifest,
  file = project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

matrix_counts <- aggregate(
  rep(1L, nrow(matrix_manifest)),
  by = list(
    dataset_id = matrix_manifest$dataset_id,
    matrix_format = matrix_manifest$matrix_format,
    analysis_role = matrix_manifest$analysis_role
  ),
  FUN = sum
)
names(matrix_counts)[4] <- "n_matrices"

ready_counts <- aggregate(
  as.integer(matrix_manifest$files_ready),
  by = list(
    dataset_id = matrix_manifest$dataset_id,
    matrix_format = matrix_manifest$matrix_format,
    analysis_role = matrix_manifest$analysis_role
  ),
  FUN = sum
)
names(ready_counts)[4] <- "n_ready"

matrix_summary <- merge(
  matrix_counts,
  ready_counts,
  by = c("dataset_id", "matrix_format", "analysis_role"),
  sort = FALSE
)

write.table(
  matrix_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_matrix_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell matrix manifest written to data/metadata/single_cell_matrix_manifest.tsv\n")
cat("Single-cell matrix summary written to res/qc/single_cell/single_cell_matrix_summary.tsv\n")
