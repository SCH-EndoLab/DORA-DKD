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

suppressPackageStartupMessages(library(Matrix))

read_lines <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(readLines(gzfile(path), warn = FALSE))
  }
  readLines(path, warn = FALSE)
}

read_sparse_mtx <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readMM(con)
}

quantile_or_na <- function(x, probs) {
  if (length(x) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, names = FALSE))
}

manifest <- read.delim(
  project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Focus on filtered/processed matrices and the explicit GSE279086 raw fallback.
manifest_qc <- subset(
  manifest,
  dataset_id %in% c("GSE266146", "GSE279086")
)

sample_qc_list <- vector("list", nrow(manifest_qc))

for (i in seq_len(nrow(manifest_qc))) {
  row_df <- manifest_qc[i, ]
  mat <- read_sparse_mtx(row_df$matrix_path)
  barcode_n <- length(read_lines(row_df$barcodes_path))
  feature_n <- length(read_lines(row_df$features_path))

  ncount <- Matrix::colSums(mat)
  nfeature <- Matrix::colSums(mat > 0)
  nonzero_cells <- sum(ncount > 0)
  q_count <- quantile_or_na(ncount, c(0.25, 0.5, 0.75, 0.9))
  q_feature <- quantile_or_na(nfeature, c(0.25, 0.5, 0.75, 0.9))

  sample_qc_list[[i]] <- data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    n_genes = nrow(mat),
    n_barcodes = ncol(mat),
    n_nonzero_cells = nonzero_cells,
    barcode_lines = barcode_n,
    feature_lines = feature_n,
    mean_nCount = mean(ncount),
    median_nCount = q_count[2],
    q25_nCount = q_count[1],
    q75_nCount = q_count[3],
    q90_nCount = q_count[4],
    mean_nFeature = mean(nfeature),
    median_nFeature = q_feature[2],
    q25_nFeature = q_feature[1],
    q75_nFeature = q_feature[3],
    q90_nFeature = q_feature[4],
    cells_nFeature_ge_200 = sum(nfeature >= 200),
    cells_nFeature_ge_500 = sum(nfeature >= 500),
    cells_nFeature_ge_1000 = sum(nfeature >= 1000),
    cells_nCount_ge_500 = sum(ncount >= 500),
    cells_nCount_ge_1000 = sum(ncount >= 1000),
    cells_nCount_ge_2000 = sum(ncount >= 2000),
    stringsAsFactors = FALSE
  )

  rm(mat, ncount, nfeature)
  gc(verbose = FALSE)
}

sample_qc <- do.call(rbind, sample_qc_list)
sample_qc <- sample_qc[order(sample_qc$dataset_id, sample_qc$sample_id, sample_qc$matrix_id), ]
row.names(sample_qc) <- NULL

dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)
write.table(
  sample_qc,
  file = project_path("res", "qc", "single_cell", "single_cell_filtered_matrix_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

group_key <- paste(sample_qc$dataset_id, sample_qc$matrix_format, sample_qc$group_std, sep = "||")
dataset_qc_list <- lapply(split(sample_qc, group_key), function(df) {
  data.frame(
    dataset_id = df$dataset_id[1],
    matrix_format = df$matrix_format[1],
    group_std = df$group_std[1],
    n_matrices = nrow(df),
    total_barcodes = sum(df$n_barcodes),
    total_nonzero_cells = sum(df$n_nonzero_cells),
    median_barcodes_per_matrix = median(df$n_barcodes),
    median_nonzero_cells_per_matrix = median(df$n_nonzero_cells),
    median_nCount = median(df$median_nCount),
    median_nFeature = median(df$median_nFeature),
    cells_nFeature_ge_200 = sum(df$cells_nFeature_ge_200),
    cells_nFeature_ge_500 = sum(df$cells_nFeature_ge_500),
    cells_nFeature_ge_1000 = sum(df$cells_nFeature_ge_1000),
    stringsAsFactors = FALSE
  )
})
dataset_qc <- do.call(rbind, dataset_qc_list)
dataset_qc <- dataset_qc[order(dataset_qc$dataset_id, dataset_qc$group_std, dataset_qc$matrix_format), ]

write.table(
  dataset_qc,
  file = project_path("res", "qc", "single_cell", "single_cell_filtered_dataset_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Filtered single-cell matrix QC written to res/qc/single_cell/single_cell_filtered_matrix_qc.tsv\n")
cat("Filtered single-cell dataset QC written to res/qc/single_cell/single_cell_filtered_dataset_qc.tsv\n")
