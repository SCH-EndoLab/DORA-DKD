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

read_feature_symbols <- function(path) {
  lines <- read_lines(path)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  gene_symbol <- vapply(
    split_lines,
    function(x) {
      if (length(x) >= 2) {
        x[2]
      } else {
        x[1]
      }
    },
    character(1)
  )
  trimws(gene_symbol)
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

safe_pct <- function(num, den) {
  out <- rep(0, length(den))
  keep <- den > 0
  out[keep] <- 100 * num[keep] / den[keep]
  out
}

manifest <- read.delim(
  project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest_qc <- subset(
  manifest,
  dataset_id %in% c("GSE266146", "GSE279086")
)

barcode_qc_list <- vector("list", nrow(manifest_qc))
sample_summary_list <- vector("list", nrow(manifest_qc))

for (i in seq_len(nrow(manifest_qc))) {
  row_df <- manifest_qc[i, ]
  features <- read_feature_symbols(row_df$features_path)
  barcodes <- read_lines(row_df$barcodes_path)
  mat <- read_sparse_mtx(row_df$matrix_path)

  stopifnot(length(features) == nrow(mat))
  stopifnot(length(barcodes) == ncol(mat))

  gene_upper <- toupper(features)
  mito_idx <- grepl("^MT-", gene_upper)
  ribo_idx <- grepl("^RPL|^RPS", gene_upper)
  hb_idx <- grepl("^HB(?!P)", gene_upper, perl = TRUE)

  ncount <- Matrix::colSums(mat)
  nfeature <- Matrix::colSums(mat > 0)
  mito_counts <- if (any(mito_idx)) Matrix::colSums(mat[mito_idx, , drop = FALSE]) else rep(0, ncol(mat))
  ribo_counts <- if (any(ribo_idx)) Matrix::colSums(mat[ribo_idx, , drop = FALSE]) else rep(0, ncol(mat))
  hb_counts <- if (any(hb_idx)) Matrix::colSums(mat[hb_idx, , drop = FALSE]) else rep(0, ncol(mat))

  keep_cells <- ncount > 0

  barcode_qc <- data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    barcode = barcodes[keep_cells],
    nCount = as.numeric(ncount[keep_cells]),
    nFeature = as.numeric(nfeature[keep_cells]),
    pct_mito = safe_pct(as.numeric(mito_counts[keep_cells]), as.numeric(ncount[keep_cells])),
    pct_ribo = safe_pct(as.numeric(ribo_counts[keep_cells]), as.numeric(ncount[keep_cells])),
    pct_hb = safe_pct(as.numeric(hb_counts[keep_cells]), as.numeric(ncount[keep_cells])),
    stringsAsFactors = FALSE
  )

  barcode_qc$core_pass_200_500_20 <- with(
    barcode_qc,
    nFeature >= 200 & nCount >= 500 & pct_mito <= 20
  )
  barcode_qc$strict_pass_500_1000_15 <- with(
    barcode_qc,
    nFeature >= 500 & nCount >= 1000 & pct_mito <= 15
  )

  q_count <- quantile_or_na(barcode_qc$nCount, c(0.25, 0.5, 0.75, 0.9))
  q_feature <- quantile_or_na(barcode_qc$nFeature, c(0.25, 0.5, 0.75, 0.9))
  q_mito <- quantile_or_na(barcode_qc$pct_mito, c(0.25, 0.5, 0.75, 0.9))

  sample_summary_list[[i]] <- data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    n_cells = nrow(barcode_qc),
    median_nCount = q_count[2],
    q25_nCount = q_count[1],
    q75_nCount = q_count[3],
    q90_nCount = q_count[4],
    median_nFeature = q_feature[2],
    q25_nFeature = q_feature[1],
    q75_nFeature = q_feature[3],
    q90_nFeature = q_feature[4],
    median_pct_mito = q_mito[2],
    q25_pct_mito = q_mito[1],
    q75_pct_mito = q_mito[3],
    q90_pct_mito = q_mito[4],
    cells_pct_mito_le_10 = sum(barcode_qc$pct_mito <= 10),
    cells_pct_mito_le_20 = sum(barcode_qc$pct_mito <= 20),
    cells_core_pass_200_500_20 = sum(barcode_qc$core_pass_200_500_20),
    cells_strict_pass_500_1000_15 = sum(barcode_qc$strict_pass_500_1000_15),
    stringsAsFactors = FALSE
  )

  barcode_qc_list[[i]] <- barcode_qc

  rm(mat, features, barcodes, ncount, nfeature, mito_counts, ribo_counts, hb_counts, barcode_qc)
  gc(verbose = FALSE)
}

barcode_qc_all <- do.call(rbind, barcode_qc_list)
sample_qc_summary <- do.call(rbind, sample_summary_list)

dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)

write.table(
  barcode_qc_all,
  file = gzfile(project_path("res", "qc", "single_cell", "single_cell_barcode_qc.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sample_qc_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_barcode_qc_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

group_key <- paste(sample_qc_summary$dataset_id, sample_qc_summary$matrix_format, sample_qc_summary$group_std, sep = "||")
group_summary <- do.call(
  rbind,
  lapply(split(sample_qc_summary, group_key), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      matrix_format = df$matrix_format[1],
      group_std = df$group_std[1],
      n_matrices = nrow(df),
      total_cells = sum(df$n_cells),
      median_cells_per_matrix = median(df$n_cells),
      median_nCount = median(df$median_nCount),
      median_nFeature = median(df$median_nFeature),
      median_pct_mito = median(df$median_pct_mito),
      cells_pct_mito_le_10 = sum(df$cells_pct_mito_le_10),
      cells_pct_mito_le_20 = sum(df$cells_pct_mito_le_20),
      cells_core_pass_200_500_20 = sum(df$cells_core_pass_200_500_20),
      cells_strict_pass_500_1000_15 = sum(df$cells_strict_pass_500_1000_15),
      stringsAsFactors = FALSE
    )
  })
)
group_summary <- group_summary[order(group_summary$dataset_id, group_summary$group_std, group_summary$matrix_format), ]

write.table(
  group_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_barcode_qc_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell barcode QC table written to res/qc/single_cell/single_cell_barcode_qc.tsv.gz\n")
cat("Single-cell barcode QC sample summary written to res/qc/single_cell/single_cell_barcode_qc_summary.tsv\n")
cat("Single-cell barcode QC group summary written to res/qc/single_cell/single_cell_barcode_qc_group_summary.tsv\n")
