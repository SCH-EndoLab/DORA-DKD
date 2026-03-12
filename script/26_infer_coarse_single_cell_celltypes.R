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
      if (length(x) >= 2) x[2] else x[1]
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

safe_pct <- function(num, den) {
  out <- rep(0, length(den))
  keep <- den > 0
  out[keep] <- 100 * num[keep] / den[keep]
  out
}

compute_gene_set_scores <- function(mat, features, lib_size, marker_list) {
  gene_upper <- toupper(features)
  score_list <- lapply(names(marker_list), function(set_id) {
    marker_genes <- unique(toupper(marker_list[[set_id]]))
    idx <- which(gene_upper %in% marker_genes)
    if (length(idx) == 0) {
      return(rep(NA_real_, ncol(mat)))
    }
    submat <- mat[idx, , drop = FALSE]
    norm_submat <- sweep(submat, 2, lib_size, "/") * 10000
    as.numeric(Matrix::colMeans(log1p(norm_submat)))
  })
  score_df <- as.data.frame(score_list, stringsAsFactors = FALSE)
  names(score_df) <- names(marker_list)
  score_df
}

assign_celltype <- function(score_df) {
  score_mat <- as.matrix(score_df)
  top_idx <- max.col(score_mat, ties.method = "first")
  top_score <- score_mat[cbind(seq_len(nrow(score_mat)), top_idx)]
  celltype <- colnames(score_mat)[top_idx]

  second_score <- rep(NA_real_, nrow(score_mat))
  for (i in seq_len(nrow(score_mat))) {
    x <- score_mat[i, ]
    ord <- order(x, decreasing = TRUE, na.last = TRUE)
    second_score[i] <- if (length(ord) >= 2) x[ord[2]] else NA_real_
  }

  margin <- top_score - second_score
  celltype[is.na(top_score) | top_score < 0.35 | margin < 0.05] <- "UNRESOLVED"

  data.frame(
    inferred_celltype = celltype,
    top_score = top_score,
    second_score = second_score,
    score_margin = margin,
    stringsAsFactors = FALSE
  )
}

marker_df <- read.delim(
  project_path("data", "metadata", "gene_sets", "kidney_celltype_marker_membership.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mechanism_df <- read.delim(
  project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

marker_list <- split(marker_df$gene_symbol, marker_df$marker_set_id)
mechanism_extended <- unique(mechanism_df$gene_symbol[mechanism_df$gene_set_id == "OXEIPTOSIS_EXTENDED"])

manifest <- read.delim(
  project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Use primary kidney raw matrices, processed external kidney matrices, and urine matrices.
# Exclude GSE279086 raw fallback matrices from the primary typing layer.
typing_manifest <- subset(
  manifest,
  dataset_id %in% c("GSE209781", "GSE266146", "GSE279086") &
    matrix_format != "10x_mtx_raw_fallback"
)

barcode_out <- list()
sample_out <- list()
celltype_out <- list()

for (i in seq_len(nrow(typing_manifest))) {
  row_df <- typing_manifest[i, ]
  features <- read_feature_symbols(row_df$features_path)
  barcodes <- read_lines(row_df$barcodes_path)
  mat <- read_sparse_mtx(row_df$matrix_path)

  ncount <- Matrix::colSums(mat)
  nfeature <- Matrix::colSums(mat > 0)
  mito_idx <- grepl("^MT-", toupper(features))
  mito_counts <- if (any(mito_idx)) Matrix::colSums(mat[mito_idx, , drop = FALSE]) else rep(0, ncol(mat))
  pct_mito <- safe_pct(as.numeric(mito_counts), as.numeric(ncount))

  if (row_df$dataset_id == "GSE209781") {
    keep <- ncount >= 500 & nfeature >= 200 & pct_mito <= 20
  } else {
    keep <- ncount >= 500 & nfeature >= 200 & pct_mito <= 20
  }

  kept_mat <- mat[, keep, drop = FALSE]
  kept_barcodes <- barcodes[keep]
  kept_ncount <- as.numeric(ncount[keep])
  kept_nfeature <- as.numeric(nfeature[keep])
  kept_pct_mito <- pct_mito[keep]

  if (ncol(kept_mat) == 0) {
    sample_out[[i]] <- data.frame(
      dataset_id = row_df$dataset_id,
      sample_id = row_df$sample_id,
      sample_title = row_df$sample_title,
      group_std = row_df$group_std,
      matrix_id = row_df$matrix_id,
      matrix_format = row_df$matrix_format,
      analysis_role = row_df$analysis_role,
      total_barcodes = length(barcodes),
      kept_cells = 0,
      median_nCount = NA_real_,
      median_nFeature = NA_real_,
      median_pct_mito = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  score_df <- compute_gene_set_scores(
    mat = kept_mat,
    features = features,
    lib_size = kept_ncount,
    marker_list = marker_list
  )
  assign_df <- assign_celltype(score_df)

  mech_scores <- compute_gene_set_scores(
    mat = kept_mat,
    features = features,
    lib_size = kept_ncount,
    marker_list = list(OXEIPTOSIS_EXTENDED = mechanism_extended)
  )

  barcode_df <- data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    barcode = kept_barcodes,
    nCount = kept_ncount,
    nFeature = kept_nfeature,
    pct_mito = kept_pct_mito,
    oxeiptosis_extended_score = mech_scores$OXEIPTOSIS_EXTENDED,
    stringsAsFactors = FALSE
  )
  barcode_df <- cbind(barcode_df, assign_df, score_df)
  barcode_out[[i]] <- barcode_df

  sample_out[[i]] <- data.frame(
    dataset_id = row_df$dataset_id,
    sample_id = row_df$sample_id,
    sample_title = row_df$sample_title,
    group_std = row_df$group_std,
    matrix_id = row_df$matrix_id,
    matrix_format = row_df$matrix_format,
    analysis_role = row_df$analysis_role,
    total_barcodes = length(barcodes),
    kept_cells = ncol(kept_mat),
    median_nCount = median(kept_ncount),
    median_nFeature = median(kept_nfeature),
    median_pct_mito = median(kept_pct_mito),
    stringsAsFactors = FALSE
  )

  ct_tab <- as.data.frame(table(barcode_df$inferred_celltype), stringsAsFactors = FALSE)
  names(ct_tab) <- c("inferred_celltype", "n_cells")
  ct_tab$dataset_id <- row_df$dataset_id
  ct_tab$sample_id <- row_df$sample_id
  ct_tab$sample_title <- row_df$sample_title
  ct_tab$group_std <- row_df$group_std
  ct_tab$matrix_id <- row_df$matrix_id
  ct_tab$matrix_format <- row_df$matrix_format
  ct_tab$analysis_role <- row_df$analysis_role
  ct_tab$fraction_cells <- ct_tab$n_cells / ncol(kept_mat)
  ct_tab$median_oxeiptosis_extended_score <- vapply(
    ct_tab$inferred_celltype,
    function(ct) median(barcode_df$oxeiptosis_extended_score[barcode_df$inferred_celltype == ct], na.rm = TRUE),
    numeric(1)
  )
  celltype_out[[i]] <- ct_tab

  rm(mat, kept_mat, score_df, assign_df, barcode_df, ct_tab, mech_scores)
  gc(verbose = FALSE)
}

barcode_res <- do.call(rbind, barcode_out)
sample_res <- do.call(rbind, sample_out)
celltype_res <- do.call(rbind, celltype_out)

dir.create(project_path("res", "tables", "single_cell"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)

write.table(
  barcode_res,
  file = gzfile(project_path("res", "tables", "single_cell", "single_cell_inferred_celltypes.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sample_res,
  file = project_path("res", "qc", "single_cell", "single_cell_typing_sample_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  celltype_res,
  file = gzfile(project_path("res", "tables", "single_cell", "single_cell_inferred_celltype_composition.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

group_key <- paste(celltype_res$dataset_id, celltype_res$group_std, celltype_res$inferred_celltype, sep = "||")
group_summary <- do.call(
  rbind,
  lapply(split(celltype_res, group_key), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      group_std = df$group_std[1],
      inferred_celltype = df$inferred_celltype[1],
      n_samples = nrow(df),
      total_cells = sum(df$n_cells),
      median_fraction_cells = median(df$fraction_cells),
      median_oxeiptosis_extended_score = median(df$median_oxeiptosis_extended_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  group_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_typing_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell inferred celltypes written to res/tables/single_cell/single_cell_inferred_celltypes.tsv.gz\n")
cat("Single-cell typing sample summary written to res/qc/single_cell/single_cell_typing_sample_summary.tsv\n")
cat("Single-cell typing group summary written to res/qc/single_cell/single_cell_typing_group_summary.tsv\n")
