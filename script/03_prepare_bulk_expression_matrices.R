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

read_sample_sheet <- function(dataset_id) {
  path <- project_path("data", "metadata", "soft_samples", paste0(dataset_id, "_samples.tsv"))
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

read_geo_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin_idx <- which(lines == "!series_matrix_table_begin")
  end_idx <- which(lines == "!series_matrix_table_end")

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop("Could not find a valid series_matrix table block in: ", path)
  }

  block <- lines[(begin_idx + 1):(end_idx - 1)]
  handle <- textConnection(block)
  on.exit(close(handle), add = TRUE)

  mat <- read.delim(
    handle,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )

  colnames(mat)[1] <- "feature_id"
  mat
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

reorder_matrix <- function(expr, sample_ids) {
  missing_ids <- setdiff(sample_ids, colnames(expr))
  if (length(missing_ids) > 0) {
    stop("Missing samples in matrix: ", paste(missing_ids, collapse = ", "))
  }

  expr[, c("feature_id", sample_ids), drop = FALSE]
}

merge_split_series_matrix <- function(paths) {
  mats <- lapply(paths, read_geo_series_matrix)
  feature_ids <- mats[[1]]$feature_id

  for (mat in mats[-1]) {
    if (!identical(feature_ids, mat$feature_id)) {
      stop("Feature IDs are not aligned across split matrices.")
    }
  }

  merged <- mats[[1]]

  for (mat in mats[-1]) {
    merged <- cbind(merged, mat[, -1, drop = FALSE])
  }

  merged
}

build_summary_row <- function(dataset_id, output_label, expr, pheno) {
  primary_samples <- sum(pheno$keep_for_primary_contrast == "yes", na.rm = TRUE)

  data.frame(
    dataset_id = dataset_id,
    output_label = output_label,
    n_features = nrow(expr),
    n_samples = ncol(expr) - 1,
    n_primary_samples = primary_samples,
    primary_groups = paste(sort(unique(pheno$group_std[pheno$keep_for_primary_contrast == "yes"])), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

bulk_dir <- project_path("data", "processed", "bulk")
summary_rows <- list()

gse30528_expr <- read_geo_series_matrix(
  project_path("data", "raw", "series_matrix", "GSE30528_series_matrix.txt.gz")
)
gse30528_pheno <- read_sample_sheet("GSE30528")
gse30528_pheno$group_std <- ifelse(
  gse30528_pheno$char_disease_state == "diabetic kidney disease (DKD)",
  "DKD",
  ifelse(gse30528_pheno$char_disease_state == "control", "Control", NA_character_)
)
gse30528_pheno$keep_for_primary_contrast <- ifelse(
  gse30528_pheno$group_std %in% c("DKD", "Control"),
  "yes",
  "no"
)
gse30528_pheno$analysis_label <- "glomerular_discovery"
gse30528_expr <- reorder_matrix(gse30528_expr, gse30528_pheno$sample_id)

write_tsv_gz(gse30528_expr, file.path(bulk_dir, "GSE30528_expr.tsv.gz"))
write_tsv(gse30528_pheno, file.path(bulk_dir, "GSE30528_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE30528", "full", gse30528_expr, gse30528_pheno)

gse96804_expr <- read_geo_series_matrix(
  project_path("data", "raw", "series_matrix", "GSE96804_series_matrix.txt.gz")
)
gse96804_pheno <- read_sample_sheet("GSE96804")
gse96804_pheno$group_std <- ifelse(
  grepl("diabetic nephropathy", gse96804_pheno$source_name_ch1, ignore.case = TRUE),
  "DKD",
  ifelse(
    grepl("tumor nephrectom", gse96804_pheno$source_name_ch1, ignore.case = TRUE),
    "Control",
    NA_character_
  )
)
gse96804_pheno$keep_for_primary_contrast <- ifelse(
  gse96804_pheno$group_std %in% c("DKD", "Control"),
  "yes",
  "no"
)
gse96804_pheno$analysis_label <- "glomerular_validation"
gse96804_expr <- reorder_matrix(gse96804_expr, gse96804_pheno$sample_id)

write_tsv_gz(gse96804_expr, file.path(bulk_dir, "GSE96804_expr.tsv.gz"))
write_tsv(gse96804_pheno, file.path(bulk_dir, "GSE96804_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE96804", "full", gse96804_expr, gse96804_pheno)

gse104954_expr <- merge_split_series_matrix(
  c(
    project_path("data", "raw", "series_matrix", "GSE104954-GPL24120_series_matrix.txt.gz"),
    project_path("data", "raw", "series_matrix", "GSE104954-GPL22945_series_matrix.txt.gz")
  )
)
gse104954_pheno <- read_sample_sheet("GSE104954")
gse104954_pheno$group_std <- ifelse(
  gse104954_pheno$char_diagnosis == "Diabetic nephropathy",
  "DKD",
  ifelse(gse104954_pheno$char_diagnosis == "Tumor nephrectomy", "Control", "Other")
)
gse104954_pheno$keep_for_primary_contrast <- ifelse(
  gse104954_pheno$group_std %in% c("DKD", "Control"),
  "yes",
  "no"
)
gse104954_pheno$analysis_label <- "tubulointerstitial_validation"
gse104954_expr <- reorder_matrix(gse104954_expr, gse104954_pheno$sample_id)

write_tsv_gz(gse104954_expr, file.path(bulk_dir, "GSE104954_full_expr.tsv.gz"))
write_tsv(gse104954_pheno, file.path(bulk_dir, "GSE104954_full_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE104954", "full", gse104954_expr, gse104954_pheno)

gse104954_primary_pheno <- gse104954_pheno[gse104954_pheno$keep_for_primary_contrast == "yes", , drop = FALSE]
gse104954_primary_expr <- reorder_matrix(gse104954_expr, gse104954_primary_pheno$sample_id)

write_tsv_gz(
  gse104954_primary_expr,
  file.path(bulk_dir, "GSE104954_primary_dkd_vs_tumor_expr.tsv.gz")
)
write_tsv(
  gse104954_primary_pheno,
  file.path(bulk_dir, "GSE104954_primary_dkd_vs_tumor_pheno.tsv")
)
summary_rows[[length(summary_rows) + 1]] <- build_summary_row(
  "GSE104954",
  "primary_dkd_vs_tumor",
  gse104954_primary_expr,
  gse104954_primary_pheno
)

summary_table <- do.call(rbind, summary_rows)

write_tsv(
  summary_table,
  project_path("res", "qc", "bulk", "bulk_matrix_prep_summary.tsv")
)

cat("Bulk expression matrices written to data/processed/bulk\n")
cat("Preparation summary written to res/qc/bulk/bulk_matrix_prep_summary.tsv\n")
print(summary_table)
