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

read_geo_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin_idx <- which(lines == "!series_matrix_table_begin")
  end_idx <- which(lines == "!series_matrix_table_end")

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop("Could not find a valid series matrix table in: ", path)
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
  colnames(mat)[1] <- "cpg_id"
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

expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE233758_series_matrix.txt.gz"))
pheno <- read.delim(
  project_path("data", "metadata", "soft_samples", "GSE233758_samples.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pheno$group_std <- ifelse(
  grepl("diabetic nephropathy", pheno$title, ignore.case = TRUE) | grepl("_TD", pheno$title, ignore.case = TRUE),
  "DKD",
  "Control"
)
pheno$analysis_label <- "proximal_tubule_methylation"

sample_ids <- pheno$sample_id
expr <- expr[, c("cpg_id", sample_ids), drop = FALSE]

write_tsv_gz(
  expr,
  project_path("data", "processed", "methylation_GSE233758_beta.tsv.gz")
)
write_tsv(
  pheno,
  project_path("data", "processed", "GSE233758_pheno.tsv")
)

summary_table <- data.frame(
  dataset_id = "GSE233758",
  n_cpg = nrow(expr),
  n_samples = ncol(expr) - 1,
  n_control = sum(pheno$group_std == "Control"),
  n_dkd = sum(pheno$group_std == "DKD"),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  project_path("res", "qc", "mechanism", "GSE233758_methylation_prep_summary.tsv")
)

cat("Processed GSE233758 beta matrix written to data/processed/methylation_GSE233758_beta.tsv.gz\n")
cat("Processed GSE233758 phenotype written to data/processed/GSE233758_pheno.tsv\n")
print(summary_table)
