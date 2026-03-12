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

read_tsv_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }

  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
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

collapse_to_gene_level <- function(expr, mapping) {
  sample_ids <- setdiff(colnames(expr), "feature_id")
  merged <- merge(expr, mapping, by = "feature_id", all.x = TRUE, sort = FALSE)
  merged <- merged[!(is.na(merged$gene_symbol) | merged$gene_symbol == ""), , drop = FALSE]

  if (nrow(merged) == 0) {
    stop("No mapped features retained after joining expression matrix to platform mapping.")
  }

  expr_only <- merged[, sample_ids, drop = FALSE]
  row_variance <- apply(expr_only, 1, var, na.rm = TRUE)
  merged$row_variance <- row_variance

  split_index <- split(seq_len(nrow(merged)), merged$gene_symbol)
  keep_index <- vapply(split_index, function(idx) {
    idx[which.max(merged$row_variance[idx])]
  }, integer(1))

  collapsed <- merged[keep_index, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]
  collapsed <- collapsed[order(collapsed$gene_symbol), , drop = FALSE]
  rownames(collapsed) <- NULL
  collapsed
}

build_summary_row <- function(dataset_id, gene_expr, original_expr) {
  data.frame(
    dataset_id = dataset_id,
    n_original_features = nrow(original_expr),
    n_gene_level_rows = nrow(gene_expr),
    n_samples = ncol(gene_expr) - 2,
    stringsAsFactors = FALSE
  )
}

mapping_dir <- project_path("data", "processed", "annotation")
bulk_dir <- project_path("data", "processed", "bulk")
bulk_gene_dir <- project_path("data", "processed", "bulk_gene")

datasets <- list(
  list(
    dataset_id = "GSE30528",
    expr_file = file.path(bulk_dir, "GSE30528_expr.tsv.gz"),
    mapping_file = file.path(mapping_dir, "GPL571_mapping.tsv.gz"),
    output_file = file.path(bulk_gene_dir, "GSE30528_gene_expr.tsv.gz")
  ),
  list(
    dataset_id = "GSE96804",
    expr_file = file.path(bulk_dir, "GSE96804_expr.tsv.gz"),
    mapping_file = file.path(mapping_dir, "GPL17586_mapping.tsv.gz"),
    output_file = file.path(bulk_gene_dir, "GSE96804_gene_expr.tsv.gz")
  ),
  list(
    dataset_id = "GSE104954_full",
    expr_file = file.path(bulk_dir, "GSE104954_full_expr.tsv.gz"),
    mapping_file = file.path(mapping_dir, "GSE104954_merged_mapping.tsv.gz"),
    output_file = file.path(bulk_gene_dir, "GSE104954_full_gene_expr.tsv.gz")
  ),
  list(
    dataset_id = "GSE104954_primary_dkd_vs_tumor",
    expr_file = file.path(bulk_dir, "GSE104954_primary_dkd_vs_tumor_expr.tsv.gz"),
    mapping_file = file.path(mapping_dir, "GSE104954_merged_mapping.tsv.gz"),
    output_file = file.path(bulk_gene_dir, "GSE104954_primary_dkd_vs_tumor_gene_expr.tsv.gz")
  )
)

summary_rows <- list()

for (dataset in datasets) {
  expr <- read_tsv_auto(dataset$expr_file)
  mapping <- read_tsv_auto(dataset$mapping_file)
  gene_expr <- collapse_to_gene_level(expr, mapping)
  write_tsv_gz(gene_expr, dataset$output_file)
  summary_rows[[length(summary_rows) + 1]] <- build_summary_row(dataset$dataset_id, gene_expr, expr)
}

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file = project_path("res", "qc", "bulk", "bulk_gene_level_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Gene-level bulk matrices written to data/processed/bulk_gene\n")
cat("Gene-level summary written to res/qc/bulk/bulk_gene_level_summary.tsv\n")
print(summary_table)
