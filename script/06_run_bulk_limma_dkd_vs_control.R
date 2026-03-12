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

suppressPackageStartupMessages(library(limma))

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

run_limma_case_control <- function(dataset_id, expr_file, pheno_file) {
  expr <- read_tsv_auto(expr_file)
  pheno <- read_tsv_auto(pheno_file)
  pheno <- pheno[pheno$group_std %in% c("Control", "DKD"), , drop = FALSE]

  sample_ids <- pheno$sample_id
  expr <- expr[, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]

  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- make.unique(expr$gene_symbol)

  group <- factor(pheno$group_std, levels = c("Control", "DKD"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(expr_matrix, design)
  contrast_matrix <- makeContrasts(DKD - Control, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

  result <- topTable(fit2, number = Inf, sort.by = "P")
  result$gene_symbol <- expr$gene_symbol[match(rownames(result), rownames(expr_matrix))]
  result$feature_id <- expr$feature_id[match(rownames(result), rownames(expr_matrix))]
  result$dataset_id <- dataset_id
  result$sign_direction <- ifelse(result$logFC > 0, "Up_in_DKD", "Down_in_DKD")
  result$sig_fdr <- ifelse(result$adj.P.Val < 0.05, "yes", "no")
  result$sig_fdr_fc <- ifelse(result$adj.P.Val < 0.05 & abs(result$logFC) >= 0.5, "yes", "no")
  result <- result[, c(
    "dataset_id", "gene_symbol", "feature_id", "logFC", "AveExpr",
    "t", "P.Value", "adj.P.Val", "B", "sign_direction", "sig_fdr", "sig_fdr_fc"
  )]

  result
}

dataset_specs <- list(
  list(
    dataset_id = "GSE30528",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE30528_gene_expr.tsv.gz"),
    pheno_file = project_path("data", "processed", "bulk", "GSE30528_pheno.tsv")
  ),
  list(
    dataset_id = "GSE96804",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE96804_gene_expr.tsv.gz"),
    pheno_file = project_path("data", "processed", "bulk", "GSE96804_pheno.tsv")
  ),
  list(
    dataset_id = "GSE104954_primary_dkd_vs_tumor",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE104954_primary_dkd_vs_tumor_gene_expr.tsv.gz"),
    pheno_file = project_path("data", "processed", "bulk", "GSE104954_primary_dkd_vs_tumor_pheno.tsv")
  )
)

result_list <- list()
summary_rows <- list()

for (spec in dataset_specs) {
  result <- run_limma_case_control(spec$dataset_id, spec$expr_file, spec$pheno_file)
  result_list[[spec$dataset_id]] <- result

  out_file <- project_path(
    "res", "tables", "bulk",
    paste0(spec$dataset_id, "_limma_dkd_vs_control.tsv.gz")
  )
  write_tsv_gz(result, out_file)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    dataset_id = spec$dataset_id,
    n_tested = nrow(result),
    n_sig_fdr = sum(result$sig_fdr == "yes"),
    n_sig_fdr_fc = sum(result$sig_fdr_fc == "yes"),
    top_gene = result$gene_symbol[[1]],
    top_logFC = result$logFC[[1]],
    top_adj_p = result$adj.P.Val[[1]],
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file = project_path("res", "qc", "bulk", "bulk_limma_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

sig_direction_lists <- lapply(result_list, function(x) {
  subset(x, adj.P.Val < 0.05, select = c("gene_symbol", "logFC", "sign_direction"))
})

common_genes <- Reduce(intersect, lapply(sig_direction_lists, function(x) unique(x$gene_symbol)))
common_rows <- lapply(names(sig_direction_lists), function(dataset_id) {
  x <- sig_direction_lists[[dataset_id]]
  x <- x[x$gene_symbol %in% common_genes, , drop = FALSE]
  names(x)[2:3] <- paste0(names(x)[2:3], "_", dataset_id)
  x
})

common_overlap <- Reduce(function(a, b) merge(a, b, by = "gene_symbol", all = FALSE), common_rows)

direction_cols <- grep("^sign_direction_", names(common_overlap), value = TRUE)
common_overlap$consistent_direction <- apply(
  common_overlap[, direction_cols, drop = FALSE],
  1,
  function(x) length(unique(x)) == 1
)

common_overlap <- common_overlap[common_overlap$consistent_direction, , drop = FALSE]
write_tsv_gz(
  common_overlap,
  project_path("res", "tables", "bulk", "bulk_limma_common_consistent_genes.tsv.gz")
)

cat("Bulk limma results written to res/tables/bulk\n")
cat("Limma summary written to res/qc/bulk/bulk_limma_summary.tsv\n")
cat("Consistent overlap written to res/tables/bulk/bulk_limma_common_consistent_genes.tsv.gz\n")
print(summary_table)
cat("\nConsistent overlap genes:", nrow(common_overlap), "\n")
