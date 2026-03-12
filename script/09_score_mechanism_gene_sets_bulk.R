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

score_gene_set_mean_z <- function(expr_matrix, genes) {
  genes_present <- intersect(genes, rownames(expr_matrix))
  if (length(genes_present) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_matrix)), genes_present = character(0)))
  }

  sub <- expr_matrix[genes_present, , drop = FALSE]
  row_means <- rowMeans(sub, na.rm = TRUE)
  row_sds <- apply(sub, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  z <- sweep(sub, 1, row_means, "-")
  z <- sweep(z, 1, row_sds, "/")

  list(score = colMeans(z, na.rm = TRUE), genes_present = genes_present)
}

bulk_specs <- list(
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

gene_sets <- read_tsv_auto(project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"))
gene_set_list <- split(gene_sets$gene_symbol, gene_sets$gene_set_id)

sample_score_rows <- list()
score_summary_rows <- list()
gene_presence_rows <- list()
core_gene_expr_rows <- list()

for (spec in bulk_specs) {
  expr <- read_tsv_auto(spec$expr_file)
  pheno <- read_tsv_auto(spec$pheno_file)
  pheno <- pheno[pheno$group_std %in% c("Control", "DKD"), , drop = FALSE]

  sample_ids <- pheno$sample_id
  expr <- expr[, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]
  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- expr$gene_symbol

  for (set_name in names(gene_set_list)) {
    score_obj <- score_gene_set_mean_z(expr_matrix, unique(gene_set_list[[set_name]]))
    score_df <- data.frame(
      dataset_id = spec$dataset_id,
      sample_id = sample_ids,
      group_std = pheno$group_std[match(sample_ids, pheno$sample_id)],
      gene_set_id = set_name,
      score = as.numeric(score_obj$score),
      n_genes_present = length(score_obj$genes_present),
      genes_present = paste(score_obj$genes_present, collapse = ";"),
      stringsAsFactors = FALSE
    )
    sample_score_rows[[length(sample_score_rows) + 1]] <- score_df

    control_scores <- score_df$score[score_df$group_std == "Control"]
    dkd_scores <- score_df$score[score_df$group_std == "DKD"]

    score_summary_rows[[length(score_summary_rows) + 1]] <- data.frame(
      dataset_id = spec$dataset_id,
      gene_set_id = set_name,
      n_control = sum(score_df$group_std == "Control"),
      n_dkd = sum(score_df$group_std == "DKD"),
      n_genes_present = length(score_obj$genes_present),
      mean_control = mean(control_scores, na.rm = TRUE),
      mean_dkd = mean(dkd_scores, na.rm = TRUE),
      delta_mean_dkd_minus_control = mean(dkd_scores, na.rm = TRUE) - mean(control_scores, na.rm = TRUE),
      median_control = median(control_scores, na.rm = TRUE),
      median_dkd = median(dkd_scores, na.rm = TRUE),
      wilcox_p = tryCatch(wilcox.test(dkd_scores, control_scores)$p.value, error = function(e) NA_real_),
      ttest_p = tryCatch(t.test(dkd_scores, control_scores)$p.value, error = function(e) NA_real_),
      stringsAsFactors = FALSE
    )

    gene_presence_rows[[length(gene_presence_rows) + 1]] <- data.frame(
      dataset_id = spec$dataset_id,
      gene_set_id = set_name,
      gene_symbol = unique(gene_set_list[[set_name]]),
      present_in_dataset = unique(gene_set_list[[set_name]]) %in% rownames(expr_matrix),
      stringsAsFactors = FALSE
    )
  }

  core_genes <- unique(gene_sets$gene_symbol[gene_sets$gene_set_id == "OXEIPTOSIS_CORE"])
  core_genes_present <- intersect(core_genes, rownames(expr_matrix))
  if (length(core_genes_present) > 0) {
    core_df <- data.frame(
      dataset_id = spec$dataset_id,
      gene_symbol = rep(core_genes_present, each = length(sample_ids)),
      sample_id = rep(sample_ids, times = length(core_genes_present)),
      group_std = rep(pheno$group_std[match(sample_ids, pheno$sample_id)], times = length(core_genes_present)),
      expr_value = as.vector(t(expr_matrix[core_genes_present, sample_ids, drop = FALSE])),
      stringsAsFactors = FALSE
    )
    core_gene_expr_rows[[length(core_gene_expr_rows) + 1]] <- core_df
  }
}

sample_scores <- do.call(rbind, sample_score_rows)
score_summary <- do.call(rbind, score_summary_rows)
score_summary$wilcox_fdr <- p.adjust(score_summary$wilcox_p, method = "BH")
score_summary$ttest_fdr <- p.adjust(score_summary$ttest_p, method = "BH")
gene_presence <- do.call(rbind, gene_presence_rows)
core_gene_expr <- do.call(rbind, core_gene_expr_rows)

write_tsv_gz(
  sample_scores,
  project_path("res", "tables", "mechanism", "bulk_mechanism_sample_scores.tsv.gz")
)
write_tsv(
  score_summary,
  project_path("res", "qc", "mechanism", "bulk_mechanism_score_summary.tsv")
)
write_tsv(
  gene_presence,
  project_path("res", "qc", "mechanism", "bulk_mechanism_gene_presence.tsv")
)
write_tsv_gz(
  core_gene_expr,
  project_path("res", "tables", "mechanism", "bulk_oxeiptosis_core_gene_expression.tsv.gz")
)

cat("Mechanism sample scores written to res/tables/mechanism/bulk_mechanism_sample_scores.tsv.gz\n")
cat("Mechanism score summary written to res/qc/mechanism/bulk_mechanism_score_summary.tsv\n")
cat("Mechanism gene presence written to res/qc/mechanism/bulk_mechanism_gene_presence.tsv\n")
cat("Core gene expression table written to res/tables/mechanism/bulk_oxeiptosis_core_gene_expression.tsv.gz\n")
print(score_summary)
