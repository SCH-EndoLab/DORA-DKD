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

candidate_ranking <- read.delim(
  gzfile(project_path("res", "tables", "mechanism", "mechanism_linked_candidate_ranking.tsv.gz")),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

methylation_df <- read.delim(
  project_path("res", "qc", "mechanism", "GSE233758_top_candidate_gene_methylation_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

single_cell_group <- read.delim(
  project_path("res", "qc", "single_cell", "single_cell_candidate_gene_group_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_genes <- sort(unique(single_cell_group$gene_symbol))

get_sc_value <- function(dataset_id, group_std, inferred_celltype, gene_symbol, field) {
  idx <- single_cell_group$dataset_id == dataset_id &
    single_cell_group$group_std == group_std &
    single_cell_group$inferred_celltype == inferred_celltype &
    single_cell_group$gene_symbol == gene_symbol
  vals <- single_cell_group[[field]][idx]
  if (length(vals) == 0) NA_real_ else vals[1]
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

rows <- lapply(candidate_genes, function(gene) {
  bulk_row <- candidate_ranking[candidate_ranking$gene_symbol == gene, , drop = FALSE]
  meth_row <- methylation_df[methylation_df$gene_symbol == gene, , drop = FALSE]

  gse209_pt_ctrl <- get_sc_value("GSE209781", "Control", "PT", gene, "median_mean_log_expr")
  gse209_pt_dkd <- get_sc_value("GSE209781", "DKD", "PT", gene, "median_mean_log_expr")
  gse209_pti_ctrl <- get_sc_value("GSE209781", "Control", "PT_INJURED", gene, "median_mean_log_expr")
  gse209_pti_dkd <- get_sc_value("GSE209781", "DKD", "PT_INJURED", gene, "median_mean_log_expr")

  gse209_tal_dkd <- get_sc_value("GSE209781", "DKD", "TAL", gene, "median_mean_log_expr")
  gse209_endo_dkd <- get_sc_value("GSE209781", "DKD", "ENDOTHELIAL", gene, "median_mean_log_expr")
  gse209_macro_dkd <- get_sc_value("GSE209781", "DKD", "MACROPHAGE", gene, "median_mean_log_expr")

  gse279_pt_ctrl <- get_sc_value("GSE279086", "Control", "PT", gene, "median_mean_log_expr")
  gse279_pt_dm <- get_sc_value("GSE279086", "Diabetes", "PT", gene, "median_mean_log_expr")
  gse279_pti_ctrl <- get_sc_value("GSE279086", "Control", "PT_INJURED", gene, "median_mean_log_expr")
  gse279_pti_dm <- get_sc_value("GSE279086", "Diabetes", "PT_INJURED", gene, "median_mean_log_expr")

  urine_pti <- get_sc_value("GSE266146", "DKD_urine_pool", "PT_INJURED", gene, "median_mean_log_expr")

  tubular_mean_dkd <- safe_mean(c(gse209_pt_dkd, gse209_pti_dkd))
  non_tubular_mean_dkd <- safe_mean(c(gse209_tal_dkd, gse209_endo_dkd, gse209_macro_dkd))
  tubular_specificity_dkd <- tubular_mean_dkd - non_tubular_mean_dkd

  pt_delta <- gse209_pt_dkd - gse209_pt_ctrl
  pti_delta <- gse209_pti_dkd - gse209_pti_ctrl
  external_tubular_delta <- safe_mean(c(gse279_pt_dm - gse279_pt_ctrl, gse279_pti_dm - gse279_pti_ctrl))

  methylation_supported <- nrow(meth_row) > 0 && any(is.finite(meth_row$best_adj_p))
  methylation_near_sig <- methylation_supported && min(meth_row$best_adj_p, na.rm = TRUE) <= 0.1

  support_points <- 0
  support_points <- support_points + ifelse(nrow(bulk_row) == 1 && bulk_row$integrated_priority_score >= 5, 1, 0)
  support_points <- support_points + ifelse(methylation_near_sig, 1, 0)
  support_points <- support_points + ifelse(is.finite(tubular_specificity_dkd) && tubular_specificity_dkd >= 0.2, 1, 0)
  support_points <- support_points + ifelse(is.finite(pt_delta) && pt_delta > 0.1, 1, 0)
  support_points <- support_points + ifelse(is.finite(pti_delta) && pti_delta > 0.1, 1, 0)
  support_points <- support_points + ifelse(is.finite(external_tubular_delta) && external_tubular_delta > 0, 1, 0)
  support_points <- support_points + ifelse(is.finite(urine_pti) && urine_pti > 0.1, 1, 0)

  data.frame(
    gene_symbol = gene,
    bulk_integrated_priority_score = if (nrow(bulk_row) == 1) bulk_row$integrated_priority_score else NA_real_,
    bulk_fisher_fdr = if (nrow(bulk_row) == 1) bulk_row$fisher_fdr else NA_real_,
    bulk_tubule_support = if (nrow(bulk_row) == 1) bulk_row$tubule_support else NA_character_,
    methylation_best_adj_p = if (nrow(meth_row) > 0) min(meth_row$best_adj_p, na.rm = TRUE) else NA_real_,
    methylation_best_delta_beta = if (nrow(meth_row) > 0) meth_row$best_delta_beta[which.min(meth_row$best_adj_p)][1] else NA_real_,
    methylation_direction = if (nrow(meth_row) > 0) meth_row$best_dm_direction[which.min(meth_row$best_adj_p)][1] else NA_character_,
    gse209781_pt_control = gse209_pt_ctrl,
    gse209781_pt_dkd = gse209_pt_dkd,
    gse209781_pt_delta = pt_delta,
    gse209781_pt_injured_control = gse209_pti_ctrl,
    gse209781_pt_injured_dkd = gse209_pti_dkd,
    gse209781_pt_injured_delta = pti_delta,
    gse209781_tubular_specificity_dkd = tubular_specificity_dkd,
    gse279086_tubular_external_delta = external_tubular_delta,
    gse266146_pt_injured_expr = urine_pti,
    support_points = support_points,
    stringsAsFactors = FALSE
  )
})

shortlist_df <- do.call(rbind, rows)
shortlist_df <- shortlist_df[order(-shortlist_df$support_points, -shortlist_df$bulk_integrated_priority_score), ]

write.table(
  shortlist_df,
  file = project_path("res", "qc", "mechanism", "multiomics_candidate_shortlist.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Multi-omics candidate shortlist written to res/qc/mechanism/multiomics_candidate_shortlist.tsv\n")
