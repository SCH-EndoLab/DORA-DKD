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

fast_spearman_gene_score <- function(expr_mat, score_vector) {
  n <- length(score_vector)
  score_rank <- rank(score_vector, ties.method = "average")
  expr_rank <- t(apply(expr_mat, 1, rank, ties.method = "average"))
  rho <- cor(t(expr_rank), score_rank, method = "pearson")
  rho <- as.numeric(rho)
  names(rho) <- rownames(expr_mat)

  denom <- pmax(1 - rho^2, 1e-12)
  t_stat <- rho * sqrt((n - 2) / denom)
  p_value <- 2 * pt(-abs(t_stat), df = n - 2)

  data.frame(
    gene_symbol = rownames(expr_mat),
    rho = rho,
    p_value = p_value,
    fdr = p.adjust(p_value, method = "BH"),
    stringsAsFactors = FALSE
  )
}

consensus <- read_tsv_auto(project_path("res", "tables", "bulk", "bulk_strong_glomerular_consensus.tsv.gz"))
sample_scores <- read_tsv_auto(project_path("res", "tables", "mechanism", "bulk_mechanism_sample_scores.tsv.gz"))
gene_sets <- read_tsv_auto(project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"))

target_genes <- unique(consensus$gene_symbol)
mechanism_members <- unique(gene_sets$gene_symbol[gene_sets$gene_set_id == "OXEIPTOSIS_EXTENDED"])

dataset_specs <- list(
  list(
    dataset_id = "GSE30528",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE30528_gene_expr.tsv.gz")
  ),
  list(
    dataset_id = "GSE96804",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE96804_gene_expr.tsv.gz")
  ),
  list(
    dataset_id = "GSE104954_primary_dkd_vs_tumor",
    expr_file = project_path("data", "processed", "bulk_gene", "GSE104954_primary_dkd_vs_tumor_gene_expr.tsv.gz")
  )
)

cor_rows <- list()

for (spec in dataset_specs) {
  expr <- read_tsv_auto(spec$expr_file)
  score_df <- sample_scores[
    sample_scores$dataset_id == spec$dataset_id &
      sample_scores$gene_set_id == "OXEIPTOSIS_EXTENDED",
    ,
    drop = FALSE
  ]

  sample_ids <- score_df$sample_id
  expr <- expr[expr$gene_symbol %in% target_genes, c("gene_symbol", sample_ids), drop = FALSE]
  expr_mat <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_mat) <- "double"
  rownames(expr_mat) <- expr$gene_symbol

  cor_df <- fast_spearman_gene_score(expr_mat, score_df$score)
  cor_df$dataset_id <- spec$dataset_id
  cor_rows[[length(cor_rows) + 1]] <- cor_df
}

cor_all <- do.call(rbind, cor_rows)

wide_merge <- Reduce(
  function(a, b) merge(a, b, by = "gene_symbol", all = TRUE),
  lapply(split(cor_all, cor_all$dataset_id), function(df) {
    dataset_id <- unique(df$dataset_id)
    df <- df[, c("gene_symbol", "rho", "p_value", "fdr")]
    names(df)[-1] <- paste0(names(df)[-1], "_", dataset_id)
    df
  })
)

ranked <- merge(consensus, wide_merge, by = "gene_symbol", all.x = TRUE)
ranked$is_mechanism_member <- ranked$gene_symbol %in% mechanism_members
ranked$same_sign_score_corr_glomerular <- sign(ranked$rho_GSE30528) == sign(ranked$rho_GSE96804)
ranked$same_sign_score_corr_all_available <- with(
  ranked,
  sign(rho_GSE30528) == sign(rho_GSE96804) &
    sign(rho_GSE30528) == sign(rho_GSE104954_primary_dkd_vs_tumor)
)

ranked$score_link_tier <- ifelse(
  ranked$same_sign_score_corr_glomerular &
    ranked$fdr_GSE30528 < 0.05 &
    ranked$fdr_GSE96804 < 0.05,
  "glomerular_fdr",
  ifelse(
    ranked$same_sign_score_corr_glomerular &
      ranked$p_value_GSE30528 < 0.05 &
      ranked$p_value_GSE96804 < 0.05,
    "glomerular_nominal",
    ifelse(ranked$same_sign_score_corr_glomerular, "same_direction_only", "discordant")
  )
)

ranked$integrated_priority_score <- (-log10(pmax(ranked$fisher_fdr, 1e-300))) *
  rowMeans(cbind(abs(ranked$rho_GSE30528), abs(ranked$rho_GSE96804)), na.rm = TRUE) *
  ifelse(ranked$tubule_support == "same_direction_fdr", 1.5,
         ifelse(ranked$tubule_support == "same_direction_nominal", 1.2, 1.0))

ranked <- ranked[order(
  factor(ranked$score_link_tier, levels = c("glomerular_fdr", "glomerular_nominal", "same_direction_only", "discordant")),
  -ranked$integrated_priority_score
), ]

top_ranked <- ranked[, c(
  "gene_symbol", "fisher_fdr", "tubule_support", "is_mechanism_member",
  "rho_GSE30528", "fdr_GSE30528", "rho_GSE96804", "fdr_GSE96804",
  "rho_GSE104954_primary_dkd_vs_tumor", "fdr_GSE104954_primary_dkd_vs_tumor",
  "score_link_tier", "integrated_priority_score"
)]

core_limma_rows <- list()
for (dataset_id in c("GSE30528", "GSE96804", "GSE104954_primary_dkd_vs_tumor")) {
  limma_df <- read_tsv_auto(project_path("res", "tables", "bulk", paste0(dataset_id, "_limma_dkd_vs_control.tsv.gz")))
  core_limma_rows[[length(core_limma_rows) + 1]] <- limma_df[limma_df$gene_symbol %in% mechanism_members, ]
}
mechanism_gene_limma <- do.call(rbind, core_limma_rows)

summary_table <- data.frame(
  metric = c(
    "consensus_genes_total",
    "score_link_glomerular_fdr",
    "score_link_glomerular_nominal",
    "score_link_same_direction_only",
    "mechanism_members_within_consensus"
  ),
  value = c(
    nrow(top_ranked),
    sum(top_ranked$score_link_tier == "glomerular_fdr", na.rm = TRUE),
    sum(top_ranked$score_link_tier == "glomerular_nominal", na.rm = TRUE),
    sum(top_ranked$score_link_tier == "same_direction_only", na.rm = TRUE),
    sum(top_ranked$is_mechanism_member, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write_tsv_gz(
  top_ranked,
  project_path("res", "tables", "mechanism", "mechanism_linked_candidate_ranking.tsv.gz")
)
write_tsv_gz(
  mechanism_gene_limma,
  project_path("res", "tables", "mechanism", "mechanism_member_limma_across_bulk.tsv.gz")
)
write_tsv(
  summary_table,
  project_path("res", "qc", "mechanism", "mechanism_candidate_summary.tsv")
)

cat("Mechanism-linked candidate ranking written to res/tables/mechanism/mechanism_linked_candidate_ranking.tsv.gz\n")
cat("Mechanism member limma table written to res/tables/mechanism/mechanism_member_limma_across_bulk.tsv.gz\n")
cat("Mechanism candidate summary written to res/qc/mechanism/mechanism_candidate_summary.tsv\n")
print(summary_table)
