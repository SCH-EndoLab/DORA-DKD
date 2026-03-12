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

suppressPackageStartupMessages({
  library(ggplot2)
  library(limma)
})

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

pheno <- read.delim(
  project_path("data", "processed", "bulk", "GSE189007_GPL23126_pheno.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE189007_GPL23126_gene_expr.tsv.gz"))

group_order <- c(
  "Control",
  "T2D_NoComp_LT5Y",
  "T2DN_LT5Y",
  "T2DR_LT5Y",
  "T2D_NoComp_GT15Y",
  "T2DN_GT15Y",
  "T2DR_GT15Y"
)

disease_axis_order <- c("Control", "T2D_NoComp", "T2DN", "T2DR")

pheno <- pheno[pheno$group_std %in% group_order, , drop = FALSE]
pheno$group_std <- factor(pheno$group_std, levels = group_order)
pheno$disease_axis_group <- factor(pheno$disease_axis_group, levels = disease_axis_order)

sample_ids <- pheno$sample_id
expr <- expr[, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]
expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
storage.mode(expr_matrix) <- "double"
rownames(expr_matrix) <- expr$gene_symbol

gene_sets <- read_tsv_auto(project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"))
gene_set_list <- split(gene_sets$gene_symbol, gene_sets$gene_set_id)

candidate_shortlist <- read.delim(
  project_path("res", "qc", "mechanism", "multiomics_candidate_shortlist.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
top_genes <- candidate_shortlist$gene_symbol[1:6]

score_rows <- list()
summary_rows <- list()

for (set_name in names(gene_set_list)) {
  score_obj <- score_gene_set_mean_z(expr_matrix, unique(gene_set_list[[set_name]]))
  score_df <- data.frame(
    sample_id = sample_ids,
    group_std = pheno$group_std[match(sample_ids, pheno$sample_id)],
    disease_axis_group = pheno$disease_axis_group[match(sample_ids, pheno$sample_id)],
    gene_set_id = set_name,
    score = as.numeric(score_obj$score),
    n_genes_present = length(score_obj$genes_present),
    stringsAsFactors = FALSE
  )
  score_rows[[length(score_rows) + 1]] <- score_df

  axis_means <- tapply(score_df$score, score_df$disease_axis_group, mean, na.rm = TRUE)
  axis_medians <- tapply(score_df$score, score_df$disease_axis_group, median, na.rm = TRUE)

  control_scores <- score_df$score[score_df$disease_axis_group == "Control"]
  t2d_scores <- score_df$score[score_df$disease_axis_group == "T2D_NoComp"]
  t2dn_scores <- score_df$score[score_df$disease_axis_group == "T2DN"]
  t2dr_scores <- score_df$score[score_df$disease_axis_group == "T2DR"]

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    gene_set_id = set_name,
    n_genes_present = length(score_obj$genes_present),
    mean_control = axis_means["Control"],
    mean_t2d_nocomp = axis_means["T2D_NoComp"],
    mean_t2dn = axis_means["T2DN"],
    mean_t2dr = axis_means["T2DR"],
    median_control = axis_medians["Control"],
    median_t2d_nocomp = axis_medians["T2D_NoComp"],
    median_t2dn = axis_medians["T2DN"],
    median_t2dr = axis_medians["T2DR"],
    p_t2dn_vs_t2d_nocomp = tryCatch(wilcox.test(t2dn_scores, t2d_scores)$p.value, error = function(e) NA_real_),
    p_t2dr_vs_t2d_nocomp = tryCatch(wilcox.test(t2dr_scores, t2d_scores)$p.value, error = function(e) NA_real_),
    p_control_vs_t2d_nocomp = tryCatch(wilcox.test(t2d_scores, control_scores)$p.value, error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}

score_table <- do.call(rbind, score_rows)
score_summary <- do.call(rbind, summary_rows)

write_tsv_gz(
  score_table,
  project_path("res", "tables", "mechanism", "GSE189007_GPL23126_mechanism_scores.tsv.gz")
)
write_tsv(
  score_summary,
  project_path("res", "qc", "mechanism", "GSE189007_GPL23126_mechanism_score_summary.tsv")
)

axis_keep <- pheno$disease_axis_group %in% disease_axis_order
axis_pheno <- pheno[axis_keep, , drop = FALSE]
axis_expr <- expr[, c("gene_symbol", "feature_id", axis_pheno$sample_id), drop = FALSE]

expr_long_rows <- list()
for (gene in intersect(top_genes, axis_expr$gene_symbol)) {
  gene_row <- axis_expr[axis_expr$gene_symbol == gene, , drop = FALSE][1, ]
  values <- as.numeric(gene_row[, axis_pheno$sample_id, drop = TRUE])
  expr_long_rows[[length(expr_long_rows) + 1]] <- data.frame(
    gene_symbol = gene,
    sample_id = axis_pheno$sample_id,
    group_std = axis_pheno$group_std,
    disease_axis_group = axis_pheno$disease_axis_group,
    expr_value = values,
    stringsAsFactors = FALSE
  )
}
expr_long <- do.call(rbind, expr_long_rows)

write_tsv_gz(
  expr_long,
  project_path("res", "tables", "mechanism", "GSE189007_GPL23126_top_candidate_expression.tsv.gz")
)

design <- model.matrix(~ 0 + disease_axis_group, data = axis_pheno)
colnames(design) <- levels(axis_pheno$disease_axis_group)
expr_limma <- as.matrix(axis_expr[, axis_pheno$sample_id, drop = FALSE])
rownames(expr_limma) <- make.unique(axis_expr$gene_symbol)

fit <- lmFit(expr_limma, design)
contrast_matrix <- makeContrasts(
  T2DN_vs_T2D_NoComp = T2DN - T2D_NoComp,
  T2DR_vs_T2D_NoComp = T2DR - T2D_NoComp,
  T2D_NoComp_vs_Control = T2D_NoComp - Control,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

limma_rows <- list()
for (coef_name in colnames(contrast_matrix)) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene_symbol <- axis_expr$gene_symbol[match(rownames(tt), rownames(expr_limma))]
  tt$feature_id <- axis_expr$feature_id[match(rownames(tt), rownames(expr_limma))]
  tt$contrast <- coef_name
  limma_rows[[coef_name]] <- tt[, c("contrast", "gene_symbol", "feature_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}
limma_table <- do.call(rbind, limma_rows)

write_tsv_gz(
  limma_table,
  project_path("res", "tables", "bulk", "GSE189007_GPL23126_diabetes_axis_limma.tsv.gz")
)

dir.create(project_path("figure", "export", "mechanism"), recursive = TRUE, showWarnings = FALSE)

score_table$gene_set_id <- factor(
  score_table$gene_set_id,
  levels = c("OXEIPTOSIS_CORE", "KEAP1_NRF2_RESPONSE", "OXEIPTOSIS_EXTENDED")
)

score_plot <- ggplot(
  score_table,
  aes(x = disease_axis_group, y = score, fill = disease_axis_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
  facet_wrap(~ gene_set_id, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Mean z-score",
    title = "Mechanism score shifts across the diabetes-complication axis"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "GSE189007_diabetes_axis_mechanism_scores.pdf"),
  plot = score_plot,
  width = 7.5,
  height = 9
)

expr_long$gene_symbol <- factor(expr_long$gene_symbol, levels = top_genes)

candidate_plot <- ggplot(
  expr_long,
  aes(x = disease_axis_group, y = expr_value, fill = disease_axis_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.7) +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Gene-level expression",
    title = "Top integrated candidates across the diabetes-complication axis"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "GSE189007_diabetes_axis_top_candidates.pdf"),
  plot = candidate_plot,
  width = 8.5,
  height = 8
)

cat("GSE189007 mechanism score table written to res/tables/mechanism/GSE189007_GPL23126_mechanism_scores.tsv.gz\n")
cat("GSE189007 mechanism score summary written to res/qc/mechanism/GSE189007_GPL23126_mechanism_score_summary.tsv\n")
cat("GSE189007 top candidate expression table written to res/tables/mechanism/GSE189007_GPL23126_top_candidate_expression.tsv.gz\n")
cat("GSE189007 diabetes-axis limma table written to res/tables/bulk/GSE189007_GPL23126_diabetes_axis_limma.tsv.gz\n")
cat("GSE189007 mechanism score plot written to figure/export/mechanism/GSE189007_diabetes_axis_mechanism_scores.pdf\n")
cat("GSE189007 top candidate plot written to figure/export/mechanism/GSE189007_diabetes_axis_top_candidates.pdf\n")
