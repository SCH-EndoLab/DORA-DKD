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

suppressPackageStartupMessages(library(ggplot2))

read_tsv_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }

  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

scores <- read_tsv_auto(project_path("res", "tables", "mechanism", "bulk_mechanism_sample_scores.tsv.gz"))
score_summary <- read_tsv_auto(project_path("res", "qc", "mechanism", "bulk_mechanism_score_summary.tsv"))
core_expr <- read_tsv_auto(project_path("res", "tables", "mechanism", "bulk_oxeiptosis_core_gene_expression.tsv.gz"))

scores$dataset_id <- factor(
  scores$dataset_id,
  levels = c("GSE30528", "GSE96804", "GSE104954_primary_dkd_vs_tumor")
)
scores$gene_set_id <- factor(
  scores$gene_set_id,
  levels = c("OXEIPTOSIS_CORE", "KEAP1_NRF2_RESPONSE", "OXEIPTOSIS_EXTENDED")
)
scores$group_std <- factor(scores$group_std, levels = c("Control", "DKD"))

score_labels <- unique(score_summary[, c("dataset_id", "gene_set_id", "delta_mean_dkd_minus_control", "wilcox_fdr")])
score_labels$label <- paste0(
  "delta=", sprintf("%.2f", score_labels$delta_mean_dkd_minus_control),
  "\nFDR=", sprintf("%.2g", score_labels$wilcox_fdr)
)

score_plot <- ggplot(scores, aes(x = group_std, y = score, fill = group_std)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.8) +
  facet_grid(gene_set_id ~ dataset_id, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#7aa6c2", DKD = "#d56f5d")) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f1efe7"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Mean z-score", title = "Mechanism gene-set scores across kidney bulk cohorts")

ggsave(
  filename = project_path("figure", "export", "mechanism", "bulk_mechanism_scores_boxplot.pdf"),
  plot = score_plot,
  width = 10,
  height = 7
)

core_expr$dataset_id <- factor(
  core_expr$dataset_id,
  levels = c("GSE30528", "GSE96804", "GSE104954_primary_dkd_vs_tumor")
)
core_expr$group_std <- factor(core_expr$group_std, levels = c("Control", "DKD"))

core_plot <- ggplot(core_expr, aes(x = group_std, y = expr_value, fill = group_std)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.75) +
  facet_grid(gene_symbol ~ dataset_id, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#7aa6c2", DKD = "#d56f5d")) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f1efe7"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Gene-level expression", title = "Core oxeiptosis gene expression across kidney bulk cohorts")

ggsave(
  filename = project_path("figure", "export", "mechanism", "bulk_oxeiptosis_core_gene_expression.pdf"),
  plot = core_plot,
  width = 10,
  height = 7
)

cat("Mechanism score plot written to figure/export/mechanism/bulk_mechanism_scores_boxplot.pdf\n")
cat("Core gene expression plot written to figure/export/mechanism/bulk_oxeiptosis_core_gene_expression.pdf\n")
