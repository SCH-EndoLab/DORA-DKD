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

sample_qc <- read.delim(
  project_path("res", "qc", "single_cell", "single_cell_barcode_qc_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

group_qc <- read.delim(
  project_path("res", "qc", "single_cell", "single_cell_barcode_qc_group_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dir.create(project_path("figure", "export", "single_cell"), recursive = TRUE, showWarnings = FALSE)

sample_qc$plot_group <- paste(sample_qc$dataset_id, sample_qc$group_std, sep = " | ")
sample_qc$matrix_class <- sample_qc$matrix_format

scatter_plot <- ggplot(
  sample_qc,
  aes(x = median_nFeature, y = median_pct_mito, color = plot_group, shape = matrix_class, size = n_cells)
) +
  geom_point(alpha = 0.85) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#8f8f8f") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Median detected features per barcode",
    y = "Median mitochondrial percentage",
    color = "Dataset | group",
    shape = "Matrix class",
    size = "Nonzero cells",
    title = "Single-cell matrix QC landscape across kidney and urine cohorts"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_qc_feature_mito_scatter.pdf"),
  plot = scatter_plot,
  width = 9.5,
  height = 6.5
)

group_qc$plot_group <- paste(group_qc$dataset_id, group_qc$group_std, sep = " | ")
ordered_levels <- unique(group_qc$plot_group[order(group_qc$median_pct_mito)])
group_qc$plot_group <- factor(group_qc$plot_group, levels = ordered_levels)

mito_plot <- ggplot(
  group_qc,
  aes(x = plot_group, y = median_pct_mito, fill = matrix_format)
) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#8f8f8f") +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Median mitochondrial percentage",
    fill = "Matrix class",
    title = "Group-level mitochondrial burden across single-cell cohorts"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_group_mito_barplot.pdf"),
  plot = mito_plot,
  width = 8.5,
  height = 5.5
)

filter_plot <- ggplot(
  group_qc,
  aes(x = plot_group, y = cells_core_pass_200_500_20 / total_cells, fill = matrix_format)
) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Fraction passing core barcode filter",
    fill = "Matrix class",
    title = "Core barcode-filter retention across single-cell cohorts"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_group_core_filter_fraction.pdf"),
  plot = filter_plot,
  width = 8.5,
  height = 5.5
)

cat("Single-cell QC scatter written to figure/export/single_cell/single_cell_qc_feature_mito_scatter.pdf\n")
cat("Single-cell group mito barplot written to figure/export/single_cell/single_cell_group_mito_barplot.pdf\n")
cat("Single-cell filter fraction barplot written to figure/export/single_cell/single_cell_group_core_filter_fraction.pdf\n")
