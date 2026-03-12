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

typing_df <- read.delim(
  project_path("res", "tables", "single_cell", "single_cell_inferred_celltypes.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

composition_df <- read.delim(
  project_path("res", "tables", "single_cell", "single_cell_inferred_celltype_composition.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dir.create(project_path("res", "qc", "single_cell"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figure", "export", "single_cell"), recursive = TRUE, showWarnings = FALSE)

sample_key <- paste(typing_df$dataset_id, typing_df$sample_id, typing_df$matrix_id, typing_df$group_std, sep = "||")

sample_score_summary <- do.call(
  rbind,
  lapply(split(typing_df, paste(sample_key, typing_df$inferred_celltype, sep = "||")), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      sample_id = df$sample_id[1],
      sample_title = df$sample_title[1],
      group_std = df$group_std[1],
      matrix_id = df$matrix_id[1],
      inferred_celltype = df$inferred_celltype[1],
      n_cells = nrow(df),
      median_oxeiptosis_extended_score = median(df$oxeiptosis_extended_score, na.rm = TRUE),
      mean_oxeiptosis_extended_score = mean(df$oxeiptosis_extended_score, na.rm = TRUE),
      median_nCount = median(df$nCount, na.rm = TRUE),
      median_nFeature = median(df$nFeature, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

sample_score_summary$dataset_group <- paste(sample_score_summary$dataset_id, sample_score_summary$group_std, sep = " | ")

group_score_summary <- do.call(
  rbind,
  lapply(split(sample_score_summary, paste(sample_score_summary$dataset_id, sample_score_summary$group_std, sample_score_summary$inferred_celltype, sep = "||")), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      group_std = df$group_std[1],
      inferred_celltype = df$inferred_celltype[1],
      n_samples = nrow(df),
      median_sample_cells = median(df$n_cells),
      median_sample_oxeiptosis_extended_score = median(df$median_oxeiptosis_extended_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

tubular_axis_df <- do.call(
  rbind,
  lapply(split(typing_df, sample_key), function(df) {
    tubular_idx <- df$inferred_celltype %in% c("PT", "PT_INJURED")
    data.frame(
      dataset_id = df$dataset_id[1],
      sample_id = df$sample_id[1],
      sample_title = df$sample_title[1],
      group_std = df$group_std[1],
      matrix_id = df$matrix_id[1],
      total_cells = nrow(df),
      tubular_axis_cells = sum(tubular_idx),
      tubular_axis_fraction = mean(tubular_idx),
      tubular_axis_median_oxeiptosis_extended_score = median(df$oxeiptosis_extended_score[tubular_idx], na.rm = TRUE),
      pt_fraction = mean(df$inferred_celltype == "PT"),
      pt_injured_fraction = mean(df$inferred_celltype == "PT_INJURED"),
      stringsAsFactors = FALSE
    )
  })
)

tubular_axis_df$dataset_group <- paste(tubular_axis_df$dataset_id, tubular_axis_df$group_std, sep = " | ")

tubular_group_summary <- do.call(
  rbind,
  lapply(split(tubular_axis_df, tubular_axis_df$dataset_group), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      group_std = df$group_std[1],
      n_samples = nrow(df),
      median_tubular_axis_fraction = median(df$tubular_axis_fraction, na.rm = TRUE),
      median_tubular_axis_score = median(df$tubular_axis_median_oxeiptosis_extended_score, na.rm = TRUE),
      median_pt_fraction = median(df$pt_fraction, na.rm = TRUE),
      median_pt_injured_fraction = median(df$pt_injured_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  sample_score_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_celltype_score_sample_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  group_score_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_celltype_score_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  tubular_axis_df,
  file = project_path("res", "qc", "single_cell", "single_cell_tubular_axis_sample_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  tubular_group_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_tubular_axis_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

dataset_group_order <- c(
  "GSE209781 | Control",
  "GSE209781 | DKD",
  "GSE279086 | Control",
  "GSE279086 | Diabetes",
  "GSE266146 | DKD_urine_pool"
)

tubular_axis_df$dataset_group <- factor(tubular_axis_df$dataset_group, levels = dataset_group_order)

fraction_plot <- ggplot(
  tubular_axis_df,
  aes(x = dataset_group, y = tubular_axis_fraction, fill = dataset_group)
) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.14, size = 2.1, alpha = 0.85) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Fraction of PT + PT_INJURED cells",
    title = "Tubular injury-axis abundance at the sample level"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_tubular_axis_fraction_boxplot.pdf"),
  plot = fraction_plot,
  width = 8.8,
  height = 4.8
)

score_plot <- ggplot(
  tubular_axis_df,
  aes(x = dataset_group, y = tubular_axis_median_oxeiptosis_extended_score, fill = dataset_group)
) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.14, size = 2.1, alpha = 0.85) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Median oxeiptosis-extended score in PT/PT_INJURED",
    title = "Tubular injury-axis mechanism score at the sample level"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_tubular_axis_score_boxplot.pdf"),
  plot = score_plot,
  width = 8.8,
  height = 4.8
)

focus_df <- subset(
  sample_score_summary,
  inferred_celltype %in% c("PT", "PT_INJURED", "TAL", "ENDOTHELIAL") & n_cells >= 20
)
focus_df$dataset_group <- factor(focus_df$dataset_group, levels = dataset_group_order)
focus_df$inferred_celltype <- factor(focus_df$inferred_celltype, levels = c("PT", "PT_INJURED", "TAL", "ENDOTHELIAL"))

focus_plot <- ggplot(
  focus_df,
  aes(x = dataset_group, y = median_oxeiptosis_extended_score, color = inferred_celltype)
) +
  geom_point(
    aes(size = n_cells),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65),
    alpha = 0.85
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Sample-level median oxeiptosis-extended score",
    color = "Cell type",
    size = "Cells",
    title = "Cell-type-resolved mechanism score shifts across groups"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_celltype_score_focus_dotplot.pdf"),
  plot = focus_plot,
  width = 10.2,
  height = 5.2
)

cat("Single-cell celltype score sample summary written to res/qc/single_cell/single_cell_celltype_score_sample_summary.tsv\n")
cat("Single-cell celltype score group summary written to res/qc/single_cell/single_cell_celltype_score_group_summary.tsv\n")
cat("Single-cell tubular axis sample summary written to res/qc/single_cell/single_cell_tubular_axis_sample_summary.tsv\n")
cat("Single-cell tubular axis group summary written to res/qc/single_cell/single_cell_tubular_axis_group_summary.tsv\n")
cat("Single-cell tubular axis fraction plot written to figure/export/single_cell/single_cell_tubular_axis_fraction_boxplot.pdf\n")
cat("Single-cell tubular axis score plot written to figure/export/single_cell/single_cell_tubular_axis_score_boxplot.pdf\n")
cat("Single-cell celltype score focus plot written to figure/export/single_cell/single_cell_celltype_score_focus_dotplot.pdf\n")
