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

celltype_order <- c(
  "PT",
  "PT_INJURED",
  "TAL",
  "DISTAL_CD",
  "ENDOTHELIAL",
  "MACROPHAGE",
  "MESANGIAL_PERICYTE",
  "PODOCYTE",
  "FIBROBLAST",
  "T_NK",
  "B_PLASMA",
  "PROLIFERATION",
  "UNRESOLVED"
)

celltype_palette <- c(
  PT = "#b35806",
  PT_INJURED = "#f1a340",
  TAL = "#998ec3",
  DISTAL_CD = "#542788",
  ENDOTHELIAL = "#5ab4ac",
  MACROPHAGE = "#d73027",
  MESANGIAL_PERICYTE = "#80cdc1",
  PODOCYTE = "#1b7837",
  FIBROBLAST = "#762a83",
  T_NK = "#4575b4",
  B_PLASMA = "#74add1",
  PROLIFERATION = "#f46d43",
  UNRESOLVED = "#bdbdbd"
)

composition_df <- read.delim(
  project_path("res", "tables", "single_cell", "single_cell_inferred_celltype_composition.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

group_df <- read.delim(
  project_path("res", "qc", "single_cell", "single_cell_typing_group_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_df <- read.delim(
  project_path("res", "qc", "single_cell", "single_cell_typing_sample_summary.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dir.create(project_path("figure", "export", "single_cell"), recursive = TRUE, showWarnings = FALSE)

composition_df$inferred_celltype <- factor(composition_df$inferred_celltype, levels = celltype_order)
group_df$inferred_celltype <- factor(group_df$inferred_celltype, levels = celltype_order)

sample_df$sample_label <- paste(sample_df$dataset_id, sample_df$group_std, sample_df$matrix_id, sep = " | ")
sample_order <- sample_df$sample_label[order(sample_df$dataset_id, sample_df$group_std, -sample_df$kept_cells, sample_df$matrix_id)]

composition_df$sample_label <- paste(composition_df$dataset_id, composition_df$group_std, composition_df$matrix_id, sep = " | ")
composition_df$sample_label <- factor(composition_df$sample_label, levels = unique(sample_order))
composition_df$dataset_group <- paste(composition_df$dataset_id, composition_df$group_std, sep = " | ")

sample_plot <- ggplot(
  composition_df,
  aes(x = sample_label, y = fraction_cells, fill = inferred_celltype)
) +
  geom_col(width = 0.88) +
  facet_grid(~ dataset_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = celltype_palette, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.background = element_rect(fill = "#f2f2f2", color = "#d0d0d0"),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "Fraction of filtered cells",
    fill = "Inferred cell type",
    title = "Sample-level coarse cell-type composition across kidney and urine cohorts"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_typing_sample_composition.pdf"),
  plot = sample_plot,
  width = 16,
  height = 6.5
)

group_df$dataset_group <- paste(group_df$dataset_id, group_df$group_std, sep = " | ")
group_order <- unique(group_df$dataset_group[order(group_df$dataset_id, group_df$group_std)])
group_df$dataset_group <- factor(group_df$dataset_group, levels = group_order)

fraction_plot <- ggplot(
  group_df,
  aes(x = inferred_celltype, y = dataset_group, fill = median_fraction_cells)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, max(group_df$median_fraction_cells, na.rm = TRUE))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Median fraction",
    title = "Median inferred cell-type fractions by dataset and group"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_typing_group_fraction_heatmap.pdf"),
  plot = fraction_plot,
  width = 9.5,
  height = 4.8
)

score_plot <- ggplot(
  group_df,
  aes(x = inferred_celltype, y = dataset_group, fill = median_oxeiptosis_extended_score)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient(low = "#fff5eb", high = "#7f2704") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Median score",
    title = "Median oxeiptosis-extended score by inferred cell type"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_typing_group_oxeiptosis_heatmap.pdf"),
  plot = score_plot,
  width = 9.5,
  height = 4.8
)

cat("Single-cell typing sample composition plot written to figure/export/single_cell/single_cell_typing_sample_composition.pdf\n")
cat("Single-cell typing group fraction heatmap written to figure/export/single_cell/single_cell_typing_group_fraction_heatmap.pdf\n")
cat("Single-cell typing group oxeiptosis heatmap written to figure/export/single_cell/single_cell_typing_group_oxeiptosis_heatmap.pdf\n")
