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

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))

read_lines <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(readLines(gzfile(path), warn = FALSE))
  }
  readLines(path, warn = FALSE)
}

read_feature_symbols <- function(path) {
  lines <- read_lines(path)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  gene_symbol <- vapply(
    split_lines,
    function(x) {
      if (length(x) >= 2) x[2] else x[1]
    },
    character(1)
  )
  trimws(gene_symbol)
}

read_sparse_mtx <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readMM(con)
}

candidate_genes <- c(
  "MACROD1",
  "HSPA1A",
  "VTI1B",
  "TALDO1",
  "APOO",
  "VPS28",
  "ETFB",
  "MRPL16",
  "DCXR",
  "SLC9A3R1",
  "ACOT7",
  "EPHX2"
)

focus_celltypes <- c("PT", "PT_INJURED", "TAL", "DISTAL_CD", "ENDOTHELIAL", "MACROPHAGE")

typing_df <- read.delim(
  project_path("res", "tables", "single_cell", "single_cell_inferred_celltypes.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest <- read.delim(
  project_path("data", "metadata", "single_cell_matrix_manifest.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest <- subset(manifest, matrix_id %in% unique(typing_df$matrix_id))

sample_gene_out <- list()
presence_out <- list()

for (i in seq_len(nrow(manifest))) {
  row_df <- manifest[i, ]
  sample_typing <- subset(typing_df, matrix_id == row_df$matrix_id)
  if (nrow(sample_typing) == 0) {
    next
  }

  features <- read_feature_symbols(row_df$features_path)
  barcodes <- read_lines(row_df$barcodes_path)
  mat <- read_sparse_mtx(row_df$matrix_path)
  feature_upper <- toupper(features)

  barcode_idx <- match(sample_typing$barcode, barcodes)
  keep <- !is.na(barcode_idx)
  if (!any(keep)) {
    next
  }

  sample_typing <- sample_typing[keep, , drop = FALSE]
  barcode_idx <- barcode_idx[keep]

  present_genes <- candidate_genes[candidate_genes %in% feature_upper]
  missing_genes <- setdiff(candidate_genes, present_genes)
  if (length(missing_genes) > 0) {
    presence_out[[length(presence_out) + 1]] <- data.frame(
      dataset_id = row_df$dataset_id,
      matrix_id = row_df$matrix_id,
      gene_symbol = missing_genes,
      present_in_matrix = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (length(present_genes) == 0) {
    next
  }

  presence_out[[length(presence_out) + 1]] <- data.frame(
    dataset_id = row_df$dataset_id,
    matrix_id = row_df$matrix_id,
    gene_symbol = present_genes,
    present_in_matrix = TRUE,
    stringsAsFactors = FALSE
  )

  for (gene in present_genes) {
    gene_idx <- which(feature_upper == gene)
    gene_counts <- if (length(gene_idx) == 1) {
      as.numeric(mat[gene_idx, barcode_idx, drop = TRUE])
    } else {
      as.numeric(Matrix::colSums(mat[gene_idx, barcode_idx, drop = FALSE]))
    }
    log_expr <- log1p(gene_counts / sample_typing$nCount * 10000)

    by_celltype <- split(seq_len(nrow(sample_typing)), sample_typing$inferred_celltype)
    for (celltype in names(by_celltype)) {
      idx <- by_celltype[[celltype]]
      sample_gene_out[[length(sample_gene_out) + 1]] <- data.frame(
        dataset_id = row_df$dataset_id,
        sample_id = row_df$sample_id,
        sample_title = row_df$sample_title,
        group_std = row_df$group_std,
        matrix_id = row_df$matrix_id,
        gene_symbol = gene,
        inferred_celltype = celltype,
        n_cells = length(idx),
        mean_log_expr = mean(log_expr[idx], na.rm = TRUE),
        median_log_expr = median(log_expr[idx], na.rm = TRUE),
        detection_fraction = mean(gene_counts[idx] > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  rm(mat, sample_typing, features, barcodes)
  gc(verbose = FALSE)
}

sample_gene_summary <- do.call(rbind, sample_gene_out)
presence_summary <- do.call(rbind, presence_out)

group_gene_summary <- do.call(
  rbind,
  lapply(split(sample_gene_summary, paste(sample_gene_summary$dataset_id, sample_gene_summary$group_std, sample_gene_summary$inferred_celltype, sample_gene_summary$gene_symbol, sep = "||")), function(df) {
    data.frame(
      dataset_id = df$dataset_id[1],
      group_std = df$group_std[1],
      inferred_celltype = df$inferred_celltype[1],
      gene_symbol = df$gene_symbol[1],
      n_samples = nrow(df),
      median_sample_cells = median(df$n_cells),
      median_mean_log_expr = median(df$mean_log_expr, na.rm = TRUE),
      median_detection_fraction = median(df$detection_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  sample_gene_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_candidate_gene_sample_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  group_gene_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_candidate_gene_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  presence_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_candidate_gene_presence_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_df <- subset(
  group_gene_summary,
  inferred_celltype %in% focus_celltypes & n_samples >= 2
)
plot_df$dataset_group <- paste(plot_df$dataset_id, plot_df$group_std, sep = " | ")
plot_df$dataset_group <- factor(
  plot_df$dataset_group,
  levels = c(
    "GSE209781 | Control",
    "GSE209781 | DKD",
    "GSE279086 | Control",
    "GSE279086 | Diabetes",
    "GSE266146 | DKD_urine_pool"
  )
)
plot_df$inferred_celltype <- factor(plot_df$inferred_celltype, levels = focus_celltypes)
plot_df$gene_symbol <- factor(plot_df$gene_symbol, levels = candidate_genes)

dir.create(project_path("figure", "export", "single_cell"), recursive = TRUE, showWarnings = FALSE)

expr_plot <- ggplot(
  plot_df,
  aes(x = gene_symbol, y = dataset_group, fill = median_mean_log_expr)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(inferred_celltype ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "#f2f2f2", color = "#d0d0d0")
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Median mean log expr",
    title = "Candidate-gene localization across inferred kidney cell types"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_candidate_gene_expression_heatmap.pdf"),
  plot = expr_plot,
  width = 9.8,
  height = 9.2
)

detection_plot <- ggplot(
  plot_df,
  aes(x = gene_symbol, y = dataset_group, fill = median_detection_fraction)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(inferred_celltype ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#fff5eb", high = "#7f2704", limits = c(0, 1)) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "#f2f2f2", color = "#d0d0d0")
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Median detection",
    title = "Detection frequency of candidate genes across inferred cell types"
  )

ggsave(
  filename = project_path("figure", "export", "single_cell", "single_cell_candidate_gene_detection_heatmap.pdf"),
  plot = detection_plot,
  width = 9.8,
  height = 9.2
)

cat("Single-cell candidate gene sample summary written to res/qc/single_cell/single_cell_candidate_gene_sample_summary.tsv\n")
cat("Single-cell candidate gene group summary written to res/qc/single_cell/single_cell_candidate_gene_group_summary.tsv\n")
cat("Single-cell candidate gene presence summary written to res/qc/single_cell/single_cell_candidate_gene_presence_summary.tsv\n")
cat("Single-cell candidate gene expression heatmap written to figure/export/single_cell/single_cell_candidate_gene_expression_heatmap.pdf\n")
cat("Single-cell candidate gene detection heatmap written to figure/export/single_cell/single_cell_candidate_gene_detection_heatmap.pdf\n")
