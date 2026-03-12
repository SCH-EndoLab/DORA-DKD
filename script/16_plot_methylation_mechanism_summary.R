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

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

clean_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

figure_dir <- project_path("figure", "export", "mechanism")
ensure_dir(figure_dir)

methylation_limma <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE233758_methylation_limma.tsv.gz"))
mechanism_gene_cpgs <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE233758_mechanism_gene_cpgs.tsv.gz"))
mechanism_summary <- read_tsv_auto(project_path("res", "qc", "mechanism", "GSE233758_mechanism_gene_methylation_summary.tsv"))
candidate_summary <- read_tsv_auto(project_path("res", "qc", "mechanism", "GSE233758_top_candidate_gene_methylation_summary.tsv"))

gene_symbol_col <- if ("primary_gene_symbol" %in% colnames(methylation_limma)) {
  "primary_gene_symbol"
} else if ("gene_symbol_primary" %in% colnames(methylation_limma)) {
  "gene_symbol_primary"
} else {
  stop("No recognized primary gene symbol column found in methylation limma table.")
}

delta_beta_col <- if ("delta_beta" %in% colnames(methylation_limma)) {
  "delta_beta"
} else if ("delta_beta_dkd_minus_control" %in% colnames(methylation_limma)) {
  "delta_beta_dkd_minus_control"
} else {
  stop("No recognized delta beta column found in methylation limma table.")
}

methylation_limma$primary_gene_symbol <- clean_gene_symbol(methylation_limma[[gene_symbol_col]])
methylation_limma$delta_beta_plot <- methylation_limma[[delta_beta_col]]
methylation_limma$label_group <- "Background"

mechanism_cpg_ids <- unique(mechanism_gene_cpgs$cpg_id)
methylation_limma$label_group[methylation_limma$cpg_id %in% mechanism_cpg_ids] <- "Mechanism-linked CpG"

top_labeled <- candidate_summary[order(candidate_summary$best_adj_p, -abs(candidate_summary$best_delta_beta)), ]
top_labeled <- head(top_labeled, 12)
top_labeled$gene_symbol <- factor(top_labeled$gene_symbol, levels = rev(top_labeled$gene_symbol))

candidate_gene_set <- unique(as.character(top_labeled$gene_symbol))
methylation_limma$label_group[
  !is.na(methylation_limma$primary_gene_symbol) &
    methylation_limma$primary_gene_symbol %in% candidate_gene_set &
    methylation_limma$label_group == "Background"
] <- "Top candidate CpG"

methylation_limma$neg_log10_adj_p <- -log10(pmax(methylation_limma$adj.P.Val, 1e-300))
methylation_limma$label_group <- factor(
  methylation_limma$label_group,
  levels = c("Background", "Mechanism-linked CpG", "Top candidate CpG")
)

volcano_labels <- methylation_limma[
  methylation_limma$label_group != "Background" &
    !is.na(methylation_limma$primary_gene_symbol),
  c("cpg_id", "primary_gene_symbol", "delta_beta_plot", "neg_log10_adj_p", "label_group")
]
volcano_labels <- volcano_labels[order(-volcano_labels$neg_log10_adj_p, -abs(volcano_labels$delta_beta_plot)), ]
volcano_labels <- unique(volcano_labels[, c("primary_gene_symbol", "cpg_id", "delta_beta_plot", "neg_log10_adj_p", "label_group")])
volcano_labels <- head(volcano_labels, 18)
volcano_labels$label <- paste0(volcano_labels$primary_gene_symbol, "\n", volcano_labels$cpg_id)

volcano_plot <- ggplot(methylation_limma, aes(x = delta_beta_plot, y = neg_log10_adj_p, color = label_group)) +
  geom_point(size = 0.9, alpha = 0.65) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = "dashed", color = "#8f8f8f") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#8f8f8f") +
  geom_text(
    data = volcano_labels,
    aes(label = label),
    size = 2.6,
    check_overlap = TRUE,
    color = "#202020",
    vjust = -0.2,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Background" = "#c8c8c8",
      "Mechanism-linked CpG" = "#d56f5d",
      "Top candidate CpG" = "#4d8db7"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = "Delta beta (DKD - control)",
    y = expression(-log[10]("adjusted p")),
    title = "GSE233758 methylation landscape with mechanism-linked highlights"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "GSE233758_methylation_volcano.pdf"),
  plot = volcano_plot,
  width = 9.5,
  height = 6.5
)

mechanism_plot_df <- mechanism_summary[order(mechanism_summary$best_adj_p, -abs(mechanism_summary$best_delta_beta)), ]
mechanism_plot_df$gene_symbol <- factor(mechanism_plot_df$gene_symbol, levels = rev(mechanism_plot_df$gene_symbol))
mechanism_plot_df$neg_log10_adj_p <- -log10(pmax(mechanism_plot_df$best_adj_p, 1e-300))

mechanism_plot <- ggplot(
  mechanism_plot_df,
  aes(x = best_delta_beta, y = gene_symbol, size = n_cpg, color = neg_log10_adj_p)
) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8f8f8f") +
  scale_color_gradient(low = "#9ecae1", high = "#cb3e38") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Best delta beta (DKD - control)",
    y = NULL,
    color = expression(-log[10]("adj p")),
    size = "CpG count",
    title = "Mechanism gene methylation support in microdissected proximal tubules"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "GSE233758_mechanism_gene_methylation_summary.pdf"),
  plot = mechanism_plot,
  width = 8.8,
  height = 5.8
)

candidate_plot_df <- candidate_summary[order(candidate_summary$best_adj_p, -abs(candidate_summary$best_delta_beta)), ]
candidate_plot_df <- head(candidate_plot_df, 20)
candidate_plot_df$gene_symbol <- factor(candidate_plot_df$gene_symbol, levels = rev(candidate_plot_df$gene_symbol))
candidate_plot_df$neg_log10_adj_p <- -log10(pmax(candidate_plot_df$best_adj_p, 1e-300))

candidate_plot <- ggplot(
  candidate_plot_df,
  aes(x = best_delta_beta, y = gene_symbol, size = n_cpg, color = neg_log10_adj_p)
) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8f8f8f") +
  scale_color_gradient(low = "#c7dcef", high = "#2c6da4") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = "Best delta beta (DKD - control)",
    y = NULL,
    color = expression(-log[10]("adj p")),
    size = "CpG count",
    title = "Top candidate genes with methylation support in GSE233758"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "GSE233758_top_candidate_methylation_summary.pdf"),
  plot = candidate_plot,
  width = 9,
  height = 6.8
)

cat("Methylation volcano plot written to figure/export/mechanism/GSE233758_methylation_volcano.pdf\n")
cat("Mechanism methylation summary plot written to figure/export/mechanism/GSE233758_mechanism_gene_methylation_summary.pdf\n")
cat("Top candidate methylation summary plot written to figure/export/mechanism/GSE233758_top_candidate_methylation_summary.pdf\n")
