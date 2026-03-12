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
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

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

write_tsvgz <- function(x, path) {
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

parse_ratio <- function(x) {
  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(
    parts,
    function(y) {
      if (length(y) != 2L) {
        return(NA_real_)
      }
      as.numeric(y[1]) / as.numeric(y[2])
    },
    numeric(1)
  )
}

wrap_text <- function(x, width = 48) {
  vapply(
    x,
    function(y) paste(strwrap(y, width = width), collapse = "\n"),
    character(1)
  )
}

bind_rows_fill <- function(df_list) {
  df_list <- df_list[!vapply(df_list, is.null, logical(1))]
  if (!length(df_list)) {
    return(NULL)
  }
  all_cols <- unique(unlist(lapply(df_list, colnames), use.names = FALSE))
  df_list <- lapply(
    df_list,
    function(df) {
      missing_cols <- setdiff(all_cols, colnames(df))
      if (length(missing_cols)) {
        for (col_name in missing_cols) {
          df[[col_name]] <- NA
        }
      }
      df[, all_cols, drop = FALSE]
    }
  )
  do.call(rbind, df_list)
}

map_symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE))
  }
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  mapped <- mapped[!is.na(mapped$SYMBOL) & !is.na(mapped$ENTREZID), c("SYMBOL", "ENTREZID"), drop = FALSE]
  mapped <- unique(mapped)
  mapped[order(mapped$SYMBOL, mapped$ENTREZID), , drop = FALSE]
}

as_result_table <- function(enrich_obj, input_id, database, ontology, n_input_genes, n_mapped_symbols, n_mapped_entrez) {
  if (is.null(enrich_obj)) {
    return(NULL)
  }
  result_tbl <- as.data.frame(enrich_obj)
  if (!nrow(result_tbl)) {
    return(NULL)
  }
  result_tbl$input_id <- input_id
  result_tbl$database <- database
  result_tbl$ontology <- ontology
  result_tbl$n_input_genes <- n_input_genes
  result_tbl$n_mapped_symbols <- n_mapped_symbols
  result_tbl$n_mapped_entrez <- n_mapped_entrez
  result_tbl$gene_ratio_numeric <- parse_ratio(result_tbl$GeneRatio)
  result_tbl$bg_ratio_numeric <- parse_ratio(result_tbl$BgRatio)
  result_tbl
}

run_one_enrichment <- function(gene_symbols, universe_entrez, input_id) {
  gene_map <- map_symbol_to_entrez(gene_symbols)
  gene_entrez <- unique(gene_map$ENTREZID)

  input_summary <- data.frame(
    input_id = input_id,
    n_input_genes = length(unique(gene_symbols)),
    n_mapped_symbols = length(unique(gene_map$SYMBOL)),
    n_mapped_entrez = length(gene_entrez),
    stringsAsFactors = FALSE
  )

  if (!length(gene_entrez)) {
    return(list(summary = input_summary, results = NULL))
  }

  go_bp <- tryCatch(
    enrichGO(
      gene = gene_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  go_mf <- tryCatch(
    enrichGO(
      gene = gene_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  go_cc <- tryCatch(
    enrichGO(
      gene = gene_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  kegg <- tryCatch(
    enrichKEGG(
      gene = gene_entrez,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) NULL
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg))) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  reactome <- tryCatch(
    enrichPathway(
      gene = gene_entrez,
      universe = universe_entrez,
      organism = "human",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  result_tbl <- bind_rows_fill(
    list(
      as_result_table(go_bp, input_id, "GO", "BP", input_summary$n_input_genes, input_summary$n_mapped_symbols, input_summary$n_mapped_entrez),
      as_result_table(go_mf, input_id, "GO", "MF", input_summary$n_input_genes, input_summary$n_mapped_symbols, input_summary$n_mapped_entrez),
      as_result_table(go_cc, input_id, "GO", "CC", input_summary$n_input_genes, input_summary$n_mapped_symbols, input_summary$n_mapped_entrez),
      as_result_table(kegg, input_id, "KEGG", NA_character_, input_summary$n_input_genes, input_summary$n_mapped_symbols, input_summary$n_mapped_entrez),
      as_result_table(reactome, input_id, "Reactome", NA_character_, input_summary$n_input_genes, input_summary$n_mapped_symbols, input_summary$n_mapped_entrez)
    )
  )

  list(summary = input_summary, results = result_tbl)
}

make_dotplot <- function(result_tbl, focus_database, output_stub, top_n = 8L) {
  if (!nrow(result_tbl)) {
    return(invisible(NULL))
  }

  plot_tbl <- result_tbl[result_tbl$database %in% focus_database & !is.na(result_tbl$p.adjust), , drop = FALSE]
  plot_tbl <- plot_tbl[plot_tbl$p.adjust < 0.05, , drop = FALSE]
  if (!nrow(plot_tbl)) {
    return(invisible(NULL))
  }

  split_key <- paste(plot_tbl$input_id, ifelse(is.na(plot_tbl$ontology), plot_tbl$database, paste0(plot_tbl$database, "_", plot_tbl$ontology)), sep = "||")
  split_list <- split(plot_tbl, split_key)
  split_list <- lapply(
    split_list,
    function(df) {
      df <- df[order(df$p.adjust, -df$Count, df$Description), , drop = FALSE]
      utils::head(df, top_n)
    }
  )
  plot_tbl <- do.call(rbind, split_list)
  if (!nrow(plot_tbl)) {
    return(invisible(NULL))
  }

  plot_tbl$facet_label <- ifelse(
    plot_tbl$database == "GO",
    paste(plot_tbl$input_id, plot_tbl$ontology, sep = " | "),
    paste(plot_tbl$input_id, plot_tbl$database, sep = " | ")
  )
  plot_tbl$description_wrapped <- wrap_text(plot_tbl$Description, width = 42)
  plot_tbl$panel_term <- paste0(plot_tbl$description_wrapped, " || ", plot_tbl$facet_label)
  ordered_terms <- unlist(
    lapply(
      split(plot_tbl, plot_tbl$facet_label),
      function(df) {
        df <- df[order(df$p.adjust, -df$Count, df$Description), , drop = FALSE]
        rev(df$panel_term)
      }
    ),
    use.names = FALSE
  )
  plot_tbl$panel_term <- factor(plot_tbl$panel_term, levels = unique(ordered_terms))
  plot_tbl$neg_log10_padj <- -log10(plot_tbl$p.adjust)

  p <- ggplot(plot_tbl, aes(x = gene_ratio_numeric, y = panel_term)) +
    geom_point(aes(size = Count, color = neg_log10_padj), alpha = 0.9) +
    facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
    scale_color_gradient(low = "#3b6c8e", high = "#c94c2c") +
    scale_y_discrete(labels = function(x) sub(" \\|\\| .*", "", x)) +
    labs(
      x = "Gene ratio",
      y = NULL,
      color = "-log10(adj.P)",
      size = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#ebf0f2", color = "#6e7f86"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "#ececec", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )

  pdf(file = paste0(output_stub, ".pdf"), width = 13, height = 10)
  print(p)
  dev.off()

  png(filename = paste0(output_stub, ".png"), width = 2600, height = 2000, res = 220)
  print(p)
  dev.off()
}

dir.create(project_path("res", "tables", "functional"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("res", "qc", "functional"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figure", "export", "functional"), recursive = TRUE, showWarnings = FALSE)

consensus <- read.delim(
  project_path("res", "tables", "bulk", "bulk_strong_glomerular_consensus.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gse30528_genes <- read.delim(
  project_path("data", "processed", "bulk_gene", "GSE30528_gene_expr.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)$gene_symbol

gse96804_genes <- read.delim(
  project_path("data", "processed", "bulk_gene", "GSE96804_gene_expr.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)$gene_symbol

universe_symbols <- sort(intersect(
  unique(gse30528_genes[!is.na(gse30528_genes) & nzchar(gse30528_genes)]),
  unique(gse96804_genes[!is.na(gse96804_genes) & nzchar(gse96804_genes)])
))

universe_map <- map_symbol_to_entrez(universe_symbols)
universe_entrez <- unique(universe_map$ENTREZID)

consensus <- consensus[
  consensus$same_direction_glomerular &
    !is.na(consensus$logFC_GSE30528) &
    !is.na(consensus$logFC_GSE96804),
  ,
  drop = FALSE
]

input_gene_sets <- list(
  strong_glomerular_down = consensus$gene_symbol[consensus$logFC_GSE30528 < 0 & consensus$logFC_GSE96804 < 0],
  strong_glomerular_up = consensus$gene_symbol[consensus$logFC_GSE30528 > 0 & consensus$logFC_GSE96804 > 0]
)

enrichment_runs <- lapply(
  names(input_gene_sets),
  function(input_id) run_one_enrichment(input_gene_sets[[input_id]], universe_entrez, input_id)
)
names(enrichment_runs) <- names(input_gene_sets)

input_summary <- do.call(rbind, lapply(enrichment_runs, `[[`, "summary"))
input_summary$universe_symbols <- length(universe_symbols)
input_summary$universe_entrez <- length(universe_entrez)

result_list <- lapply(enrichment_runs, `[[`, "results")
all_results <- bind_rows_fill(result_list)

if (!nrow(all_results)) {
  stop("No enrichment results were generated.")
}

all_results <- all_results[order(all_results$input_id, all_results$database, all_results$ontology, all_results$p.adjust, decreasing = FALSE), , drop = FALSE]
significant_results <- all_results[!is.na(all_results$p.adjust) & all_results$p.adjust < 0.05, , drop = FALSE]

top_summary <- bind_rows_fill(
  lapply(
    split(significant_results, paste(significant_results$input_id, significant_results$database, ifelse(is.na(significant_results$ontology), "NA", significant_results$ontology), sep = "||")),
    function(df) utils::head(df[order(df$p.adjust, -df$Count, df$Description), c("input_id", "database", "ontology", "ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")], 5)
  )
)

write_tsvgz(
  all_results,
  project_path("res", "tables", "functional", "glomerular_functional_enrichment_full.tsv.gz")
)
write_tsvgz(
  significant_results,
  project_path("res", "tables", "functional", "glomerular_functional_enrichment_significant.tsv.gz")
)
write_tsv(
  input_summary,
  project_path("res", "qc", "functional", "glomerular_functional_input_summary.tsv")
)
write_tsv(
  top_summary,
  project_path("res", "qc", "functional", "glomerular_functional_top_terms_summary.tsv")
)

make_dotplot(
  result_tbl = all_results,
  focus_database = "GO",
  output_stub = project_path("figure", "export", "functional", "glomerular_go_enrichment_dotplot"),
  top_n = 8L
)

make_dotplot(
  result_tbl = all_results,
  focus_database = c("KEGG", "Reactome"),
  output_stub = project_path("figure", "export", "functional", "glomerular_pathway_enrichment_dotplot"),
  top_n = 8L
)
