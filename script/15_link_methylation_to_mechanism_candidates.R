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

meth <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE233758_methylation_limma.tsv.gz"))
gene_sets <- read_tsv_auto(project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"))
candidate_rank <- read_tsv_auto(project_path("res", "tables", "mechanism", "mechanism_linked_candidate_ranking.tsv.gz"))

mechanism_genes <- unique(gene_sets$gene_symbol)
top_candidate_genes <- unique(candidate_rank$gene_symbol[1:min(100, nrow(candidate_rank))])

split_gene_symbols <- function(x) {
  x <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE)))
  unique(x[x != ""])
}

gene_symbol_list <- lapply(meth$UCSC_RefGene_Name, split_gene_symbols)
nonempty_idx <- lengths(gene_symbol_list) > 0
expanded_map <- data.frame(
  cpg_id = rep(meth$cpg_id[nonempty_idx], lengths(gene_symbol_list[nonempty_idx])),
  gene_symbol = unlist(gene_symbol_list[nonempty_idx], use.names = FALSE),
  stringsAsFactors = FALSE
)

mechanism_map <- expanded_map[expanded_map$gene_symbol %in% mechanism_genes, , drop = FALSE]
candidate_map <- expanded_map[expanded_map$gene_symbol %in% top_candidate_genes, , drop = FALSE]

mechanism_cpg <- merge(meth, unique(mechanism_map[, "cpg_id", drop = FALSE]), by = "cpg_id", all = FALSE)
candidate_cpg <- merge(meth, unique(candidate_map[, "cpg_id", drop = FALSE]), by = "cpg_id", all = FALSE)

gene_level_summary <- function(df, gene_map, gene_pool, label) {
  rows <- list()
  for (gene in gene_pool) {
    cpg_ids <- unique(gene_map$cpg_id[gene_map$gene_symbol == gene])
    sub <- df[df$cpg_id %in% cpg_ids, , drop = FALSE]
    if (nrow(sub) == 0) {
      next
    }
    best <- sub[order(sub$adj.P.Val, -abs(sub$delta_beta_dkd_minus_control)), ][1, , drop = FALSE]
    rows[[length(rows) + 1]] <- data.frame(
      table_type = label,
      gene_symbol = gene,
      n_cpg = nrow(sub),
      best_cpg = best$cpg_id,
      best_adj_p = best$adj.P.Val,
      best_delta_beta = best$delta_beta_dkd_minus_control,
      best_dm_direction = best$dm_direction,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

mechanism_gene_summary <- gene_level_summary(meth, mechanism_map, sort(unique(mechanism_genes)), "mechanism_gene")
candidate_gene_summary <- gene_level_summary(meth, candidate_map, sort(unique(top_candidate_genes)), "top_candidate_gene")

write_tsv_gz(
  mechanism_cpg,
  project_path("res", "tables", "mechanism", "GSE233758_mechanism_gene_cpgs.tsv.gz")
)
write_tsv_gz(
  candidate_cpg,
  project_path("res", "tables", "mechanism", "GSE233758_top_candidate_gene_cpgs.tsv.gz")
)
write_tsv(
  mechanism_gene_summary,
  project_path("res", "qc", "mechanism", "GSE233758_mechanism_gene_methylation_summary.tsv")
)
write_tsv(
  candidate_gene_summary,
  project_path("res", "qc", "mechanism", "GSE233758_top_candidate_gene_methylation_summary.tsv")
)

summary_table <- data.frame(
  metric = c(
    "mechanism_related_cpg_rows",
    "top_candidate_related_cpg_rows",
    "mechanism_genes_with_at_least_one_cpg",
    "top_candidate_genes_with_at_least_one_cpg"
  ),
  value = c(
    nrow(mechanism_cpg),
    nrow(candidate_cpg),
    nrow(mechanism_gene_summary),
    nrow(candidate_gene_summary)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  project_path("res", "qc", "mechanism", "GSE233758_methylation_link_summary.tsv")
)

cat("Mechanism-related methylation tables written to res/tables/mechanism\n")
print(summary_table)
