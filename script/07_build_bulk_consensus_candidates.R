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

gse30528 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE30528_limma_dkd_vs_control.tsv.gz"))
gse96804 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE96804_limma_dkd_vs_control.tsv.gz"))
gse104954 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE104954_primary_dkd_vs_tumor_limma_dkd_vs_control.tsv.gz"))

g1 <- gse30528[, c("gene_symbol", "feature_id", "logFC", "P.Value", "adj.P.Val")]
names(g1)[-1] <- paste0(names(g1)[-1], "_GSE30528")

g2 <- gse96804[, c("gene_symbol", "feature_id", "logFC", "P.Value", "adj.P.Val")]
names(g2)[-1] <- paste0(names(g2)[-1], "_GSE96804")

g3 <- gse104954[, c("gene_symbol", "feature_id", "logFC", "P.Value", "adj.P.Val")]
names(g3)[-1] <- paste0(names(g3)[-1], "_GSE104954")

glomerular_consensus <- merge(g1, g2, by = "gene_symbol", all = FALSE)
glomerular_consensus <- glomerular_consensus[
  !(is.na(glomerular_consensus$P.Value_GSE30528) | is.na(glomerular_consensus$P.Value_GSE96804)),
  ,
  drop = FALSE
]

glomerular_consensus$same_direction_glomerular <- sign(glomerular_consensus$logFC_GSE30528) == sign(glomerular_consensus$logFC_GSE96804)
glomerular_consensus$fisher_stat <- -2 * (
  log(glomerular_consensus$P.Value_GSE30528) +
    log(glomerular_consensus$P.Value_GSE96804)
)
glomerular_consensus$fisher_p <- pchisq(glomerular_consensus$fisher_stat, df = 4, lower.tail = FALSE)
glomerular_consensus$fisher_fdr <- p.adjust(glomerular_consensus$fisher_p, method = "BH")

consensus_with_tubule <- merge(glomerular_consensus, g3, by = "gene_symbol", all.x = TRUE)
consensus_with_tubule$tubule_same_direction <- sign(consensus_with_tubule$logFC_GSE30528) == sign(consensus_with_tubule$logFC_GSE104954)
consensus_with_tubule$tubule_support <- ifelse(
  is.na(consensus_with_tubule$logFC_GSE104954),
  "missing_in_tubule",
  ifelse(
    consensus_with_tubule$tubule_same_direction & consensus_with_tubule$adj.P.Val_GSE104954 < 0.05,
    "same_direction_fdr",
    ifelse(
      consensus_with_tubule$tubule_same_direction & consensus_with_tubule$P.Value_GSE104954 < 0.05,
      "same_direction_nominal",
      ifelse(
        consensus_with_tubule$tubule_same_direction,
        "same_direction_only",
        "opposite_direction"
      )
    )
  )
)

consensus_with_tubule <- consensus_with_tubule[order(consensus_with_tubule$fisher_fdr, -abs(consensus_with_tubule$logFC_GSE96804)), ]

strong_glomerular <- consensus_with_tubule[
  consensus_with_tubule$same_direction_glomerular & consensus_with_tubule$fisher_fdr < 0.05,
  ,
  drop = FALSE
]

write_tsv_gz(
  consensus_with_tubule,
  project_path("res", "tables", "bulk", "bulk_glomerular_consensus_with_tubule_support.tsv.gz")
)

write_tsv_gz(
  strong_glomerular,
  project_path("res", "tables", "bulk", "bulk_strong_glomerular_consensus.tsv.gz")
)

summary_table <- data.frame(
  metric = c(
    "all_glomerular_overlap_genes",
    "same_direction_glomerular",
    "same_direction_glomerular_fdr_lt_0.05",
    "same_direction_glomerular_fdr_lt_0.05_and_tubule_fdr",
    "same_direction_glomerular_fdr_lt_0.05_and_tubule_nominal",
    "same_direction_glomerular_fdr_lt_0.05_and_tubule_same_direction_only"
  ),
  value = c(
    nrow(consensus_with_tubule),
    sum(consensus_with_tubule$same_direction_glomerular, na.rm = TRUE),
    nrow(strong_glomerular),
    sum(strong_glomerular$tubule_support == "same_direction_fdr", na.rm = TRUE),
    sum(strong_glomerular$tubule_support == "same_direction_nominal", na.rm = TRUE),
    sum(strong_glomerular$tubule_support == "same_direction_only", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = project_path("res", "qc", "bulk", "bulk_consensus_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Consensus tables written to res/tables/bulk\n")
cat("Consensus summary written to res/qc/bulk/bulk_consensus_summary.tsv\n")
print(summary_table)
