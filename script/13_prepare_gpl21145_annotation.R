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

suppressPackageStartupMessages(library(data.table))

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

annotation_path <- project_path(
  "data", "raw", "platform_annotation", "GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
)

ann <- fread(annotation_path, skip = 7, data.table = FALSE, fill = TRUE)
colnames(ann)[1] <- sub("^\\ufeff", "", colnames(ann)[1])

if ("V47" %in% colnames(ann)) {
  ann$V47 <- NULL
}

keep_cols <- intersect(
  c(
    "IlmnID", "Name", "CHR", "MAPINFO", "UCSC_RefGene_Name",
    "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island",
    "UCSC_CpG_Islands_Name", "Phantom4_Enhancers", "Phantom5_Enhancers",
    "Regulatory_Feature_Name", "Regulatory_Feature_Group"
  ),
  colnames(ann)
)
ann <- ann[, keep_cols, drop = FALSE]
ann$cpg_id <- ann$Name
ann$gene_symbol_primary <- vapply(
  strsplit(ifelse(is.na(ann$UCSC_RefGene_Name), "", ann$UCSC_RefGene_Name), ";", fixed = TRUE),
  function(x) {
    x <- trimws(x)
    x <- x[x != ""]
    if (length(x) == 0) NA_character_ else x[[1]]
  },
  character(1)
)

ann_out <- ann[, c(
  "cpg_id", "IlmnID", "CHR", "MAPINFO", "gene_symbol_primary",
  "UCSC_RefGene_Name", "UCSC_RefGene_Group",
  "Relation_to_UCSC_CpG_Island", "UCSC_CpG_Islands_Name",
  "Phantom4_Enhancers", "Phantom5_Enhancers",
  "Regulatory_Feature_Name", "Regulatory_Feature_Group"
)]

write_tsv_gz(
  ann_out,
  project_path("data", "processed", "GPL21145_annotation.tsv.gz")
)

summary_table <- data.frame(
  n_rows = nrow(ann_out),
  n_nonempty_gene_symbol = sum(!(is.na(ann_out$gene_symbol_primary) | ann_out$gene_symbol_primary == "")),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = project_path("res", "qc", "mechanism", "GPL21145_annotation_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("GPL21145 annotation written to data/processed/GPL21145_annotation.tsv.gz\n")
print(summary_table)
