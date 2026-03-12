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
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

read_geo_platform_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin_idx <- which(lines == "!platform_table_begin")
  end_idx <- which(lines == "!platform_table_end")

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop("Could not find a valid platform table block in: ", path)
  }

  block <- lines[(begin_idx + 1):(end_idx - 1)]
  handle <- textConnection(block)
  on.exit(close(handle), add = TRUE)

  read.delim(
    handle,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
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

first_non_empty <- function(values) {
  values <- trimws(values)
  values <- values[!(is.na(values) | values == "" | values == "---" | values == "NULL")]
  if (length(values) == 0) {
    return(NA_character_)
  }
  values[[1]]
}

collapse_unique_values <- function(values) {
  values <- trimws(values)
  values <- unique(values[!(is.na(values) | values == "" | values == "---" | values == "NULL")])
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(values, collapse = " /// ")
}

parse_assignment_field <- function(field_value) {
  if (is.na(field_value) || field_value == "" || field_value == "---") {
    return(c(
      gene_symbol = NA_character_,
      gene_symbol_all = NA_character_,
      entrez_id = NA_character_,
      entrez_id_all = NA_character_,
      gene_title = NA_character_,
      gene_title_all = NA_character_
    ))
  }

  entries <- strsplit(field_value, " /// ", fixed = TRUE)[[1]]
  parts <- lapply(entries, function(entry) {
    trimws(strsplit(entry, " // ", fixed = TRUE)[[1]])
  })

  symbols <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
  titles <- vapply(parts, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1))
  entrez <- vapply(parts, function(x) if (length(x) >= 5) x[[5]] else NA_character_, character(1))

  c(
    gene_symbol = first_non_empty(symbols),
    gene_symbol_all = collapse_unique_values(symbols),
    entrez_id = first_non_empty(entrez),
    entrez_id_all = collapse_unique_values(entrez),
    gene_title = first_non_empty(titles),
    gene_title_all = collapse_unique_values(titles)
  )
}

pick_first_token <- function(x) {
  ifelse(
    is.na(x) | x == "",
    NA_character_,
    trimws(vapply(strsplit(x, " /// ", fixed = TRUE), `[`, character(1), 1))
  )
}

load_platform <- function(gpl_id) {
  path <- project_path("data", "raw", "platform_tables", paste0(gpl_id, "_platform_table.txt"))
  read_geo_platform_table(path)
}

safe_load_platform <- function(gpl_id, required = TRUE) {
  tryCatch(
    load_platform(gpl_id),
    error = function(e) {
      if (required) {
        stop(e)
      }

      message("Skipping optional platform ", gpl_id, ": ", conditionMessage(e))
      NULL
    }
  )
}

gpl571 <- load_platform("GPL571")
gpl571_map <- data.frame(
  platform_id = "GPL571",
  feature_id = gpl571$ID,
  gene_symbol = pick_first_token(gpl571$`Gene Symbol`),
  gene_symbol_all = gpl571$`Gene Symbol`,
  entrez_id = pick_first_token(gpl571$ENTREZ_GENE_ID),
  entrez_id_all = gpl571$ENTREZ_GENE_ID,
  gene_title = pick_first_token(gpl571$`Gene Title`),
  gene_title_all = gpl571$`Gene Title`,
  raw_annotation = "Gene Symbol / ENTREZ_GENE_ID / Gene Title",
  stringsAsFactors = FALSE
)

gpl17586 <- load_platform("GPL17586")
gpl17586_parsed <- t(vapply(gpl17586$gene_assignment, parse_assignment_field, character(6)))
gpl17586_map <- data.frame(
  platform_id = "GPL17586",
  feature_id = gpl17586$ID,
  gene_symbol = gpl17586_parsed[, "gene_symbol"],
  gene_symbol_all = gpl17586_parsed[, "gene_symbol_all"],
  entrez_id = gpl17586_parsed[, "entrez_id"],
  entrez_id_all = gpl17586_parsed[, "entrez_id_all"],
  gene_title = gpl17586_parsed[, "gene_title"],
  gene_title_all = gpl17586_parsed[, "gene_title_all"],
  raw_annotation = gpl17586$gene_assignment,
  stringsAsFactors = FALSE
)

gpl24120 <- load_platform("GPL24120")
gpl24120_entrez <- as.character(gpl24120$ENTREZ_GENE_ID)
gpl24120_symbol <- mapIds(
  org.Hs.eg.db,
  keys = gpl24120_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
gpl24120_title <- mapIds(
  org.Hs.eg.db,
  keys = gpl24120_entrez,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)
gpl24120_map <- data.frame(
  platform_id = "GPL24120",
  feature_id = gpl24120$ID,
  gene_symbol = as.character(gpl24120_symbol),
  gene_symbol_all = as.character(gpl24120_symbol),
  entrez_id = gpl24120_entrez,
  entrez_id_all = gpl24120_entrez,
  gene_title = ifelse(is.na(gpl24120_title), gpl24120$Description, as.character(gpl24120_title)),
  gene_title_all = ifelse(is.na(gpl24120_title), gpl24120$Description, as.character(gpl24120_title)),
  raw_annotation = gpl24120$Description,
  stringsAsFactors = FALSE
)

gpl22945 <- load_platform("GPL22945")
gpl22945_map <- data.frame(
  platform_id = "GPL22945",
  feature_id = gpl22945$ID,
  gene_symbol = gpl22945$Symbol,
  gene_symbol_all = gpl22945$Symbol,
  entrez_id = as.character(gpl22945$ENTREZ_GENE_ID),
  entrez_id_all = as.character(gpl22945$ENTREZ_GENE_ID),
  gene_title = gpl22945$SPOT_ID,
  gene_title_all = gpl22945$SPOT_ID,
  raw_annotation = paste(gpl22945$Symbol, gpl22945$SPOT_ID, sep = " | "),
  stringsAsFactors = FALSE
)

gpl23126 <- safe_load_platform("GPL23126", required = FALSE)
gpl23126_map <- NULL

if (!is.null(gpl23126)) {
  gpl23126_parsed <- t(vapply(gpl23126$gene_assignment, parse_assignment_field, character(6)))
  gpl23126_map <- data.frame(
    platform_id = "GPL23126",
    feature_id = gpl23126$ID,
    gene_symbol = gpl23126_parsed[, "gene_symbol"],
    gene_symbol_all = gpl23126_parsed[, "gene_symbol_all"],
    entrez_id = gpl23126_parsed[, "entrez_id"],
    entrez_id_all = gpl23126_parsed[, "entrez_id_all"],
    gene_title = gpl23126_parsed[, "gene_title"],
    gene_title_all = gpl23126_parsed[, "gene_title_all"],
    raw_annotation = gpl23126$gene_assignment,
    stringsAsFactors = FALSE
  )
}

gse104954_combined <- merge(
  gpl24120_map[, c("feature_id", "gene_symbol", "entrez_id", "gene_title")],
  gpl22945_map[, c("feature_id", "gene_symbol", "entrez_id", "gene_title")],
  by = "feature_id",
  all = TRUE,
  suffixes = c("_gpl24120", "_gpl22945")
)

gse104954_combined$gene_symbol <- ifelse(
  is.na(gse104954_combined$gene_symbol_gpl22945) | gse104954_combined$gene_symbol_gpl22945 == "",
  gse104954_combined$gene_symbol_gpl24120,
  gse104954_combined$gene_symbol_gpl22945
)
gse104954_combined$entrez_id <- ifelse(
  is.na(gse104954_combined$entrez_id_gpl22945) | gse104954_combined$entrez_id_gpl22945 == "",
  gse104954_combined$entrez_id_gpl24120,
  gse104954_combined$entrez_id_gpl22945
)
gse104954_combined$gene_title <- ifelse(
  is.na(gse104954_combined$gene_title_gpl22945) | gse104954_combined$gene_title_gpl22945 == "",
  gse104954_combined$gene_title_gpl24120,
  gse104954_combined$gene_title_gpl22945
)
gse104954_combined$platform_id <- "GSE104954_merged"

mapping_list <- list(
  GPL571 = gpl571_map,
  GPL17586 = gpl17586_map,
  GPL24120 = gpl24120_map,
  GPL22945 = gpl22945_map,
  GSE104954_merged = gse104954_combined[, c("platform_id", "feature_id", "gene_symbol", "entrez_id", "gene_title")]
)

if (!is.null(gpl23126_map)) {
  mapping_list$GPL23126 <- gpl23126_map
}

for (mapping_name in names(mapping_list)) {
  mapping <- mapping_list[[mapping_name]]
  write_tsv_gz(
    mapping,
    project_path("data", "processed", "annotation", paste0(mapping_name, "_mapping.tsv.gz"))
  )
}

summary_table <- do.call(
  rbind,
  lapply(names(mapping_list), function(mapping_name) {
    mapping <- mapping_list[[mapping_name]]
    data.frame(
      mapping_id = mapping_name,
      n_rows = nrow(mapping),
      n_nonempty_gene_symbol = sum(!(is.na(mapping$gene_symbol) | mapping$gene_symbol == "")),
      n_nonempty_entrez_id = sum(!(is.na(mapping$entrez_id) | mapping$entrez_id == "")),
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  summary_table,
  file = project_path("res", "qc", "bulk", "platform_annotation_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Processed platform mappings written to data/processed/annotation\n")
cat("Annotation summary written to res/qc/bulk/platform_annotation_summary.tsv\n")
print(summary_table)
