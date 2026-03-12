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

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

normalize_text <- function(x) {
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  x <- trimws(x)
  x
}

extract_field_value <- function(line) {
  parts <- strsplit(line, " = ", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return("")
  }

  normalize_text(paste(parts[-1], collapse = " = "))
}

finalize_record <- function(rec, out_list) {
  if (is.null(rec$sample_id) || identical(rec$sample_id, "")) {
    return(out_list)
  }

  supp_files <- rec$supplementary_files
  if (length(supp_files) == 0) {
    supp_files <- NA_character_
  }

  relations <- rec$relations
  if (length(relations) == 0) {
    relations <- NA_character_
  }

  characteristics <- rec$characteristics
  if (length(characteristics) == 0) {
    characteristics <- NA_character_
  }

  data.frame(
    dataset_id = rec$dataset_id,
    sample_id = rec$sample_id,
    sample_title = rec$sample_title %||% NA_character_,
    sample_type = rec$sample_type %||% NA_character_,
    source_name = rec$source_name %||% NA_character_,
    platform_id = rec$platform_id %||% NA_character_,
    characteristic_summary = paste(characteristics, collapse = " | "),
    relation_summary = paste(relations, collapse = " | "),
    supplementary_file = supp_files,
    file_name = basename(supp_files),
    stringsAsFactors = FALSE
  ) |>
    split(seq_len(length(supp_files))) |>
    c(out_list)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) {
    return(y)
  }
  x
}

parse_soft_manifest <- function(dataset_id) {
  soft_path <- project_path("data", "raw", "geo_soft", paste0(dataset_id, "_family.soft.gz"))
  lines <- readLines(gzfile(soft_path), warn = FALSE, encoding = "UTF-8")

  out <- list()
  rec <- list(
    dataset_id = dataset_id,
    supplementary_files = character(),
    relations = character(),
    characteristics = character()
  )

  for (line in lines) {
    if (startsWith(line, "^SAMPLE = ")) {
      out <- finalize_record(rec, out)
      rec <- list(
        dataset_id = dataset_id,
        sample_id = extract_field_value(line),
        supplementary_files = character(),
        relations = character(),
        characteristics = character()
      )
      next
    }

    if (startsWith(line, "!Sample_title = ")) {
      rec$sample_title <- extract_field_value(line)
      next
    }

    if (startsWith(line, "!Sample_type = ")) {
      rec$sample_type <- extract_field_value(line)
      next
    }

    if (startsWith(line, "!Sample_source_name_ch1 = ")) {
      rec$source_name <- extract_field_value(line)
      next
    }

    if (startsWith(line, "!Sample_platform_id = ")) {
      rec$platform_id <- extract_field_value(line)
      next
    }

    if (startsWith(line, "!Sample_characteristics_ch1 = ")) {
      rec$characteristics <- c(rec$characteristics, extract_field_value(line))
      next
    }

    if (startsWith(line, "!Sample_relation = ")) {
      rec$relations <- c(rec$relations, extract_field_value(line))
      next
    }

    if (grepl("^!Sample_supplementary_file_[0-9]+ = ", line)) {
      rec$supplementary_files <- c(rec$supplementary_files, extract_field_value(line))
      next
    }
  }

  out <- finalize_record(rec, out)
  if (length(out) == 0) {
    return(data.frame())
  }

  manifest <- do.call(rbind, out)
  manifest$file_type <- ifelse(
    grepl("processed", manifest$file_name, ignore.case = TRUE),
    "processed",
    ifelse(grepl("matrix|barcodes|features|tar\\.gz", manifest$file_name, ignore.case = TRUE), "matrix_or_archive", "other")
  )
  manifest
}

datasets <- c("GSE209781", "GSE279086", "GSE266146")
manifest_list <- lapply(datasets, parse_soft_manifest)
manifest <- do.call(rbind, manifest_list)
manifest <- manifest[order(manifest$dataset_id, manifest$sample_id, manifest$file_name), ]
row.names(manifest) <- NULL

ensure_dir(project_path("data", "metadata"))
ensure_dir(project_path("res", "qc", "single_cell"))

write.table(
  manifest,
  file = project_path("data", "metadata", "single_cell_supplementary_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_df <- aggregate(
  list(file_count = manifest$file_name),
  by = list(dataset_id = manifest$dataset_id, sample_id = manifest$sample_id, file_type = manifest$file_type),
  FUN = length
)
names(summary_df)[names(summary_df) == "file_count"] <- "n_files"

dataset_summary <- aggregate(
  list(
    n_samples = manifest$sample_id,
    n_files = manifest$file_name
  ),
  by = list(dataset_id = manifest$dataset_id),
  FUN = function(x) length(unique(x))
)
names(dataset_summary)[2:3] <- c("n_samples", "n_files")

write.table(
  summary_df,
  file = project_path("res", "qc", "single_cell", "single_cell_supplementary_file_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  dataset_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_supplementary_dataset_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell supplementary manifest written to data/metadata/single_cell_supplementary_manifest.tsv\n")
cat("Single-cell supplementary summaries written to res/qc/single_cell/\n")
