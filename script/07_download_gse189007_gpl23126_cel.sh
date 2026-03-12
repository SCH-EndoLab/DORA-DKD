#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
SAMPLE_SHEET="${ROOT_DIR}/data/metadata/soft_samples/GSE189007_samples.tsv"
OUT_DIR="${ROOT_DIR}/data/raw/cel/GSE189007_GPL23126"
LOG_FILE="${ROOT_DIR}/Log/07_download_gse189007_gpl23126_cel.log"

mkdir -p "${OUT_DIR}"

Rscript - <<'RS' "${SAMPLE_SHEET}" "${OUT_DIR}" "${LOG_FILE}"
args <- commandArgs(trailingOnly = TRUE)
sample_sheet <- args[1]
out_dir <- args[2]
log_file <- args[3]

df <- read.delim(sample_sheet, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
sub <- subset(df, platform_id == "GPL23126")

extract_cel_url <- function(x) {
  parts <- strsplit(x, " | ", fixed = TRUE)[[1]]
  hits <- parts[grepl("\\.CEL\\.gz$", parts, ignore.case = TRUE)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[1]
}

sub$cel_url <- vapply(sub$supplementary_file, extract_cel_url, character(1))
sub$out_file <- file.path(out_dir, basename(sub$cel_url))

con <- file(log_file, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("[gse189007_gpl23126_cel_download] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), con)

for (i in seq_len(nrow(sub))) {
  url <- sub$cel_url[i]
  out_file <- sub$out_file[i]
  sample_id <- sub$sample_id[i]

  if (is.na(url) || url == "") {
    writeLines(sprintf("SKIP\t%s\tmissing_cel_url", sample_id), con)
    next
  }

  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    writeLines(sprintf("EXISTS\t%s\t%s", sample_id, out_file), con)
    next
  }

  writeLines(sprintf("DOWNLOAD\t%s\t%s", sample_id, url), con)
  status <- system2(
    "curl",
    c("-L", "--fail", "--silent", "--show-error", url, "-o", out_file),
    stdout = FALSE,
    stderr = FALSE
  )

  if (!identical(status, 0L) || !file.exists(out_file) || file.info(out_file)$size <= 0) {
    if (file.exists(out_file)) {
      unlink(out_file)
    }
    writeLines(sprintf("ERROR\t%s\t%s", sample_id, url), con)
    stop("Failed downloading ", sample_id, " from ", url)
  }

  writeLines(sprintf("DONE\t%s\t%s", sample_id, out_file), con)
}
RS

