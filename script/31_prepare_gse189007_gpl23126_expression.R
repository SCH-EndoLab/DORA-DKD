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
  library(oligo)
  library(pd.clariom.d.human)
  library(clariomdhumantranscriptcluster.db)
  library(AnnotationDbi)
})

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

derive_diabetes_group <- function(title_text) {
  text <- tolower(title_text)

  if (grepl("^control", text)) {
    return("Control")
  }
  if (grepl("without complications within 5 years", text)) {
    return("T2D_NoComp_LT5Y")
  }
  if (grepl("without complications above 15 years", text)) {
    return("T2D_NoComp_GT15Y")
  }
  if (grepl("nephropathy within 5 years", text)) {
    return("T2DN_LT5Y")
  }
  if (grepl("nephropathy above 15 years", text)) {
    return("T2DN_GT15Y")
  }
  if (grepl("retinopathy within 5 years", text)) {
    return("T2DR_LT5Y")
  }
  if (grepl("retinopathy above 15 years", text)) {
    return("T2DR_GT15Y")
  }
  NA_character_
}

derive_disease_axis <- function(group_label) {
  if (is.na(group_label)) {
    return(NA_character_)
  }
  if (group_label == "Control") {
    return("Control")
  }
  if (grepl("^T2D_NoComp", group_label)) {
    return("T2D_NoComp")
  }
  if (grepl("^T2DN_", group_label)) {
    return("T2DN")
  }
  if (grepl("^T2DR_", group_label)) {
    return("T2DR")
  }
  NA_character_
}

collapse_to_gene_level <- function(expr_df, mapping_df) {
  sample_ids <- setdiff(colnames(expr_df), c("feature_id", "gene_symbol"))
  merged <- merge(expr_df, mapping_df, by = "feature_id", all.x = TRUE, sort = FALSE)
  merged <- merged[!(is.na(merged$gene_symbol) | merged$gene_symbol == ""), , drop = FALSE]

  if (nrow(merged) == 0) {
    stop("No mapped transcript clusters retained for GSE189007 GPL23126.")
  }

  row_variance <- apply(merged[, sample_ids, drop = FALSE], 1, var, na.rm = TRUE)
  merged$row_variance <- row_variance
  split_index <- split(seq_len(nrow(merged)), merged$gene_symbol)
  keep_index <- vapply(split_index, function(idx) idx[which.max(merged$row_variance[idx])], integer(1))

  collapsed <- merged[keep_index, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]
  collapsed <- collapsed[order(collapsed$gene_symbol), , drop = FALSE]
  rownames(collapsed) <- NULL
  collapsed
}

sample_sheet <- read.delim(
  project_path("data", "metadata", "soft_samples", "GSE189007_samples.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pheno <- subset(sample_sheet, platform_id == "GPL23126")
pheno$group_std <- vapply(pheno$title, derive_diabetes_group, character(1))
pheno$disease_axis_group <- vapply(pheno$group_std, derive_disease_axis, character(1))
pheno$keep_diabetes_axis <- ifelse(pheno$disease_axis_group %in% c("Control", "T2D_NoComp", "T2DN", "T2DR"), "yes", "no")

cel_dir <- project_path("data", "raw", "cel", "GSE189007_GPL23126")
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)

if (length(cel_files) == 0) {
  stop("No GPL23126 CEL.gz files found under ", cel_dir)
}

sample_ids_from_file <- sub("_.*$", "", basename(cel_files))
cel_manifest <- data.frame(
  sample_id = sample_ids_from_file,
  cel_file = cel_files,
  stringsAsFactors = FALSE
)

pheno <- merge(pheno, cel_manifest, by = "sample_id", all.x = FALSE, sort = FALSE)
pheno <- pheno[match(sample_ids_from_file, pheno$sample_id), , drop = FALSE]

raw_data <- read.celfiles(pheno$cel_file, pkgname = "pd.clariom.d.human")
norm_eset <- rma(raw_data, target = "core")
expr_mat <- exprs(norm_eset)

expr_df <- data.frame(
  feature_id = rownames(expr_mat),
  expr_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(expr_df)[-1] <- sampleNames(norm_eset)

sample_lookup <- data.frame(
  sample_name = sampleNames(norm_eset),
  cel_base = basename(pheno$cel_file),
  sample_id = pheno$sample_id,
  stringsAsFactors = FALSE
)

col_order_ids <- sample_lookup$sample_id[match(colnames(expr_df)[-1], sample_lookup$sample_name)]
colnames(expr_df) <- c("feature_id", col_order_ids)

mapping <- AnnotationDbi::select(
  clariomdhumantranscriptcluster.db,
  keys = rownames(expr_mat),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)
mapping <- mapping[!duplicated(mapping$PROBEID), , drop = FALSE]
colnames(mapping) <- c("feature_id", "gene_symbol")

gene_expr <- collapse_to_gene_level(expr_df, mapping)

out_bulk_dir <- project_path("data", "processed", "bulk")
out_gene_dir <- project_path("data", "processed", "bulk_gene")

write_tsv_gz(expr_df, file.path(out_bulk_dir, "GSE189007_GPL23126_expr.tsv.gz"))
write_tsv(pheno, file.path(out_bulk_dir, "GSE189007_GPL23126_pheno.tsv"))
write_tsv_gz(gene_expr, file.path(out_gene_dir, "GSE189007_GPL23126_gene_expr.tsv.gz"))

summary_df <- data.frame(
  dataset_id = "GSE189007_GPL23126",
  n_samples = ncol(expr_df) - 1,
  n_core_transcript_clusters = nrow(expr_df),
  n_gene_level_rows = nrow(gene_expr),
  groups = paste(sort(unique(pheno$group_std)), collapse = ";"),
  disease_axis_groups = paste(sort(unique(pheno$disease_axis_group)), collapse = ";"),
  stringsAsFactors = FALSE
)

write_tsv(summary_df, project_path("res", "qc", "bulk", "GSE189007_GPL23126_prep_summary.tsv"))

cat("GSE189007 GPL23126 probe-level matrix written to data/processed/bulk/GSE189007_GPL23126_expr.tsv.gz\n")
cat("GSE189007 GPL23126 phenotype table written to data/processed/bulk/GSE189007_GPL23126_pheno.tsv\n")
cat("GSE189007 GPL23126 gene-level matrix written to data/processed/bulk_gene/GSE189007_GPL23126_gene_expr.tsv.gz\n")
cat("GSE189007 GPL23126 prep summary written to res/qc/bulk/GSE189007_GPL23126_prep_summary.tsv\n")
