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

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

normalize_text <- function(x) {
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  trimws(x)
}

derive_group_from_title <- function(x) {
  x <- toupper(normalize_text(x))
  ifelse(
    grepl("^DKD", x),
    "DKD",
    ifelse(grepl("^NM", x), "Control", NA_character_)
  )
}

derive_group_from_disease_state <- function(x) {
  x <- tolower(normalize_text(x))
  ifelse(
    grepl("diabetic kidney disease", x),
    "DKD",
    ifelse(grepl("non-diabetic", x) | grepl("normal", x) | grepl("control", x), "Control", NA_character_)
  )
}

gse209781 <- read_tsv(project_path("data", "metadata", "soft_samples", "GSE209781_samples.tsv"))
gse279086 <- read_tsv(project_path("data", "metadata", "soft_samples", "GSE279086_samples.tsv"))
gse266146 <- read_tsv(project_path("data", "metadata", "soft_samples", "GSE266146_samples.tsv"))

gse209781_manifest <- data.frame(
  dataset_id = "GSE209781",
  sample_id = gse209781$sample_id,
  sample_title = gse209781$title,
  tissue = normalize_text(gse209781$char_tissue),
  disease_label_raw = normalize_text(gse209781$char_disease_state),
  group_from_title = derive_group_from_title(gse209781$title),
  group_from_metadata = derive_group_from_disease_state(gse209781$char_disease_state),
  stringsAsFactors = FALSE
)
gse209781_manifest$group_std <- ifelse(
  !is.na(gse209781_manifest$group_from_title),
  gse209781_manifest$group_from_title,
  gse209781_manifest$group_from_metadata
)
gse209781_manifest$diabetes_context <- "DKD_kidney_cortex"
gse209781_manifest$sample_origin <- "Kidney cortex"
gse209781_manifest$analysis_role <- "primary_single_cell_kidney"
gse209781_manifest$title_metadata_conflict <- with(
  gse209781_manifest,
  !is.na(group_from_title) &
    !is.na(group_from_metadata) &
    group_from_title != group_from_metadata
)
gse209781_manifest$notes <- ifelse(
  gse209781_manifest$title_metadata_conflict,
  "Title prefix conflicts with GEO disease-state field; using title-based grouping pending expression-level verification.",
  "Title and metadata are concordant."
)

gse279086_manifest <- data.frame(
  dataset_id = "GSE279086",
  sample_id = gse279086$sample_id,
  sample_title = gse279086$title,
  tissue = normalize_text(gse279086$char_tissue),
  disease_label_raw = normalize_text(gse279086$char_disease),
  group_from_title = NA_character_,
  group_from_metadata = ifelse(
    normalize_text(gse279086$char_disease) == "Lean Control",
    "Control",
    ifelse(normalize_text(gse279086$char_disease) == "Type 1 Diabetes", "Diabetes", NA_character_)
  ),
  stringsAsFactors = FALSE
)
gse279086_manifest$group_std <- gse279086_manifest$group_from_metadata
gse279086_manifest$diabetes_context <- ifelse(
  gse279086_manifest$group_std == "Diabetes",
  "T1D_kidney_cortex",
  "non_diabetic_control_kidney_cortex"
)
gse279086_manifest$sample_origin <- "Kidney cortex"
gse279086_manifest$analysis_role <- "external_diabetes_kidney_validation"
gse279086_manifest$title_metadata_conflict <- FALSE
gse279086_manifest$notes <- "Cortex scRNA-seq from young adults with type 1 diabetes or lean controls."

gse266146_manifest <- data.frame(
  dataset_id = "GSE266146",
  sample_id = gse266146$sample_id,
  sample_title = gse266146$title,
  tissue = normalize_text(gse266146$char_tissue),
  disease_label_raw = "Early and advanced DKD urine sediment pooled by batch",
  group_from_title = NA_character_,
  group_from_metadata = "DKD",
  stringsAsFactors = FALSE
)
gse266146_manifest$group_std <- "DKD_urine_pool"
gse266146_manifest$diabetes_context <- "DKD_urine_early_advanced_mixed"
gse266146_manifest$sample_origin <- "Urine sediment"
gse266146_manifest$analysis_role <- "external_urine_validation"
gse266146_manifest$title_metadata_conflict <- FALSE
gse266146_manifest$notes <- "Sample-level labels are batch identifiers only; series metadata indicates urine sediment from early and advanced DKD."

single_cell_manifest <- rbind(
  gse209781_manifest,
  gse279086_manifest,
  gse266146_manifest
)

single_cell_manifest$priority_rank <- ifelse(
  single_cell_manifest$dataset_id == "GSE209781",
  1,
  ifelse(single_cell_manifest$dataset_id == "GSE266146", 2, 3)
)
single_cell_manifest <- single_cell_manifest[
  order(single_cell_manifest$priority_rank, single_cell_manifest$dataset_id, single_cell_manifest$sample_id),
]
row.names(single_cell_manifest) <- NULL

ensure_dir(project_path("data", "metadata"))
ensure_dir(project_path("res", "qc", "single_cell"))

write.table(
  single_cell_manifest,
  file = project_path("data", "metadata", "single_cell_study_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

group_summary <- aggregate(
  list(n_samples = single_cell_manifest$sample_id),
  by = list(
    dataset_id = single_cell_manifest$dataset_id,
    group_std = single_cell_manifest$group_std,
    analysis_role = single_cell_manifest$analysis_role
  ),
  FUN = length
)

write.table(
  group_summary,
  file = project_path("res", "qc", "single_cell", "single_cell_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

conflict_table <- single_cell_manifest[single_cell_manifest$title_metadata_conflict, ]
write.table(
  conflict_table,
  file = project_path("res", "qc", "single_cell", "GSE209781_label_conflict.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Single-cell study manifest written to data/metadata/single_cell_study_manifest.tsv\n")
cat("Single-cell group summary written to res/qc/single_cell/single_cell_group_summary.tsv\n")
cat("GSE209781 label conflict table written to res/qc/single_cell/GSE209781_label_conflict.tsv\n")
