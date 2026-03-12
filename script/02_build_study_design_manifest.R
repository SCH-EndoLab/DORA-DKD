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

read_sample_sheet <- function(dataset_id) {
  path <- project_path("data", "metadata", "soft_samples", paste0(dataset_id, "_samples.tsv"))
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

derive_gse189007_group <- function(title_text) {
  text <- tolower(title_text)

  if (grepl("^control", text)) {
    return("Control")
  }

  if (grepl("without complications", text) && grepl("5 years", text)) {
    return("T2D_NoComp_LT5Y")
  }

  if (grepl("without complications", text) && grepl("15 years", text)) {
    return("T2D_NoComp_GT15Y")
  }

  if (grepl("nephropathy", text) && grepl("5 years", text)) {
    return("T2DN_LT5Y")
  }

  if (grepl("nephropathy", text) && grepl("15 years", text)) {
    return("T2DN_GT15Y")
  }

  if (grepl("retinopathy", text) && grepl("5 years", text)) {
    return("T2DR_LT5Y")
  }

  if (grepl("retinopathy", text) && grepl("15 years", text)) {
    return("T2DR_GT15Y")
  }

  NA_character_
}

derive_gse189007_axis <- function(group_label) {
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

build_manifest_row <- function(
  dataset_id,
  analysis_tier,
  tissue,
  expression_source,
  matrix_file,
  planned_contrast,
  case_definition,
  control_definition,
  total_samples,
  case_samples,
  control_samples,
  status,
  note
) {
  data.frame(
    dataset_id = dataset_id,
    analysis_tier = analysis_tier,
    tissue = tissue,
    expression_source = expression_source,
    matrix_file = matrix_file,
    planned_contrast = planned_contrast,
    case_definition = case_definition,
    control_definition = control_definition,
    total_samples = total_samples,
    case_samples = case_samples,
    control_samples = control_samples,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
}

gse30528 <- read_sample_sheet("GSE30528")
gse30528$group_std <- ifelse(
  gse30528$char_disease_state == "diabetic kidney disease (DKD)",
  "DKD",
  ifelse(gse30528$char_disease_state == "control", "Control", NA_character_)
)

gse96804 <- read_sample_sheet("GSE96804")
gse96804$group_std <- ifelse(
  grepl("diabetic nephropathy", gse96804$source_name_ch1, ignore.case = TRUE),
  "DKD",
  ifelse(
    grepl("tumor nephrectom", gse96804$source_name_ch1, ignore.case = TRUE),
    "Control",
    NA_character_
  )
)

gse104954 <- read_sample_sheet("GSE104954")
gse104954$group_std <- ifelse(
  gse104954$char_diagnosis == "Diabetic nephropathy",
  "DKD",
  ifelse(gse104954$char_diagnosis == "Tumor nephrectomy", "Control", "Other")
)

gse189007 <- read_sample_sheet("GSE189007")
gse189007_clariom <- gse189007[gse189007$platform_id == "GPL23126", , drop = FALSE]
gse189007_clariom$group_std <- vapply(gse189007_clariom$title, derive_gse189007_group, character(1))
gse189007_clariom$disease_axis_group <- vapply(
  gse189007_clariom$group_std,
  derive_gse189007_axis,
  character(1)
)

study_manifest <- rbind(
  build_manifest_row(
    dataset_id = "GSE30528",
    analysis_tier = "primary_discovery",
    tissue = "glomeruli",
    expression_source = "series_matrix",
    matrix_file = "data/raw/series_matrix/GSE30528_series_matrix.txt.gz",
    planned_contrast = "DKD vs Control",
    case_definition = "char_disease_state == diabetic kidney disease (DKD)",
    control_definition = "char_disease_state == control",
    total_samples = nrow(gse30528),
    case_samples = sum(gse30528$group_std == "DKD", na.rm = TRUE),
    control_samples = sum(gse30528$group_std == "Control", na.rm = TRUE),
    status = "expression_ready",
    note = "Microdissected glomeruli; use as the main discovery cohort."
  ),
  build_manifest_row(
    dataset_id = "GSE96804",
    analysis_tier = "glomerular_validation",
    tissue = "glomeruli",
    expression_source = "series_matrix",
    matrix_file = "data/raw/series_matrix/GSE96804_series_matrix.txt.gz",
    planned_contrast = "DKD vs Control",
    case_definition = "source_name_ch1 contains diabetic nephropathy",
    control_definition = "source_name_ch1 contains unaffected portion of tumor nephrectomies",
    total_samples = nrow(gse96804),
    case_samples = sum(gse96804$group_std == "DKD", na.rm = TRUE),
    control_samples = sum(gse96804$group_std == "Control", na.rm = TRUE),
    status = "expression_ready",
    note = "Glomerular validation cohort with larger sample size."
  ),
  build_manifest_row(
    dataset_id = "GSE104954",
    analysis_tier = "tubulointerstitial_validation",
    tissue = "tubulointerstitium",
    expression_source = "series_matrix_split_by_platform",
    matrix_file = paste(
      "data/raw/series_matrix/GSE104954-GPL24120_series_matrix.txt.gz;",
      "data/raw/series_matrix/GSE104954-GPL22945_series_matrix.txt.gz"
    ),
    planned_contrast = "DKD vs Tumor nephrectomy control",
    case_definition = "char_diagnosis == Diabetic nephropathy",
    control_definition = "char_diagnosis == Tumor nephrectomy",
    total_samples = nrow(gse104954),
    case_samples = sum(gse104954$group_std == "DKD", na.rm = TRUE),
    control_samples = sum(gse104954$group_std == "Control", na.rm = TRUE),
    status = "expression_ready",
    note = "Use the full merged matrix, but restrict the primary contrast to DKD and tumor-nephrectomy controls."
  ),
  build_manifest_row(
    dataset_id = "GSE189007",
    analysis_tier = "blood_diabetes_axis_validation",
    tissue = "whole_blood",
    expression_source = "GPL23126_supplementary_or_raw",
    matrix_file = "series_matrix lacks numeric expression; use GPL23126 CEL/CHP supplementary files",
    planned_contrast = "T2D without complication vs T2DN; optional T2DR branch",
    case_definition = "GPL23126 titles containing nephropathy",
    control_definition = "GPL23126 titles containing without complications or Control",
    total_samples = nrow(gse189007_clariom),
    case_samples = sum(gse189007_clariom$disease_axis_group == "T2DN", na.rm = TRUE),
    control_samples = sum(gse189007_clariom$disease_axis_group %in% c("Control", "T2D_NoComp"), na.rm = TRUE),
    status = "metadata_ready_expression_pending",
    note = "Blood diabetes-axis cohort; useful after Clariom D expression extraction."
  )
)

study_counts <- rbind(
  data.frame(dataset_id = "GSE30528", group = sort(unique(gse30528$group_std)), n = as.integer(table(gse30528$group_std)[sort(unique(gse30528$group_std))]), stringsAsFactors = FALSE),
  data.frame(dataset_id = "GSE96804", group = sort(unique(gse96804$group_std)), n = as.integer(table(gse96804$group_std)[sort(unique(gse96804$group_std))]), stringsAsFactors = FALSE),
  data.frame(dataset_id = "GSE104954", group = sort(unique(gse104954$group_std)), n = as.integer(table(gse104954$group_std)[sort(unique(gse104954$group_std))]), stringsAsFactors = FALSE),
  data.frame(
    dataset_id = "GSE189007_GPL23126",
    group = names(sort(table(gse189007_clariom$group_std))),
    n = as.integer(sort(table(gse189007_clariom$group_std))),
    stringsAsFactors = FALSE
  )
)

write_tsv(
  study_manifest,
  project_path("data", "metadata", "study_design_manifest.tsv")
)

write_tsv(
  study_counts,
  project_path("res", "tables", "bulk", "study_group_counts.tsv")
)

cat("Study design manifest written to data/metadata/study_design_manifest.tsv\n")
cat("Group count table written to res/tables/bulk/study_group_counts.tsv\n")
print(study_manifest)
