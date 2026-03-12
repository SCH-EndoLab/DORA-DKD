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

marker_sets <- list(
  PT = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "ALDOB", "PDZK1"),
  PT_INJURED = c("HAVCR1", "LCN2", "VCAM1", "KRT8", "KRT18", "SPP1"),
  TAL = c("UMOD", "SLC12A1", "KCNJ1", "CLDN16"),
  DISTAL_CD = c("SLC12A3", "CALB1", "AQP2", "FXYD4", "ATP6V1B1", "SLC4A1"),
  PODOCYTE = c("NPHS1", "NPHS2", "PODXL", "SYNPO", "WT1"),
  ENDOTHELIAL = c("EMCN", "KDR", "PECAM1", "ESAM", "FLT1", "VWF", "EGFL7"),
  MESANGIAL_PERICYTE = c("PDGFRB", "RGS5", "MCAM", "CSPG4", "DES", "COL4A1", "NOTCH3"),
  FIBROBLAST = c("DCN", "COL1A1", "COL1A2", "LUM", "PDGFRA"),
  MACROPHAGE = c("LYZ", "FCER1G", "TYROBP", "CTSS", "C1QC", "AIF1"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "NKG7", "KLRD1", "IL7R"),
  B_PLASMA = c("MS4A1", "CD79A", "MZB1", "JCHAIN", "SDC1"),
  PROLIFERATION = c("MKI67", "TOP2A", "UBE2C", "CENPF")
)

marker_df <- do.call(
  rbind,
  lapply(names(marker_sets), function(set_id) {
    data.frame(
      marker_set_id = set_id,
      gene_symbol = marker_sets[[set_id]],
      stringsAsFactors = FALSE
    )
  })
)

dir.create(project_path("data", "metadata", "gene_sets"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("manuscript", "notes"), recursive = TRUE, showWarnings = FALSE)

write.table(
  marker_df,
  file = project_path("data", "metadata", "gene_sets", "kidney_celltype_marker_membership.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_df <- aggregate(
  rep(1L, nrow(marker_df)),
  by = list(marker_set_id = marker_df$marker_set_id),
  FUN = sum
)
names(summary_df)[2] <- "n_markers"

write.table(
  summary_df,
  file = project_path("res", "qc", "single_cell", "kidney_celltype_marker_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

note_lines <- c(
  "# Kidney Celltype Marker Note",
  "",
  "This marker panel is intended for coarse cell-state and cell-type inference in kidney and urine single-cell matrices.",
  "",
  "Included marker sets:",
  "",
  paste0("- `", summary_df$marker_set_id, "`: ", summary_df$n_markers, " markers")
)

writeLines(
  note_lines,
  con = project_path("manuscript", "notes", "kidney_celltype_marker_note.md")
)

cat("Kidney celltype marker membership written to data/metadata/gene_sets/kidney_celltype_marker_membership.tsv\n")
cat("Kidney celltype marker note written to manuscript/notes/kidney_celltype_marker_note.md\n")
