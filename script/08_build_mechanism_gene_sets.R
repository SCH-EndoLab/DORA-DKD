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

gene_set_records <- rbind(
  data.frame(
    gene_set_id = "OXEIPTOSIS_CORE",
    gene_symbol = c("KEAP1", "PGAM5", "AIFM1"),
    evidence_tier = "core",
    evidence_note = c(
      "ROS sensor linked to oxeiptosis execution",
      "Phosphatase released from KEAP1 under high ROS",
      "Execution factor dephosphorylated by PGAM5"
    ),
    source_ref = "PMID:29255269",
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene_set_id = "KEAP1_NRF2_RESPONSE",
    gene_symbol = c(
      "KEAP1", "NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM",
      "TXN", "TXNRD1", "SOD1", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2"
    ),
    evidence_tier = c("context", "context", rep("context", 12)),
    evidence_note = c(
      "Sensor node shared by antioxidant buffering and oxeiptosis thresholding",
      "Canonical oxidative stress response transcription factor downstream of KEAP1",
      rep("Canonical oxidative stress response effector", 12)
    ),
    source_ref = c(
      "PMID:29255269",
      "PMID:29255269",
      rep("Canonical KEAP1-NRF2 oxidative stress response context", 12)
    ),
    stringsAsFactors = FALSE
  )
)

extended_genes <- unique(gene_set_records$gene_symbol)
extended_set <- data.frame(
  gene_set_id = "OXEIPTOSIS_EXTENDED",
  gene_symbol = extended_genes,
  evidence_tier = ifelse(
    extended_genes %in% c("KEAP1", "PGAM5", "AIFM1"),
    "core_or_context_union",
    "context_union"
  ),
  evidence_note = "Union set for transcriptome-level scoring in kidney bulk cohorts",
  source_ref = "PMID:29255269 + KEAP1/NRF2 oxidative stress context",
  stringsAsFactors = FALSE
)

all_sets <- rbind(gene_set_records, extended_set)

set_summary <- aggregate(
  gene_symbol ~ gene_set_id,
  data = all_sets,
  FUN = length
)
names(set_summary)[2] <- "n_genes"

write_tsv(
  all_sets,
  project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv")
)

write_tsv(
  set_summary,
  project_path("res", "qc", "mechanism", "mechanism_gene_set_summary.tsv")
)

cat("Mechanism gene set membership written to data/metadata/gene_sets/mechanism_gene_set_membership.tsv\n")
cat("Mechanism gene set summary written to res/qc/mechanism/mechanism_gene_set_summary.tsv\n")
print(set_summary)
