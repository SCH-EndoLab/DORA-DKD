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

suppressPackageStartupMessages(library(limma))

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

beta_df <- read_tsv_auto(project_path("data", "processed", "methylation_GSE233758_beta.tsv.gz"))
pheno <- read_tsv_auto(project_path("data", "processed", "GSE233758_pheno.tsv"))
ann <- read_tsv_auto(project_path("data", "processed", "GPL21145_annotation.tsv.gz"))

sample_ids <- pheno$sample_id
beta_df <- beta_df[, c("cpg_id", sample_ids), drop = FALSE]
beta_mat <- as.matrix(beta_df[, sample_ids, drop = FALSE])
storage.mode(beta_mat) <- "double"
rownames(beta_mat) <- beta_df$cpg_id

beta_mat[is.na(beta_mat)] <- 0.5
beta_mat[beta_mat <= 1e-6] <- 1e-6
beta_mat[beta_mat >= 1 - 1e-6] <- 1 - 1e-6
m_mat <- log2(beta_mat / (1 - beta_mat))

group <- factor(pheno$group_std, levels = c("Control", "DKD"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(m_mat, design)
contrast_matrix <- makeContrasts(DKD - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

result <- topTable(fit2, number = Inf, sort.by = "P")
result$cpg_id <- rownames(result)
control_beta <- rowMeans(beta_mat[, group == "Control", drop = FALSE], na.rm = TRUE)
dkd_beta <- rowMeans(beta_mat[, group == "DKD", drop = FALSE], na.rm = TRUE)
result$mean_beta_control <- control_beta[result$cpg_id]
result$mean_beta_dkd <- dkd_beta[result$cpg_id]
result$delta_beta_dkd_minus_control <- result$mean_beta_dkd - result$mean_beta_control
result$dm_direction <- ifelse(result$delta_beta_dkd_minus_control > 0, "Hypermethylated_in_DKD", "Hypomethylated_in_DKD")
result$sig_fdr <- ifelse(result$adj.P.Val < 0.05, "yes", "no")
result$sig_fdr_dbeta <- ifelse(result$adj.P.Val < 0.05 & abs(result$delta_beta_dkd_minus_control) >= 0.10, "yes", "no")

result <- merge(result, ann, by = "cpg_id", all.x = TRUE, sort = FALSE)
result <- result[, c(
  "cpg_id", "gene_symbol_primary", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
  "CHR", "MAPINFO", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group",
  "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
  "mean_beta_control", "mean_beta_dkd", "delta_beta_dkd_minus_control",
  "dm_direction", "sig_fdr", "sig_fdr_dbeta"
)]

write_tsv_gz(
  result,
  project_path("res", "tables", "mechanism", "GSE233758_methylation_limma.tsv.gz")
)

summary_table <- data.frame(
  n_tested = nrow(result),
  n_sig_fdr = sum(result$sig_fdr == "yes"),
  n_sig_fdr_dbeta = sum(result$sig_fdr_dbeta == "yes"),
  top_cpg = result$cpg_id[[1]],
  top_gene = result$gene_symbol_primary[[1]],
  top_adj_p = result$adj.P.Val[[1]],
  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  project_path("res", "qc", "mechanism", "GSE233758_methylation_limma_summary.tsv")
)

cat("GSE233758 methylation limma table written to res/tables/mechanism/GSE233758_methylation_limma.tsv.gz\n")
print(summary_table)
