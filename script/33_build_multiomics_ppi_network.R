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
  library(STRINGdb)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
})

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

candidate_shortlist <- read.delim(
  project_path("res", "qc", "mechanism", "multiomics_candidate_shortlist.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mechanism_membership <- read.delim(
  project_path("data", "metadata", "gene_sets", "mechanism_gene_set_membership.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_genes <- candidate_shortlist$gene_symbol
primary_candidates <- candidate_shortlist$gene_symbol[candidate_shortlist$support_points >= 5]
core_genes <- mechanism_membership$gene_symbol[mechanism_membership$gene_set_id == "OXEIPTOSIS_CORE"]
context_genes <- setdiff(
  mechanism_membership$gene_symbol[mechanism_membership$gene_set_id == "OXEIPTOSIS_EXTENDED"],
  core_genes
)
seed_genes <- unique(c(candidate_genes, core_genes, context_genes))

string_cache_dir <- project_path("ref", "stringdb")
dir.create(string_cache_dir, recursive = TRUE, showWarnings = FALSE)

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = string_cache_dir
)

mapped <- string_db$map(
  data.frame(gene_symbol = seed_genes, stringsAsFactors = FALSE),
  "gene_symbol",
  removeUnmappedRows = TRUE
)

unmapped_genes <- setdiff(seed_genes, mapped$gene_symbol)

seed_interactions <- string_db$get_interactions(mapped$STRING_id)
seed_interactions$key <- apply(
  seed_interactions[, c("from", "to")],
  1,
  function(x) paste(sort(x), collapse = "||")
)
seed_interactions <- seed_interactions[!duplicated(seed_interactions$key), , drop = FALSE]
seed_interactions <- seed_interactions[seed_interactions$combined_score >= 600, , drop = FALSE]

id_to_symbol <- setNames(mapped$gene_symbol, mapped$STRING_id)
edge_df <- data.frame(
  source = id_to_symbol[seed_interactions$from],
  target = id_to_symbol[seed_interactions$to],
  source_string_id = seed_interactions$from,
  target_string_id = seed_interactions$to,
  combined_score = seed_interactions$combined_score,
  stringsAsFactors = FALSE
)
edge_df <- edge_df[!(is.na(edge_df$source) | is.na(edge_df$target) | edge_df$source == edge_df$target), , drop = FALSE]

connected_genes <- unique(c(edge_df$source, edge_df$target))
node_df <- merge(
  mapped[, c("gene_symbol", "STRING_id")],
  candidate_shortlist,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

node_df$node_role <- ifelse(
  node_df$gene_symbol %in% primary_candidates,
  "Candidate_Primary",
  ifelse(
    node_df$gene_symbol %in% candidate_genes,
    "Candidate_Secondary",
    ifelse(
      node_df$gene_symbol %in% core_genes,
      "Mechanism_Core",
      "Mechanism_Context"
    )
  )
)

node_df$node_role <- factor(
  node_df$node_role,
  levels = c("Candidate_Primary", "Candidate_Secondary", "Mechanism_Core", "Mechanism_Context")
)

node_df$node_size <- ifelse(
  node_df$node_role == "Candidate_Primary",
  8.8,
  ifelse(
    node_df$node_role == "Candidate_Secondary",
    7.1,
    ifelse(node_df$node_role == "Mechanism_Core", 6.2, 4.8)
  )
)

node_df$node_fill <- c(
  Candidate_Primary = "#c94c2c",
  Candidate_Secondary = "#f1a340",
  Mechanism_Core = "#2b8cbe",
  Mechanism_Context = "#94a9b7"
)[as.character(node_df$node_role)]

node_df$degree_seed600 <- vapply(
  node_df$gene_symbol,
  function(g) sum(edge_df$source == g | edge_df$target == g),
  integer(1)
)

label_genes <- unique(c(
  primary_candidates,
  core_genes,
  "NFE2L2",
  "HMOX1",
  "NQO1"
))
node_df$label_gene <- node_df$gene_symbol %in% label_genes & node_df$gene_symbol %in% connected_genes
node_df$label_face <- ifelse(node_df$gene_symbol %in% primary_candidates, "bold", "plain")

edge_export <- edge_df[order(-edge_df$combined_score, edge_df$source, edge_df$target), , drop = FALSE]
node_export <- node_df[order(node_df$node_role, -node_df$degree_seed600, node_df$gene_symbol), , drop = FALSE]

dir.create(project_path("res", "tables", "mechanism"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figure", "export", "mechanism"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("res", "qc", "mechanism"), recursive = TRUE, showWarnings = FALSE)

write_tsv(
  node_export,
  project_path("res", "tables", "mechanism", "multiomics_ppi_nodes.tsv")
)
write_tsv(
  edge_export,
  project_path("res", "tables", "mechanism", "multiomics_ppi_edges.tsv")
)
write_tsv(
  data.frame(
    metric = c(
      "n_seed_genes",
      "n_mapped_seed_genes",
      "n_unmapped_seed_genes",
      "n_connected_nodes",
      "n_edges_combined_score_ge_600"
    ),
    value = c(
      length(seed_genes),
      nrow(mapped),
      length(unmapped_genes),
      length(connected_genes),
      nrow(edge_export)
    ),
    stringsAsFactors = FALSE
  ),
  project_path("res", "qc", "mechanism", "multiomics_ppi_summary.tsv")
)
write_tsv(
  data.frame(gene_symbol = unmapped_genes, stringsAsFactors = FALSE),
  project_path("res", "qc", "mechanism", "multiomics_ppi_unmapped_genes.tsv")
)

graph_nodes <- node_export[node_export$gene_symbol %in% connected_genes, , drop = FALSE]
graph_obj <- graph_from_data_frame(
  d = edge_export[, c("source", "target", "combined_score")],
  directed = FALSE,
  vertices = graph_nodes
)

V(graph_obj)$gene_symbol <- V(graph_obj)$name
layout_tbl <- create_layout(as_tbl_graph(graph_obj), layout = "stress")

ppi_plot <- ggraph(layout_tbl) +
  geom_edge_link(
    aes(width = combined_score, alpha = combined_score),
    color = "#7d8a90",
    lineend = "round"
  ) +
  scale_edge_width(range = c(0.35, 1.8), guide = "none") +
  scale_edge_alpha(range = c(0.22, 0.75), guide = "none") +
  geom_node_point(
    aes(size = node_size, fill = node_role),
    shape = 21,
    color = "#1d1d1d",
    stroke = 0.35
  ) +
  scale_size_identity() +
  scale_fill_manual(
    values = c(
      Candidate_Primary = "#c94c2c",
      Candidate_Secondary = "#f1a340",
      Mechanism_Core = "#2b8cbe",
      Mechanism_Context = "#94a9b7"
    )
  ) +
  geom_node_text(
    aes(
      label = ifelse(label_gene, gene_symbol, ""),
      fontface = label_face
    ),
    repel = TRUE,
    size = 3.1,
    family = "sans",
    color = "#111111",
    max.overlaps = Inf,
    point.padding = unit(0.12, "lines"),
    box.padding = unit(0.18, "lines")
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Multi-omics candidate-centered PPI network",
    subtitle = "STRING v12.0, human, combined score >= 600",
    fill = "Node role"
  )

ggsave(
  filename = project_path("figure", "export", "mechanism", "multiomics_candidate_ppi_network.pdf"),
  plot = ppi_plot,
  width = 10.2,
  height = 8.2
)

ggsave(
  filename = project_path("figure", "export", "mechanism", "multiomics_candidate_ppi_network.png"),
  plot = ppi_plot,
  width = 10.2,
  height = 8.2,
  dpi = 320
)

cat("Multi-omics PPI node table written to res/tables/mechanism/multiomics_ppi_nodes.tsv\n")
cat("Multi-omics PPI edge table written to res/tables/mechanism/multiomics_ppi_edges.tsv\n")
cat("Multi-omics PPI summary written to res/qc/mechanism/multiomics_ppi_summary.tsv\n")
cat("Multi-omics PPI unmapped genes written to res/qc/mechanism/multiomics_ppi_unmapped_genes.tsv\n")
cat("Multi-omics PPI plot written to figure/export/mechanism/multiomics_candidate_ppi_network.pdf\n")
