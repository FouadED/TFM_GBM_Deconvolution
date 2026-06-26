# ==============================================================================
# PPI Networks, GO Enrichment and Cytoscape Export
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 13_Networks
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   For each cell type with one-vs-rest markers, run three complementary
#   analyses:
#     A) PPI network (STRINGdb v12.0 + igraph): hub genes by degree and
#        betweenness; hub_score = degree * log2FC (connectivity x specificity).
#     B) GO Biological Process enrichment (clusterProfiler) over the marker set.
#     C) Cytoscape export: one GraphML per cell type with the network structure
#        AND all node attributes embedded (degree, log2FC, hub_score, node_type),
#        avoiding the key-matching issues of a two-file CSV import.
#
# HUB DEFINITION
#   A high-LFC marker is specific but not necessarily important. Hubs are defined
#   by network degree (top quintile, minimum 2). hub_score combines both axes
#   (degree * log2FC) so a gene ranks highly only when it is both specific to the
#   cell type AND well connected. Degree and hub_score are reported separately so
#   it is clear which axis drives each hub. All hubs are whatever the network
#   analysis returns from the input markers; no gene is added or removed by hand.
#
# MARKER SOURCE
#   Markers are read from MARKERS_DIR; this must be the same directory the
#   marker-calling script (11) and the 11 figures script point to.
#
# REQUIRES INTERNET: STRINGdb downloads interaction data on first run; run
#   locally, not on an offline HPC node.
#
# INPUTS
#   <MARKERS_DIR>/markers_<type>.csv   (columns: gene, log2FC, padj)
#
# OUTPUTS
#   scripts/13_Networks/results/hubs_<type>.csv
#   scripts/13_Networks/results/all_hubs_summary.csv
#   scripts/13_Networks/results/GO/GO_BP_<type>.csv
#   scripts/13_Networks/results/cytoscape/<type>.graphml
#   scripts/13_Networks/figures/network_<type>.png
#   scripts/13_Networks/figures/GO/GO_dotplot_<type>.png
#
# REFERENCES
#   Szklarczyk et al. (2023) STRING. Nucleic Acids Res 51:D638.
#   Csardi & Nepusz (2006) igraph. InterJournal Complex Systems 1695.
#   Wu et al. (2021) clusterProfiler 4.0. Innovation 2:100141.
# ==============================================================================
set.seed(123)
# ==============================================================================
# DEPENDENCIES
# ==============================================================================
suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})
source("scripts/utils/tfm_theme.R")  # tfm_colors, theme_tfm
# clusterProfiler masks dplyr::select/filter; reassign explicitly.
select <- dplyr::select
filter <- dplyr::filter
# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Canonical marker directory: written by 11_markers.R, also read by the
# 11_Marcadores figures script. All three must point here.
MARKERS_DIR <- "scripts/11_Marcadores/results"
OUT_DIR     <- "scripts/13_Networks/results"
FIG_DIR     <- "scripts/13_Networks/figures"
GO_OUT_DIR  <- file.path(OUT_DIR, "GO")
GO_FIG_DIR  <- file.path(FIG_DIR, "GO")
CYTO_DIR    <- file.path(OUT_DIR, "cytoscape")
for (d in c(OUT_DIR, FIG_DIR, GO_OUT_DIR, GO_FIG_DIR, CYTO_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
CELL_TYPES <- c("Glioma", "Oligo", "Endo", "CD8_TCells", "TCells",
                "a_microglia", "Pericytes", "Proliferating",
                "s_mac_1", "s_mac_2")
STRING_SCORE_THRESHOLD <- 400   # medium confidence
HUB_DEGREE_QUANTILE    <- 0.80  # hubs = top quintile by degree (min 2)
MIN_MARKERS            <- 5     # minimum markers to attempt analysis
# Network node colours (igraph base plot).
HUB_COLOR        <- tfm_colors[["Glioma"]]
PERIPHERAL_COLOR <- "#A9C6DA"
# ==============================================================================
# INITIALISE STRING
# ==============================================================================
message(strrep("=", 60))
message(" Networks + GO enrichment + Cytoscape GraphML export")
message(strrep("=", 60))
message("Initialising STRINGdb v12.0 (downloads on first run)...")
string_db <- STRINGdb$new(
  version         = "12.0",
  species         = 9606,
  score_threshold = STRING_SCORE_THRESHOLD,
  network_type    = "full",
  input_directory = OUT_DIR
)
message("STRINGdb ready.")
# ==============================================================================
# PART A: PPI NETWORK ANALYSIS
# ==============================================================================
message("\n--- PART A: PPI network analysis ---")
analyze_network <- function(ct) {
  f <- file.path(MARKERS_DIR, paste0("markers_", ct, ".csv"))
  if (!file.exists(f)) { message(sprintf("  [skip] %s -- file not found", ct)); return(NULL) }
  mk <- read.csv(f)
  if (nrow(mk) < MIN_MARKERS) {
    message(sprintf("  [skip] %s -- only %d markers", ct, nrow(mk))); return(NULL)
  }
  message(sprintf("\n  --- %s (%d markers) ---", ct, nrow(mk)))
  mapped <- tryCatch(string_db$map(mk, "gene", removeUnmappedRows = TRUE),
                     error = function(e) NULL)
  if (is.null(mapped) || nrow(mapped) < MIN_MARKERS) {
    message("    [skip] too few genes mapped"); return(NULL)
  }
  message(sprintf("    Mapped: %d / %d", nrow(mapped), nrow(mk)))
  g <- tryCatch(string_db$get_subnetwork(mapped$STRING_id), error = function(e) NULL)
  if (is.null(g) || vcount(g) < MIN_MARKERS) {
    message("    [skip] network too small"); return(NULL)
  }
  message(sprintf("    Network: %d nodes | %d edges", vcount(g), ecount(g)))
  id2sym <- setNames(mapped$gene, mapped$STRING_id)
  cent <- data.frame(
    STRING_id   = names(degree(g)),
    gene        = id2sym[names(degree(g))],
    degree      = as.integer(degree(g)),
    betweenness = round(betweenness(g, normalized = TRUE), 4),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(mapped %>% dplyr::select(STRING_id, log2FC, padj), by = "STRING_id") %>%
    dplyr::arrange(desc(degree), desc(betweenness))
  thr  <- as.integer(max(quantile(cent$degree, HUB_DEGREE_QUANTILE), 2))
  hubs <- cent %>%
    dplyr::filter(degree >= thr) %>%
    dplyr::mutate(hub_score = degree * log2FC, cell_type = ct,
                  n_nodes = vcount(g), n_edges = ecount(g)) %>%
    dplyr::arrange(desc(hub_score))
  message(sprintf("    Hubs (degree >= %d): %d", thr, nrow(hubs)))
  if (nrow(hubs) > 0) print(hubs %>% dplyr::select(gene, degree, log2FC, hub_score) %>% head(5))
  write.csv(hubs, file.path(OUT_DIR, paste0("hubs_", ct, ".csv")), row.names = FALSE)
  # Network figure
  V(g)$symbol <- id2sym[V(g)$name]
  V(g)$is_hub <- V(g)$name %in% hubs$STRING_id
  V(g)$color  <- ifelse(V(g)$is_hub, HUB_COLOR, PERIPHERAL_COLOR)
  V(g)$size   <- 4 + (degree(g) / max(max(degree(g)), 1)) * 14
  V(g)$label  <- ifelse(V(g)$is_hub, V(g)$symbol, NA)
  png(file.path(FIG_DIR, paste0("network_", ct, ".png")),
      width = 2200, height = 1900, res = 220)
  set.seed(123)
  lay <- layout_with_fr(g)
  plot(g, layout = lay, vertex.label = V(g)$label, vertex.label.cex = 0.72,
       vertex.label.color = "black", vertex.label.font = 2,
       vertex.frame.color = "white", edge.color = "grey75", edge.width = 0.9,
       main = sprintf("%s  |  %d nodes  |  %d edges  |  %d hubs (red)",
                      ct, vcount(g), ecount(g), nrow(hubs)))
  dev.off()
  list(hubs = hubs, graph = g, mapped = mapped, id2sym = id2sym, thr = thr)
}
network_results <- lapply(CELL_TYPES, function(ct) {
  tryCatch(analyze_network(ct),
           error = function(e) { message(sprintf("  ERROR in %s: %s", ct, e$message)); NULL })
})
names(network_results) <- CELL_TYPES
all_hubs <- dplyr::bind_rows(lapply(network_results, function(r) r$hubs))
write.csv(all_hubs, file.path(OUT_DIR, "all_hubs_summary.csv"), row.names = FALSE)
message("\n=== Network summary ===")
message(sprintf("Total hub genes: %d", nrow(all_hubs)))
if (nrow(all_hubs) > 0) {
  print(all_hubs %>% dplyr::group_by(cell_type) %>%
          dplyr::summarise(n_hubs = dplyr::n(),
                           top_hub = gene[which.max(hub_score)],
                           max_degree = max(degree), .groups = "drop"))
}
# ==============================================================================
# PART B: GO BIOLOGICAL PROCESS ENRICHMENT
# ==============================================================================
message("\n--- PART B: GO Biological Process enrichment ---")
run_GO <- function(ct) {
  f <- file.path(MARKERS_DIR, paste0("markers_", ct, ".csv"))
  if (!file.exists(f)) { message(sprintf("  [skip] %s -- not found", ct)); return(NULL) }
  mk <- read.csv(f)
  if (nrow(mk) < MIN_MARKERS) {
    message(sprintf("  [skip] %s -- only %d markers", ct, nrow(mk))); return(NULL)
  }
  message(sprintf("\n  --- GO: %s (%d genes) ---", ct, nrow(mk)))
  conv <- tryCatch(bitr(mk$gene, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db), error = function(e) NULL)
  if (is.null(conv) || nrow(conv) < MIN_MARKERS) {
    message("    [skip] too few Entrez IDs"); return(NULL)
  }
  message(sprintf("    Entrez mapped: %d / %d", nrow(conv), nrow(mk)))
  ego <- tryCatch(enrichGO(gene = conv$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                           pAdjustMethod = "BH", pvalueCutoff = 0.05,
                           qvalueCutoff = 0.10, readable = TRUE),
                  error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("    [none] no significant GO terms"); return(NULL)
  }
  n_terms <- nrow(as.data.frame(ego))
  message(sprintf("    Significant GO BP terms: %d", n_terms))
  write.csv(as.data.frame(ego), file.path(GO_OUT_DIR, paste0("GO_BP_", ct, ".csv")),
            row.names = FALSE)
  print(as.data.frame(ego) %>% dplyr::select(Description, GeneRatio, p.adjust) %>% head(5))
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(sprintf("GO Biological Process -- %s", ct)) +
    labs(subtitle = sprintf("%d marker genes | BH correction | %d significant terms",
                            nrow(mk), n_terms),
         caption = "Markers: BayesPrism-purified expression | one-vs-rest (LFC>1, padj<0.05)") +
    theme_tfm(show_x_text = TRUE)
  ggsave(file.path(GO_FIG_DIR, paste0("GO_dotplot_", ct, ".png")),
         p, width = 10, height = 7, dpi = 300, bg = "white")
  message(sprintf("    Dotplot saved: GO_dotplot_%s.png", ct))
  invisible(ego)
}
go_results <- lapply(CELL_TYPES, function(ct) {
  tryCatch(run_GO(ct),
           error = function(e) { message(sprintf("  ERROR in %s: %s", ct, e$message)); NULL })
})
names(go_results) <- CELL_TYPES
# ==============================================================================
# PART C: CYTOSCAPE EXPORT (single-file GraphML with embedded attributes)
# ==============================================================================
message("\n--- PART C: Cytoscape GraphML export ---")
export_graphml <- function(ct, net_result) {
  if (is.null(net_result)) return(invisible(NULL))
  g <- net_result$graph; hubs <- net_result$hubs
  mapped <- net_result$mapped; id2sym <- net_result$id2sym
  sym_names <- id2sym[V(g)$name]
  lfc_map   <- setNames(mapped$log2FC,  mapped$gene)
  score_map <- setNames(hubs$hub_score, hubs$gene)
  V(g)$label     <- sym_names
  V(g)$degree    <- as.integer(degree(g))
  V(g)$is_hub    <- as.integer(V(g)$name %in% hubs$STRING_id)
  V(g)$node_type <- ifelse(V(g)$name %in% hubs$STRING_id, "hub", "peripheral")
  V(g)$log2FC    <- round(as.numeric(lfc_map[sym_names]), 4)
  V(g)$hub_score <- round(as.numeric(score_map[sym_names]), 4)
  V(g)$hub_score[is.na(V(g)$hub_score)] <- 0
  V(g)$cell_type <- ct
  V(g)$name      <- sym_names   # gene symbols as node IDs for Cytoscape
  out_file <- file.path(CYTO_DIR, paste0(ct, ".graphml"))
  write_graph(g, out_file, format = "graphml")
  message(sprintf("  %-14s -> %s  (%d nodes | %d edges | %d hubs)",
                  ct, basename(out_file), vcount(g), ecount(g), sum(V(g)$is_hub)))
}
for (ct in CELL_TYPES) {
  tryCatch(export_graphml(ct, network_results[[ct]]),
           error = function(e) message(sprintf("  ERROR %s: %s", ct, e$message)))
}
graphml_files <- list.files(CYTO_DIR, pattern = "\\.graphml$")
message(sprintf("  %d GraphML files in %s", length(graphml_files), CYTO_DIR))
# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
message("\n", strrep("=", 60))
message(" FINAL SUMMARY")
message(strrep("=", 60))
message("\nPART A -- PPI networks:")
for (ct in CELL_TYPES) {
  r <- network_results[[ct]]
  if (!is.null(r) && nrow(r$hubs) > 0) {
    message(sprintf("  %-16s: %d hubs | top=%s (deg=%d, LFC=%.1f) | %d edges",
                    ct, nrow(r$hubs), r$hubs$gene[1], r$hubs$degree[1],
                    r$hubs$log2FC[1], r$hubs$n_edges[1]))
  } else {
    message(sprintf("  %-16s: skipped or no hubs", ct))
  }
}
message("\nPART B -- GO BP enrichment:")
for (ct in CELL_TYPES) {
  r <- go_results[[ct]]
  if (!is.null(r)) {
    message(sprintf("  %-16s: %d terms | top=%s", ct, nrow(as.data.frame(r)),
                    substr(as.data.frame(r)$Description[1], 1, 50)))
  } else {
    message(sprintf("  %-16s: no significant terms", ct))
  }
}
message(sprintf("\nPART C -- Cytoscape GraphML: %d files in %s",
                length(graphml_files), CYTO_DIR))
message("  Import: File -> Import -> Network from File -> <type>.graphml")
message(strrep("=", 60))