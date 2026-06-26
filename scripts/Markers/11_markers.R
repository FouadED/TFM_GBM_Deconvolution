# ==============================================================================
# Cell-Type Marker Identification (one-vs-rest on BayesPrism-purified expression)
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 11_Marcadores (computation)
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Identify cell-type marker genes from the BayesPrism-purified expression
#   (the Z array of the initial cell-type posterior), which removes cross-cell-
#   type contamination. Two comparisons are run:
#     - one-vs-rest for non-myeloid cell types (each type vs all others);
#     - within-myeloid for the myeloid subtypes (each vs the other myeloid
#       subtypes), since the broad one-vs-rest cannot separate similar myeloid
#       programs.
#   All markers therefore come from the SAME purified source as the stemness
#   score (resolving the marker-source ambiguity).
#
# LOG-FOLD-CHANGE AND THE PSEUDO-COUNT
#   LFC = log2((mean_target + PSEUDO) / (mean_rest + PSEUDO)).
#   A near-zero PSEUDO inflates the LFC when a gene is essentially absent from
#   the 'rest' group (mean_rest ~ 0): dividing by ~0 produces an artefactually
#   huge LFC (this is what made ADAMTS20 reach LFC~14). PSEUDO is therefore set
#   to the detection threshold MIN_EXPR, so a gene undetected in the 'rest'
#   group cannot drive an inflated ratio, while genuine markers (clearly above
#   the floor in the target) keep a realistic LFC. mean_target and mean_rest are
#   written to the output so any remaining low-expression cases are transparent.
#
# INPUT
#   scripts/02_BayesPrism/results/bayesprism_full_result.rds
#
# OUTPUTS (canonical marker location, shared with 11_Marcadores figures and 13)
#   scripts/11_Marcadores/results/markers_<type>.csv
#   scripts/11_Marcadores/results/markers_withinMyeloid_<type>.csv
#   scripts/11_Marcadores/results/markers_summary_top10.csv
#   scripts/11_Marcadores/results/purified_expression/exp_purified_<type>.csv
#
# REFERENCE
#   Chu et al. (2022) BayesPrism. Nat Cancer 3:505-517.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BP_RDS      <- "scripts/02_BayesPrism/results/bayesprism_full_result.rds"
OUT_DIR     <- "scripts/11_Marcadores/results"
EXP_DIR     <- file.path(OUT_DIR, "purified_expression")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EXP_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N       <- 50    # top markers reported per cell type
LFC_GLOBAL  <- 1.0   # one-vs-rest LFC threshold
LFC_MYELOID <- 0.5   # within-myeloid LFC threshold
MIN_EXPR    <- 0.5   # minimum mean expression in the target type
PADJ_CUT    <- 0.05  # BH-adjusted p-value cutoff

# Pseudo-count for the LFC: set to the detection threshold so that genes absent
# from the 'rest' group cannot produce inflated fold changes (see header).
PSEUDO      <- MIN_EXPR

MYELOID_TYPES <- c("a-microglia", "AP-microglia", "h-microglia", "i-microglia",
                   "MDSC", "Myeloid", "s-mac 1", "s-mac 2")

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(BayesPrism)
})

# ==============================================================================
# HELPER: one-vs-rest markers for a single cell type
#   mat_target, mat_rest: samples x genes. Returns a data frame of significant
#   markers (LFC > threshold, mean_target > MIN_EXPR, padj < PADJ_CUT), top N.
# ==============================================================================
get_markers <- function(mat_target, mat_rest, genes, lfc_thresh, top_n, label) {
  mean_t <- colMeans(mat_target)
  mean_r <- colMeans(mat_rest)
  lfc    <- log2((mean_t + PSEUDO) / (mean_r + PSEUDO))

  candidates <- genes[lfc > lfc_thresh & mean_t > MIN_EXPR]
  if (length(candidates) == 0) return(NULL)

  pvals <- vapply(candidates, function(g) {
    tryCatch(wilcox.test(mat_target[, g], mat_rest[, g],
                         alternative = "greater")$p.value,
             error = function(e) NA_real_)
  }, numeric(1))
  padj <- p.adjust(pvals, method = "BH")

  df <- data.frame(
    gene        = candidates,
    mean_target = mean_t[candidates],
    mean_rest   = mean_r[candidates],
    log2FC      = lfc[candidates],
    pvalue      = pvals,
    padj        = padj,
    comparison  = label,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$padj) & df$padj < PADJ_CUT, ]
  head(df[order(df$log2FC, decreasing = TRUE), ], top_n)
}

# ==============================================================================
# 1. LOAD BAYESPRISM RESULT AND THE PURIFIED Z ARRAY
# ==============================================================================
message("\n[1] Loading BayesPrism result...")
bp <- readRDS(BP_RDS)
Z  <- bp@posterior.initial.cellType@Z      # samples x genes x cell types
genes      <- dimnames(Z)[[2]]
cell_types <- dimnames(Z)[[3]]
message(sprintf("  Z array: %s (samples x genes x cell types)",
                paste(dim(Z), collapse = " x ")))
message(sprintf("  Cell types: %s", paste(cell_types, collapse = ", ")))

# ==============================================================================
# 2. EXPORT PURIFIED EXPRESSION PER CELL TYPE
# ==============================================================================
message("\n[2] Exporting purified expression per cell type...")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  write.csv(Z[, , ct], file.path(EXP_DIR, paste0("exp_purified_", ct_safe, ".csv")))
}
message(sprintf("  Wrote %d purified expression files", length(cell_types)))

# ==============================================================================
# 3. ONE-VS-REST MARKERS (non-myeloid cell types)
# ==============================================================================
message("\n[3] One-vs-rest markers (non-myeloid)...")
non_myeloid    <- setdiff(cell_types, MYELOID_TYPES)
markers_global <- list()

for (ct in non_myeloid) {
  rest  <- setdiff(cell_types, ct)
  mat_r <- apply(Z[, , rest], c(1, 2), mean)
  df    <- get_markers(Z[, , ct], mat_r, genes, LFC_GLOBAL, TOP_N, "one-vs-rest")
  if (!is.null(df) && nrow(df) > 0) {
    markers_global[[ct]] <- df
    ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
    write.csv(df, file.path(OUT_DIR, paste0("markers_", ct_safe, ".csv")),
              row.names = FALSE)
    message(sprintf("  %-16s: %d markers (top LFC=%.1f, %s)",
                    ct, nrow(df), df$log2FC[1], df$gene[1]))
  } else {
    message(sprintf("  %-16s: no significant markers", ct))
  }
}

# ==============================================================================
# 4. WITHIN-MYELOID MARKERS
# ==============================================================================
message("\n[4] Within-myeloid markers...")
myeloid_ok     <- intersect(MYELOID_TYPES, cell_types)
markers_myel   <- list()

for (ct in myeloid_ok) {
  rest_m <- setdiff(myeloid_ok, ct)
  mat_r  <- apply(Z[, , rest_m], c(1, 2), mean)
  df     <- get_markers(Z[, , ct], mat_r, genes, LFC_MYELOID, TOP_N, "within-myeloid")
  if (!is.null(df) && nrow(df) > 0) {
    markers_myel[[ct]] <- df
    ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
    write.csv(df, file.path(OUT_DIR, paste0("markers_withinMyeloid_", ct_safe, ".csv")),
              row.names = FALSE)
    message(sprintf("  %-16s: %d markers (top LFC=%.1f, %s)",
                    ct, nrow(df), df$log2FC[1], df$gene[1]))
  } else {
    message(sprintf("  %-16s: no significant markers", ct))
  }
}

# ==============================================================================
# 5. SUMMARY TABLE (top 10 per characterised type)
# ==============================================================================
message("\n[5] Writing summary table...")
all_markers <- c(markers_global, markers_myel)
if (length(all_markers) > 0) {
  summary_tbl <- do.call(rbind, lapply(names(all_markers), function(ct) {
    df <- head(all_markers[[ct]], 10)
    if (nrow(df) == 0) return(NULL)
    data.frame(cell_type = ct, df, stringsAsFactors = FALSE)
  }))
  write.csv(summary_tbl, file.path(OUT_DIR, "markers_summary_top10.csv"),
            row.names = FALSE)
  message("  Wrote markers_summary_top10.csv")
}

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 60))
message("  SCRIPT 11 (marker computation) COMPLETED")
message(strrep("=", 60))
message(sprintf("  One-vs-rest  : %d / %d types with markers",
                length(markers_global), length(non_myeloid)))
message(sprintf("  Within-myeloid: %d / %d types with markers",
                length(markers_myel), length(myeloid_ok)))
message(sprintf("  Pseudo-count : %.3f (= MIN_EXPR; prevents inflated LFCs)", PSEUDO))
message(sprintf("  Markers in   : %s", OUT_DIR))
message("  NOTE: figures are produced separately by the 11_Marcadores figures script.")
message(strrep("=", 60))
