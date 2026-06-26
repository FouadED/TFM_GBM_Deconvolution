# ==============================================================================
# Immune/Stromal Scoring with xCell 2.0 (ROL2, supplementary)
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 04_xCell
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Compute enrichment scores for immune and stromal cell types in the bulk
#   cohort using xCell 2.0. This is a reference-free (ROL2) method used only
#   for SUPPLEMENTARY, directional validation of immune trends; its scores are
#   on an arbitrary enrichment scale and are NOT combined with the ROL1
#   consensus proportions.
#
# REFERENCE PANEL
#   BlueprintEncode is used as the cell-type signature reference: it is the
#   recommended panel for adult human tumour microenvironment (Angel et al.,
#   2025), covering immune and stromal populations with external validation.
#
# INPUT NORMALISATION
#   xCell 2.0 expects TPM. The bulk TPM is obtained, in order of preference:
#     1. reused from the ESTIMATE step if available; otherwise
#     2. recomputed from raw counts and the EPIC gene-length cache.
#   (xCell 2.0 does NOT use the LGG-04 scRNA-seq reference, so it is unaffected
#    by the raw-count update applied to the reference object.)
#
# INPUTS
#   data/bulk/oligo_counts_protein.coding.xlsx          (raw counts)
#   scripts/06_ESTIMATE/results/bulk_tpm_for_estimate.txt (TPM, if present)
#   scripts/03_EPIC/results/gene_lengths_cache.csv        (gene lengths)
#
# OUTPUT
#   scripts/04_xCell/results/xcell2_scores.csv  (samples x cell types)
#
# NOTE
#   xCell 2.0 is installed from GitHub on first run; this requires internet
#   access (run on a machine with connectivity, not on an offline login node).
#
# REFERENCE
#   Angel et al. (2025) xCell 2.0. Genome Biol 26:335.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
COUNTS_PATH   <- "data/bulk/oligo_counts_protein.coding.xlsx"
TPM_EST_PATH  <- "scripts/06_ESTIMATE/results/bulk_tpm_for_estimate.txt"
LENGTHS_PATH  <- "scripts/03_EPIC/results/gene_lengths_cache.csv"
OUT_DIR       <- "scripts/04_xCell/results"

MIN_OVERLAP_GENES <- 5000   # warn below this bulk/reference gene overlap
MIN_SHARED_FRAC   <- 0.85   # xCell2 minSharedGenes
SPILLOVER_ALPHA   <- 0.5    # xCell2 spillover correction strength

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# ==============================================================================
# 1. INSTALL xCell 2.0 (first run only; requires internet)
# ==============================================================================
message("\n[1/5] Checking xCell 2.0 installation...")
if (!requireNamespace("xCell2", quietly = TRUE)) {
  message("  xCell2 not found. Installing from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
  tryCatch(
    devtools::install_github("AlmogAngel/xCell2", upgrade = "never"),
    error = function(e) {
      stop("xCell2 installation failed: ", conditionMessage(e),
           "\nInstall manually with devtools::install_github('AlmogAngel/xCell2') ",
           "and re-run.", call. = FALSE)
    }
  )
}
suppressPackageStartupMessages(library(xCell2))
message("  xCell 2.0 ready.")

# ==============================================================================
# HELPER: TPM from raw counts and gene lengths
#   The bp-vs-kb scaling cancels in the per-sample normalisation, so lengths
#   may be supplied in base pairs.
# ==============================================================================
compute_tpm <- function(counts, gene_lengths_bp) {
  shared       <- intersect(rownames(counts), names(gene_lengths_bp))
  counts_f     <- counts[shared, , drop = FALSE]
  rpk          <- counts_f / gene_lengths_bp[shared]
  t(t(rpk) / colSums(rpk) * 1e6)
}

# ==============================================================================
# 2. OBTAIN BULK TPM
# ==============================================================================
message("\n[2/5] Preparing bulk TPM...")
bulk_tpm   <- NULL
tpm_source <- NA_character_

# -- Option 1: reuse ESTIMATE's TPM -------------------------------------------
if (file.exists(TPM_EST_PATH)) {
  tpm_tab <- read.table(TPM_EST_PATH, header = TRUE, row.names = 1,
                        sep = "\t", check.names = FALSE)
  if (all(vapply(tpm_tab, is.numeric, logical(1)))) {
    bulk_tpm   <- as.matrix(tpm_tab)
    tpm_source <- "ESTIMATE pipeline"
    message(sprintf("  Reused TPM from ESTIMATE: %d genes x %d samples",
                    nrow(bulk_tpm), ncol(bulk_tpm)))
  }
}

# -- Option 2: recompute from raw counts + EPIC gene lengths ------------------
if (is.null(bulk_tpm)) {
  message("  Recomputing TPM from raw counts and EPIC gene lengths...")
  if (!file.exists(LENGTHS_PATH)) stop("Missing gene length cache: ", LENGTHS_PATH)
  if (!file.exists(COUNTS_PATH))  stop("Missing bulk counts: ", COUNTS_PATH)
  
  gl <- read.csv(LENGTHS_PATH, stringsAsFactors = FALSE)
  names(gl)[1:2] <- c("gene", "length")
  gl <- gl[gl$length > 0, ]
  gene_lengths <- setNames(gl$length, gl$gene)
  
  raw <- read_excel(COUNTS_PATH, sheet = 1)
  counts <- as.matrix(raw[, -1])
  rownames(counts) <- raw[[1]]
  mode(counts) <- "numeric"
  counts[is.na(counts)] <- 0
  
  bulk_tpm   <- compute_tpm(counts, gene_lengths)
  tpm_source <- "recomputed (EPIC gene lengths)"
  
  col_sums <- colSums(bulk_tpm)
  message(sprintf("  TPM per-sample sums: min %.0f, max %.0f (expected ~1e6)",
                  min(col_sums), max(col_sums)))
}

message(sprintf("  TPM source: %s | %d genes x %d samples",
                tpm_source, nrow(bulk_tpm), ncol(bulk_tpm)))

# ==============================================================================
# 3. LOAD REFERENCE PANEL
# ==============================================================================
message("\n[3/5] Loading BlueprintEncode reference...")
data("BlueprintEncode.xCell2Ref", package = "xCell2")
ref_bp <- BlueprintEncode.xCell2Ref

ct_labels <- rownames(ref_bp@spill_mat)
message(sprintf("  Cell types in reference: %d", length(ct_labels)))

overlap <- length(intersect(rownames(bulk_tpm), ref_bp@genes_used))
message(sprintf("  Bulk/reference gene overlap: %d / %d (%.1f%%)",
                overlap, length(ref_bp@genes_used),
                100 * overlap / length(ref_bp@genes_used)))
if (overlap < MIN_OVERLAP_GENES) {
  warning("Low gene overlap (< ", MIN_OVERLAP_GENES,
          "). Verify HGNC gene symbols.")
}

# ==============================================================================
# 4. RUN xCell 2.0
# ==============================================================================
message("\n[4/5] Running xCell2Analysis (this can take several minutes)...")
xcell2_scores <- tryCatch(
  xCell2Analysis(
    mix            = bulk_tpm,
    xcell2object   = ref_bp,
    spilloverAlpha = SPILLOVER_ALPHA,
    minSharedGenes = MIN_SHARED_FRAC
  ),
  error = function(e) {
    stop("xCell2Analysis failed: ", conditionMessage(e),
         "\nCheck: HGNC gene symbols, >= 5 samples, numeric matrix.",
         call. = FALSE)
  }
)
message(sprintf("  Done: %d cell types x %d samples",
                nrow(xcell2_scores), ncol(xcell2_scores)))

# ==============================================================================
# 5. CLEAN NAMES AND SAVE
#   xCell2 returns cell types in rows; transpose so samples are in rows,
#   consistent with the rest of the pipeline.
# ==============================================================================
message("\n[5/5] Saving scores...")

clean_names <- function(x) {
  x <- gsub(",\\s*", " ", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("-", "_", x)
  gsub("[^A-Za-z0-9_]", "", x)
}

type_names <- clean_names(rownames(xcell2_scores))
if (anyDuplicated(type_names)) {
  type_names <- make.unique(type_names, sep = "_v")
}

scores_df <- as.data.frame(t(xcell2_scores))
colnames(scores_df) <- type_names
scores_df <- cbind(sample = rownames(scores_df), scores_df)

out_csv <- file.path(OUT_DIR, "xcell2_scores.csv")
write.csv(scores_df, out_csv, row.names = FALSE)
message(sprintf("  Written: %s (%d samples x %d cell types)",
                out_csv, nrow(scores_df), ncol(scores_df) - 1))

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 04 (xCell 2.0, ROL2 supplementary) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Reference panel : BlueprintEncode (%d cell types)",
                length(ct_labels)))
message(sprintf("  TPM source      : %s", tpm_source))
message(sprintf("  Gene overlap    : %d genes", overlap))
message(sprintf("  Output          : %s", out_csv))
message(strrep("=", 70))