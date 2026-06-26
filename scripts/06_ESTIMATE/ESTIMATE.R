# ==============================================================================
# Tumour Purity Estimation with ESTIMATE
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 06_ESTIMATE
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Estimate stromal/immune infiltration and tumour purity from bulk expression
#   using ESTIMATE. ESTIMATE is NOT part of the ROL1 consensus: it serves as an
#   orthogonal purity reference for validation (Bland-Altman, script 09) and as
#   the anchor for the optional purity-corrected consensus (script 08).
#   It produces StromalScore, ImmuneScore and ESTIMATEScore per sample, and a
#   derived TumorPurity.
#
# IMPORTANT CAVEAT (purity calibration)
#   The TumorPurity transform
#     TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
#   was calibrated by Yoshihara et al. (2013) on AFFYMETRIX microarray data.
#   Applied to RNA-seq it is an APPROXIMATE estimate and is interpreted as a
#   relative purity proxy. Validation of the deconvolution against ESTIMATE
#   therefore uses the UNCORRECTED consensus (see script 08).
#
# INPUT NORMALISATION
#   Bulk is TPM-normalised (raw counts + EPIC gene-length cache). The TPM
#   matrix is also written to disk and reused by other ROL2 scripts (e.g. 04).
#   If the gene-length cache is unavailable, the script falls back to CPM.
#
# INPUTS
#   data/bulk/oligo_counts_protein.coding.xlsx     (raw counts)
#   scripts/03_EPIC/results/gene_lengths_cache.csv (gene lengths; optional)
#
# OUTPUTS
#   scripts/06_ESTIMATE/results/estimate_scores.csv          (used in 08/09)
#   scripts/06_ESTIMATE/results/bulk_tpm_for_estimate.txt    (reused by 04)
#   scripts/06_ESTIMATE/results/bulk_estimate_filtered.gct
#   scripts/06_ESTIMATE/results/estimate_scores.gct
#   scripts/06_ESTIMATE/figures/*.png
#
# REFERENCE
#   Yoshihara et al. (2013) ESTIMATE. Nat Commun 4:2612.
# ==============================================================================
set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BULK_PATH    <- "data/bulk/oligo_counts_protein.coding.xlsx"
LENGTHS_PATH <- "scripts/03_EPIC/results/gene_lengths_cache.csv"
OUT_DIR      <- "scripts/06_ESTIMATE/results"
FIG_DIR      <- "scripts/06_ESTIMATE/figures"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(estimate)
  library(readxl)
  library(ggplot2)
  library(reshape2)
})

# ESTIMATE's filterCommonGenes() merges the input matrix against the internal
# 'common_genes' signature table (by.x = "GeneSymbol"). That object is a lazy
# package dataset and is NOT attached automatically: filterCommonGenes() then
# fails with "'by' must specify a uniquely valid column". Loading it explicitly
# here is required for the script to run in a clean session or on the cluster.
data(common_genes, package = "estimate")

# Shared figure styling: palette, theme and scale helpers.
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# HELPER: TPM from raw counts and gene lengths
# ==============================================================================
compute_tpm <- function(counts, gene_lengths_bp) {
  shared   <- intersect(rownames(counts), names(gene_lengths_bp))
  counts_f <- counts[shared, , drop = FALSE]
  rpk      <- counts_f / gene_lengths_bp[shared]
  t(t(rpk) / colSums(rpk) * 1e6)
}

# ==============================================================================
# 1. LOAD BULK AND NORMALISE (TPM, CPM fallback)
# ==============================================================================
message("\n[1/6] Loading bulk and normalising...")
bulk_df <- as.data.frame(read_excel(BULK_PATH))
rownames(bulk_df) <- bulk_df[[1]]
bulk_df <- bulk_df[, -1, drop = FALSE]
bulk_mat <- as.matrix(bulk_df)
storage.mode(bulk_mat) <- "numeric"
bulk_mat[is.na(bulk_mat)] <- 0
message(sprintf("  Bulk: %d genes x %d samples", nrow(bulk_mat), ncol(bulk_mat)))

if (file.exists(LENGTHS_PATH)) {
  gl <- read.csv(LENGTHS_PATH, stringsAsFactors = FALSE)
  names(gl)[1:2] <- c("gene", "length")
  gl <- gl[gl$length > 0, ]
  gene_lengths <- setNames(gl$length, gl$gene)
  bulk_norm   <- compute_tpm(bulk_mat, gene_lengths)
  norm_method <- "TPM"
  message(sprintf("  TPM normalisation | per-sample sums ~%.0f",
                  mean(colSums(bulk_norm))))
} else {
  bulk_norm   <- t(t(bulk_mat) / colSums(bulk_mat) * 1e6)
  norm_method <- "CPM (fallback)"
  message("  WARNING: gene length cache not found; using CPM fallback.")
}

# ---- 2. Identify sample / Glioma / Oligo columns ----------------------------
GLIOMA_COL <- "Glioma"
OLIGO_COL  <- "Oligo"

stopifnot(GLIOMA_COL %in% colnames(cons), OLIGO_COL %in% colnames(cons))

# Sample barcodes are the ROW NAMES (first CSV column was read as row.names).
cons$sample_id <- substr(gsub("\\.", "-", rownames(cons)), 1, 12)

# ==============================================================================
# 3. FILTER GENES SHARED WITH ESTIMATE SIGNATURES
# ==============================================================================
message("\n[3/6] Filtering genes against ESTIMATE signatures...")
filterCommonGenes(
  input.f  = tpm_path,
  output.f = filtered_path,
  id       = "GeneSymbol"
)

# ==============================================================================
# 4. COMPUTE ESTIMATE SCORES
#   platform = "illumina": appropriate for the score computation on RNA-seq.
#   (The purity transform below carries the Affymetrix-calibration caveat.)
# ==============================================================================
message("\n[4/6] Computing ESTIMATE scores...")
estimateScore(
  input.ds  = filtered_path,
  output.ds = scores_path,
  platform  = "illumina"
)

# ==============================================================================
# 5. PARSE SCORES AND DERIVE TUMOUR PURITY
#   The .gct has 2 header lines; row 1 = gene/score names, col 1 = NAME,
#   col 2 = Description. Parse defensively and select score rows by name so the
#   code does not depend on their exact order or on the platform variant.
# ==============================================================================
message("\n[5/6] Parsing scores and deriving tumour purity...")
scores_raw <- read.table(scores_path, skip = 2, header = TRUE,
                         row.names = 1, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
scores_raw <- scores_raw[, -1, drop = FALSE]   # drop the Description column
scores_df  <- as.data.frame(t(scores_raw))

# Keep only the canonical ESTIMATE score columns, by name (robust to ordering).
expected <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
                      colnames(scores_df))
if (length(expected) < 3L) {
  stop("Expected StromalScore/ImmuneScore/ESTIMATEScore in the .gct; found: ",
       paste(colnames(scores_df), collapse = ", "))
}
scores_df <- scores_df[, expected, drop = FALSE]
scores_df[] <- lapply(scores_df, function(x) as.numeric(as.character(x)))

# Yoshihara et al. (2013) purity transform (Affymetrix-calibrated; see header).
scores_df$TumorPurity <- cos(0.6049872018 +
                               0.0001467884 * scores_df$ESTIMATEScore)
scores_df$TumorPurity <- pmax(0, pmin(1, scores_df$TumorPurity))

message("  Tumour purity summary:")
print(summary(scores_df$TumorPurity))

# ESTIMATE writes sample IDs with dots; restore the TCGA hyphen convention so
# downstream merges (scripts 08/09) match correctly.
rownames(scores_df) <- gsub("\\.", "-", rownames(scores_df))

# ==============================================================================
# 6. SAVE RESULTS
# ==============================================================================
out_csv <- file.path(OUT_DIR, "estimate_scores.csv")
write.csv(scores_df, out_csv)
message(sprintf("\n[6/6] Scores written: %s", out_csv))
message(sprintf("  Columns: %s", paste(colnames(scores_df), collapse = ", ")))

# ==============================================================================
# FIGURES
# ==============================================================================
message("\nGenerating figures...")
glioma_col <- tfm_colors[["Glioma"]]
scores_df$sample <- rownames(scores_df)
sample_order <- scores_df[order(scores_df$TumorPurity,
                                decreasing = TRUE), "sample"]
scores_df$sample <- factor(scores_df$sample, levels = sample_order)
purity_median <- median(scores_df$TumorPurity)

# -- Figure 1: tumour purity distribution -------------------------------------
p_dist <- ggplot(scores_df, aes(x = TumorPurity)) +
  geom_histogram(bins = 30, fill = glioma_col, color = "white", alpha = 0.85) +
  geom_vline(xintercept = purity_median, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  annotate("text", x = purity_median + 0.01, y = Inf, vjust = 1.5,
           label = sprintf("Median = %.3f", purity_median), size = 3.5) +
  labs(
    title    = "Tumour purity distribution (ESTIMATE)",
    subtitle = "Yoshihara et al. 2013 | 122 TCGA oligodendroglioma samples",
    x        = "Tumour purity (approximate; see methods)",
    y        = "Number of samples"
  ) +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(FIG_DIR, "estimate_tumor_purity_distribution.png"),
       p_dist, width = 8, height = 5, dpi = 300)

# -- Figure 2: tumour purity per sample ---------------------------------------
p_bar <- ggplot(scores_df, aes(x = sample, y = TumorPurity)) +
  geom_col(fill = glioma_col, alpha = 0.85, width = 1) +
  geom_hline(yintercept = purity_median, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title    = "Tumour purity per sample (ESTIMATE)",
    subtitle = "Ordered by decreasing purity | dashed line = median",
    x        = "Sample", y = "Tumour purity"
  ) +
  theme_tfm()
ggsave(file.path(FIG_DIR, "estimate_tumor_purity_barplot.png"),
       p_bar, width = 14, height = 5, dpi = 300)

# -- Figure 3: stromal vs immune score ----------------------------------------
p_scatter <- ggplot(scores_df, aes(x = StromalScore, y = ImmuneScore,
                                   color = TumorPurity)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
                        midpoint = purity_median, name = "Tumour\npurity") +
  labs(
    title    = "Stromal vs immune score (ESTIMATE)",
    subtitle = "Lower purity corresponds to higher stromal/immune infiltration",
    x        = "Stromal score", y = "Immune score"
  ) +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(FIG_DIR, "estimate_stromal_vs_immune.png"),
       p_scatter, width = 8, height = 6, dpi = 300)

message("Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 06 (ESTIMATE) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Role          : orthogonal purity reference (not in consensus)"))
message(sprintf("  Normalisation : %s", norm_method))
message(sprintf("  Samples       : %d", nrow(scores_df)))
message(sprintf("  Tumour purity : median %.3f, range [%.3f, %.3f]",
                purity_median, min(scores_df$TumorPurity),
                max(scores_df$TumorPurity)))
message(sprintf("  Main output   : %s", out_csv))
message("  NOTE: purity transform is Affymetrix-calibrated (approximate on RNA-seq).")
message(strrep("=", 70))