# ==============================================================================
# Immune Deconvolution with quanTIseq (ROL2, supplementary)
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 07_quanTIseq
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Estimate immune cell fractions from bulk RNA-seq with quanTIseq, via the
#   immunedeconv wrapper. quanTIseq is a ROL2 (reference-free) method: it is NOT
#   part of the ROL1 consensus and serves only as orthogonal, supplementary
#   immune-context evidence for the cold-tumour narrative.
#
# INPUT NORMALISATION
#   quanTIseq requires TPM. TPM is computed from raw CellRanger-style counts and
#   the EPIC gene-length cache (counts / gene_length_kb, then per-sample scaling
#   to 1e6). This matches the compute_tpm() used in scripts 05/06; a naive
#   counts/1000 transform is NOT TPM (it is CPM in disguise and biases long vs
#   short genes), so it is deliberately avoided here.
#
# INPUTS
#   data/bulk/oligo_counts_protein.coding.xlsx     (raw counts)
#   scripts/03_EPIC/results/gene_lengths_cache.csv (gene lengths)
#
# OUTPUTS
#   scripts/07_quanTIseq/results/quantiseq_proportions.csv
#   scripts/07_quanTIseq/results/quantiseq_session_info.txt  (for Table S1)
#   scripts/07_quanTIseq/figures/*.png
#
# REFERENCE
#   Finotello et al. (2019) quanTIseq. Genome Med 11:34.
# ==============================================================================
set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
COUNTS_PATH  <- "data/bulk/oligo_counts_protein.coding.xlsx"
LENGTHS_PATH <- "scripts/03_EPIC/results/gene_lengths_cache.csv"
OUT_DIR      <- "scripts/07_quanTIseq/results"
FIG_DIR      <- "scripts/07_quanTIseq/figures"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(immunedeconv)
  library(readxl)
  library(ggplot2)
  library(reshape2)
})

# Shared figure styling (theme only; quanTIseq cell types are not in the
# cell-type palette, so a qualitative fill scale is used).
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# HELPER: TPM from raw counts and gene lengths
#   TPM = (counts / length_kb) / sum(counts / length_kb) * 1e6
#   gene_lengths are in base pairs; division by 1000 converts to kb.
# ==============================================================================
compute_tpm <- function(counts, gene_lengths_bp) {
  shared   <- intersect(rownames(counts), names(gene_lengths_bp))
  counts_f <- counts[shared, , drop = FALSE]
  rpk      <- counts_f / (gene_lengths_bp[shared] / 1000)
  t(t(rpk) / colSums(rpk) * 1e6)
}

# ==============================================================================
# 1. LOAD BULK AND NORMALISE TO TPM
# ==============================================================================
message("\n[1/5] Loading bulk and normalising to TPM...")
if (!file.exists(LENGTHS_PATH)) stop("Missing gene length cache: ", LENGTHS_PATH)

raw <- read_excel(COUNTS_PATH, sheet = 1)
counts <- as.matrix(raw[, -1])
rownames(counts) <- raw[[1]]
mode(counts) <- "numeric"
counts[is.na(counts)] <- 0
message(sprintf("  Bulk: %d genes x %d samples", nrow(counts), ncol(counts)))

gl <- read.csv(LENGTHS_PATH, stringsAsFactors = FALSE)
names(gl)[1:2] <- c("gene", "length")
gl <- gl[gl$length > 0, ]
gene_lengths <- setNames(gl$length, gl$gene)

bulk_tpm <- compute_tpm(counts, gene_lengths)
message(sprintf("  TPM: %d genes x %d samples | per-sample sums ~%.0f",
                nrow(bulk_tpm), ncol(bulk_tpm), mean(colSums(bulk_tpm))))

# ==============================================================================
# 2. DECONVOLUTION WITH quanTIseq
#   tumor = TRUE       : applies the tumour-aware signature (uncharacterised
#                        cells are reported as 'uncharacterized cell').
#   scale_mrna = TRUE  : corrects for cell-type mRNA content differences.
# ==============================================================================
message("\n[2/5] Running quanTIseq...")
result_quantiseq <- deconvolute(
  gene_expression = bulk_tpm,
  method          = "quantiseq",
  tumor           = TRUE,
  scale_mrna      = TRUE
)
message(sprintf("  Estimated cell types: %s",
                paste(result_quantiseq$cell_type, collapse = ", ")))

# ==============================================================================
# 3. RESHAPE TO SAMPLES x CELL TYPES (pipeline convention)
# ==============================================================================
message("\n[3/5] Reshaping results...")
quant_df           <- as.data.frame(t(result_quantiseq[, -1]))
colnames(quant_df) <- result_quantiseq$cell_type
quant_df           <- cbind(sample = rownames(quant_df), quant_df)

# ==============================================================================
# 4. SAVE RESULTS AND VERSION INFO (Table S1)
# ==============================================================================
message("\n[4/5] Saving results...")
out_csv <- file.path(OUT_DIR, "quantiseq_proportions.csv")
write.csv(quant_df, out_csv, row.names = FALSE)
message(sprintf("  Written: %s (%d samples x %d cell types)",
                out_csv, nrow(quant_df), ncol(quant_df) - 1))

# Record package versions for the supplementary version table (Table S1).
ver_path <- file.path(OUT_DIR, "quantiseq_session_info.txt")
writeLines(
  c(sprintf("immunedeconv: %s",
            as.character(utils::packageVersion("immunedeconv"))),
    sprintf("Run date    : %s", Sys.Date()),
    "----------------------------------------------------------------------",
    capture.output(sessionInfo())),
  ver_path
)
message(sprintf("  Version info written: %s", ver_path))
message(sprintf("  immunedeconv version: %s",
                as.character(utils::packageVersion("immunedeconv"))))

# ==============================================================================
# 5. FIGURES
# ==============================================================================
message("\n[5/5] Generating figures...")
quant_long <- reshape2::melt(quant_df, id.vars = "sample",
                             variable.name = "cell_type",
                             value.name    = "proportion")

# -- Figure 1: stacked composition per sample ---------------------------------
p_bar <- ggplot(quant_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(
    title    = "Immune composition (quanTIseq)",
    subtitle = "Finotello et al. 2019 | 122 TCGA oligodendroglioma samples",
    x        = "Sample", y = "Proportion", fill = "Cell type"
  ) +
  theme_tfm()
ggsave(file.path(FIG_DIR, "quantiseq_composition.png"),
       p_bar, width = 14, height = 6, dpi = 300)

# -- Figure 2: per cell-type distribution -------------------------------------
p_box <- ggplot(quant_long, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  geom_jitter(width = 0.2, alpha = 0.35, size = 0.8) +
  labs(
    title    = "Proportion distribution per cell type (quanTIseq)",
    subtitle = "ROL2 supplementary immune-context evidence",
    x        = "Cell type", y = "Proportion"
  ) +
  theme_tfm(show_x_text = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave(file.path(FIG_DIR, "quantiseq_boxplot.png"),
       p_box, width = 10, height = 6, dpi = 300)

message("  Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 07 (quanTIseq, ROL2 supplementary) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Method      : quanTIseq via immunedeconv"))
message(sprintf("  Cell types  : %d", ncol(quant_df) - 1))
message(sprintf("  Samples     : %d", nrow(quant_df)))
message(sprintf("  Output      : %s", out_csv))
message(strrep("=", 70))