# ==============================================================================
# Immune/Stromal Scoring with MCP-counter (ROL2, supplementary)
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 05_MCP_Counter
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Compute abundance scores for immune and stromal populations in the bulk
#   cohort using MCP-counter. This is a reference-free (ROL2) method used only
#   for SUPPLEMENTARY, directional validation; its scores are arbitrary
#   abundance estimates (NOT proportions) and are not combined with the ROL1
#   consensus. MCP-counter uses its own internal marker panel, so it does not
#   depend on the LGG-04 scRNA-seq reference.
#
# INPUT NORMALISATION
#   Expression is TPM-normalised using raw counts and the EPIC gene-length
#   cache (corrects per-gene length, unlike a constant-factor CPM).
#
# INPUTS
#   data/bulk/oligo_counts_protein.coding.xlsx     (raw counts)
#   scripts/03_EPIC/results/gene_lengths_cache.csv (gene lengths)
#
# OUTPUTS
#   scripts/05_MCP_Counter/results/mcpcounter_scores.csv
#   scripts/05_MCP_Counter/figures/*.png
#
# NOTE
#   MCP-counter is installed from GitHub on first run (requires internet).
#
# REFERENCE
#   Becht et al. (2016) MCP-counter. Genome Biol 17:218.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
COUNTS_PATH  <- "data/bulk/oligo_counts_protein.coding.xlsx"
LENGTHS_PATH <- "scripts/03_EPIC/results/gene_lengths_cache.csv"
OUT_DIR      <- "scripts/05_MCP_Counter/results"
FIG_DIR      <- "scripts/05_MCP_Counter/figures"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(reshape2)
})

# Shared figure styling (theme only; MCP-counter cell types are not in the
# cell-type palette, so a qualitative scale is used for fills).
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# 1. INSTALL MCP-counter (first run only; requires internet)
# ==============================================================================
message("\n[1/5] Checking MCP-counter installation...")
if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  message("  MCPcounter not found. Installing from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
  tryCatch(
    devtools::install_github("ebecht/MCPcounter", ref = "master",
                             subdir = "Source", upgrade = "never"),
    error = function(e) {
      stop("MCPcounter installation failed: ", conditionMessage(e),
           call. = FALSE)
    }
  )
}
suppressPackageStartupMessages(library(MCPcounter))
message("  MCP-counter ready.")

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
# 2. LOAD BULK AND NORMALISE TO TPM
# ==============================================================================
message("\n[2/5] Loading bulk and normalising to TPM...")
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
col_sums <- colSums(bulk_tpm)
message(sprintf("  TPM: %d genes x %d samples | per-sample sums ~%.0f",
                nrow(bulk_tpm), ncol(bulk_tpm), mean(col_sums)))

# ==============================================================================
# 3. RUN MCP-counter
#   Returns abundance scores (cell types x samples), not proportions.
# ==============================================================================
message("\n[3/5] Running MCP-counter...")
mcp_scores <- MCPcounter.estimate(
  expression   = bulk_tpm,
  featuresType = "HUGO_symbols"
)
message(sprintf("  Estimated %d populations x %d samples",
                nrow(mcp_scores), ncol(mcp_scores)))

# ==============================================================================
# 4. SAVE SCORES (samples in rows, consistent with the pipeline)
# ==============================================================================
message("\n[4/5] Saving scores...")
scores_df <- as.data.frame(t(mcp_scores))
scores_df <- cbind(sample = rownames(scores_df), scores_df)

out_csv <- file.path(OUT_DIR, "mcpcounter_scores.csv")
write.csv(scores_df, out_csv, row.names = FALSE)
message(sprintf("  Written: %s (%d samples x %d populations)",
                out_csv, nrow(scores_df), ncol(scores_df) - 1))

# ==============================================================================
# 5. FIGURES
# ==============================================================================
message("\n[5/5] Generating figures...")

scores_long <- reshape2::melt(scores_df, id.vars = "sample",
                              variable.name = "cell_type",
                              value.name = "score")

# -- Figure 1: per-sample abundance scores ------------------------------------
p_bar <- ggplot(scores_long, aes(x = sample, y = score, fill = cell_type)) +
  geom_col(width = 0.8) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title    = "MCP-counter abundance scores (ROL2, supplementary)",
    subtitle = "122 TCGA oligodendroglioma samples | TPM input",
    x        = "Sample", y = "Abundance score", fill = "Population"
  ) +
  theme_tfm()

ggsave(file.path(FIG_DIR, "mcpcounter_composition.png"),
       p_bar, width = 14, height = 6, dpi = 300)

# -- Figure 2: score distribution per population ------------------------------
p_box <- ggplot(scores_long, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title    = "MCP-counter score distribution (ROL2, supplementary)",
    subtitle = "122 TCGA oligodendroglioma samples",
    x        = "Population", y = "Abundance score"
  ) +
  theme_tfm(show_x_text = TRUE) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "mcpcounter_boxplot.png"),
       p_box, width = 10, height = 6, dpi = 300)

message("  Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 05 (MCP-counter, ROL2 supplementary) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Populations : %d", nrow(mcp_scores)))
message(sprintf("  Samples     : %d", ncol(mcp_scores)))
message(sprintf("  Output      : %s", out_csv))
message(strrep("=", 70))