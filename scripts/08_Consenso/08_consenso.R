# ==============================================================================
# ROL1 Consensus: MuSiC + BayesPrism
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 08_Consensus
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Aggregate the two reference-based deconvolution methods (MuSiC, BayesPrism)
#   into a single consensus estimate at 15-cell-type resolution, and report
#   their per-cell-type concordance as evidence of robustness.
#
# RATIONALE (why these two methods, and why their average)
#   MuSiC and BayesPrism are the two top-performing reference-based methods in
#   independent benchmarks (Hu & Chikina, 2024) and the only two run at full
#   15-type resolution in this pipeline. EPIC is NOT included: at 15 types its
#   reference matrix is near-singular (collinear myeloid subtypes), so it is
#   run separately at 4-category resolution and contributes only a reference-
#   coverage metric (otherCells; script 03). With two methods the consensus is
#   their per-sample, per-cell-type mean; robustness is therefore assessed
#   directly through the Spearman concordance between the two methods rather
#   than via outlier-robust aggregation.
#
# PURITY CORRECTION AND CIRCULARITY
#   An optional purity-anchored consensus is produced by setting the malignant
#   (Glioma) fraction to the ESTIMATE-derived tumour purity and rescaling the
#   remaining types to fill (1 - purity). Because this overwrites Glioma with
#   ESTIMATE, any validation against ESTIMATE (Bland-Altman, script 09) MUST
#   use the UNCORRECTED consensus (consensus_proportions.csv) to avoid
#   circularity. The corrected file is for downstream uses that need an
#   absolute malignant fraction only.
#
# INPUTS
#   scripts/01_MuSiC/results/music_proportions.csv
#   scripts/02_BayesPrism/results/bayesprism_proportions.csv
#   scripts/06_ESTIMATE/results/estimate_scores.csv   (optional)
#
# OUTPUTS
#   scripts/08_Consenso/results/consensus_proportions.csv                 (main)
#   scripts/08_Consenso/results/consensus_proportions_purity_corrected.csv
#   scripts/08_Consenso/results/consensus_n_methods_per_cell.csv
#   scripts/08_Consenso/results/spearman_concordance_methods.csv
#   scripts/08_Consenso/figures/*.png
#
# REFERENCES
#   Hu & Chikina (2024)      Deconvolution benchmark. Genome Biol 25:169.
#   Wang et al. (2019)       MuSiC. Nat Commun 10:380.
#   Chu et al. (2022)        BayesPrism. Nat Cancer 3:505-517.
#   Yoshihara et al. (2013)  ESTIMATE. Nat Commun 4:2612.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PATH_MUSIC <- "scripts/01_MuSiC/results/music_proportions.csv"
PATH_BAYES <- "scripts/02_BayesPrism/results/bayesprism_proportions.csv"
PATH_ESTIM <- "scripts/06_ESTIMATE/results/estimate_scores.csv"

TAXONOMY_LEVEL <- "SubAssignment"  # cell-type resolution used in the pipeline

OUT_DIR <- "scripts/08_Consenso/results"
FIG_DIR <- "scripts/08_Consenso/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(dplyr)
})

# Shared figure styling: palettes, theme and scale helpers.
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# 1. LOAD DECONVOLUTION RESULTS
# ==============================================================================
message("\n[1/5] Loading MuSiC and BayesPrism results...")

read_props <- function(path, label) {
  m <- read.csv(path, row.names = 1, check.names = FALSE)
  # Ensure samples are in rows (transpose if cell types outnumber samples).
  if (ncol(m) > nrow(m)) {
    m <- as.data.frame(t(m))
    message(sprintf("  %s transposed: samples now in rows", label))
  }
  m
}

music <- read_props(PATH_MUSIC, "MuSiC")
bayes <- read_props(PATH_BAYES, "BayesPrism")

message(sprintf("  MuSiC:      %d samples x %d cell types",
                nrow(music), ncol(music)))
message(sprintf("  BayesPrism: %d samples x %d cell types",
                nrow(bayes), ncol(bayes)))

# ==============================================================================
# 2. ALIGN SAMPLES AND CELL TYPES
# ==============================================================================
message("\n[2/5] Aligning samples and cell types...")
common_samples <- intersect(rownames(music), rownames(bayes))
common_types   <- intersect(colnames(music), colnames(bayes))

if (length(common_samples) == 0) {
  stop("No common samples across methods. Check row-name formats.")
}
if (length(common_types) == 0) {
  stop("No common cell types across methods. Check column-name formats.")
}

message(sprintf("  Common samples:    %d", length(common_samples)))
message(sprintf("  Common cell types: %d", length(common_types)))

# Report any method-exclusive cell types (excluded from the consensus).
only_music <- setdiff(colnames(music), common_types)
only_bayes <- setdiff(colnames(bayes), common_types)
if (length(only_music) > 0) {
  message("  Types only in MuSiC (excluded): ",
          paste(only_music, collapse = ", "))
}
if (length(only_bayes) > 0) {
  message("  Types only in BayesPrism (excluded): ",
          paste(only_bayes, collapse = ", "))
}

music_mat <- as.matrix(music[common_samples, common_types])
bayes_mat <- as.matrix(bayes[common_samples, common_types])

# ==============================================================================
# 3. CONSENSUS = MEAN OF THE TWO METHODS
#   With two methods the consensus is their average; a cell is only assigned a
#   consensus value when both methods provide an estimate (otherwise NA).
# ==============================================================================
message("\n[3/5] Computing consensus (mean of the two methods)...")

both_present  <- !is.na(music_mat) & !is.na(bayes_mat)
consensus_mat <- (music_mat + bayes_mat) / 2
consensus_mat[!both_present] <- NA

n_methods_mat <- (!is.na(music_mat)) + (!is.na(bayes_mat))

n_incomplete <- sum(apply(consensus_mat, 1, function(x) any(is.na(x))))
if (n_incomplete > 0) {
  message(sprintf("  WARNING: %d samples have cell types missing in one method",
                  n_incomplete))
}

# Renormalise each sample to sum to 1.
consensus_norm <- t(apply(consensus_mat, 1, function(x) {
  s <- sum(x, na.rm = TRUE)
  if (!is.na(s) && s > 0) x / s else x
}))

message(sprintf("  Per-sample sums == 1: %s",
                all(abs(rowSums(consensus_norm, na.rm = TRUE) - 1) < 1e-9)))
message(sprintf("  NA cells in consensus: %d", sum(is.na(consensus_norm))))

# ==============================================================================
# 4. OPTIONAL PURITY-ANCHORED CONSENSUS (ESTIMATE)
#   See header: do NOT use this file to validate against ESTIMATE.
# ==============================================================================
message("\n[4/5] Building purity-anchored consensus (ESTIMATE)...")
consensus_corr <- consensus_norm

if (file.exists(PATH_ESTIM)) {
  estim <- read.csv(PATH_ESTIM, row.names = 1, check.names = FALSE)
  # ESTIMATE writes sample IDs with dots; restore the TCGA hyphen convention.
  rownames(estim) <- gsub("\\.", "-", rownames(estim))
  
  col_estimate <- grep("ESTIMATEScore", colnames(estim),
                       value = TRUE, ignore.case = TRUE)[1]
  col_tumor    <- grep("Glioma", colnames(consensus_corr), value = TRUE)[1]
  
  if (!is.na(col_estimate) && !is.na(col_tumor)) {
    common_estim <- intersect(common_samples, rownames(estim))
    message(sprintf("  Samples with ESTIMATE scores: %d", length(common_estim)))
    
    # setNames() is required: estim[rows, col] returns an UNNAMED vector, so
    # without explicit names purity[sample] would return NA.
    purity <- setNames(
      cos(0.6049872018 + 0.0001467884 * estim[common_estim, col_estimate]),
      common_estim
    )
    purity <- pmax(pmin(purity, 1), 0)
    message(sprintf("  Tumour purity: median %.3f, range [%.3f, %.3f]",
                    median(purity), min(purity), max(purity)))
    
    consensus_corr <- consensus_norm[common_estim, , drop = FALSE]
    cols_other     <- setdiff(colnames(consensus_corr), col_tumor)
    
    for (s in common_estim) {
      p         <- purity[s]
      sum_other <- sum(consensus_corr[s, cols_other], na.rm = TRUE)
      consensus_corr[s, col_tumor] <- p
      if (sum_other > 0) {
        consensus_corr[s, cols_other] <-
          consensus_corr[s, cols_other] * (1 - p) / sum_other
      }
    }
    message(sprintf("  Purity anchoring applied to %d samples (tumour column: %s)",
                    length(common_estim), col_tumor))
  } else {
    message("  WARNING: ESTIMATEScore or Glioma column not found; ",
            "purity-anchored consensus = uncorrected consensus.")
  }
} else {
  message("  ESTIMATE file not found; purity-anchored consensus = uncorrected.")
}

# ==============================================================================
# 5. METHOD CONCORDANCE (robustness evidence)
#   Per-cell-type Spearman correlation between MuSiC and BayesPrism.
# ==============================================================================
message("\n[5/5] Computing MuSiC vs BayesPrism concordance...")

concordance <- do.call(rbind, lapply(common_types, function(ct) {
  v1 <- music_mat[, ct]; v2 <- bayes_mat[, ct]
  if (sum(!is.na(v1) & !is.na(v2)) < 5) return(NULL)
  data.frame(
    pair     = "MuSiC vs BayesPrism",
    celltype = ct,
    spearman = round(cor(v1, v2, method = "spearman",
                         use = "complete.obs"), 3)
  )
}))

message("Spearman concordance by cell type:")
print(concordance)
message(sprintf("Median concordance across cell types: %.3f",
                median(concordance$spearman, na.rm = TRUE)))

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
message("\nSaving results...")
write.csv(consensus_norm,
          file.path(OUT_DIR, "consensus_proportions.csv"))
write.csv(consensus_corr,
          file.path(OUT_DIR, "consensus_proportions_purity_corrected.csv"))
write.csv(as.data.frame(n_methods_mat),
          file.path(OUT_DIR, "consensus_n_methods_per_cell.csv"))
write.csv(concordance,
          file.path(OUT_DIR, "spearman_concordance_methods.csv"),
          row.names = FALSE)

# ==============================================================================
# FIGURES
# ==============================================================================
message("\nGenerating figures...")

# -- Figure 1: consensus cell-type composition --------------------------------
#   Stacked barplot: Y axis in %, samples ordered by malignant fraction, and a
#   biological stacking order (malignant -> stromal-normal -> myeloid lineage ->
#   lymphoid -> vascular) so the cold-tumour structure reads top to bottom.
comp_wide        <- as.data.frame(consensus_norm)
comp_wide$sample <- rownames(comp_wide)
comp_long <- reshape2::melt(comp_wide, id.vars = "sample",
                            variable.name = "cell_type",
                            value.name = "proportion")
comp_long$cell_type <- tfm_normalise_types(comp_long$cell_type)

# Order samples by malignant (Glioma) fraction.
if ("Glioma" %in% colnames(comp_wide)) {
  sample_order <- comp_wide[order(comp_wide$Glioma, decreasing = TRUE), "sample"]
  comp_long$sample <- factor(comp_long$sample, levels = sample_order)
}

# Biological stacking order. rev() puts Glioma at the bottom of the stack.
stack_order <- c("Glioma", "Oligo", "Proliferating",
                 "h.microglia", "i.microglia", "a.microglia", "AP.microglia",
                 "s.mac.1", "s.mac.2", "MDSC", "Myeloid",
                 "TCells", "CD8.TCells", "Endo", "Pericytes")
stack_order <- intersect(stack_order, unique(comp_long$cell_type))
comp_long$cell_type <- factor(comp_long$cell_type, levels = rev(stack_order))

p_composition <- ggplot(comp_long,
                        aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.005),
                     labels = function(x) paste0(x * 100, "%")) +
  tfm_scale_fill(drop = TRUE,
                 guide = ggplot2::guide_legend(reverse = TRUE, ncol = 1)) +
  labs(
    title    = "ROL1 consensus cell-type composition",
    subtitle = sprintf("MuSiC + BayesPrism (mean) | 122 TCGA samples, ordered by malignant fraction"),
    x        = "Sample (n = 122)",
    y        = "Proportion", fill = "Cell type"
  ) +
  theme_tfm() +
  theme(axis.line.x = ggplot2::element_blank())
ggsave(file.path(FIG_DIR, "consensus_composition.png"),
       p_composition, width = 16, height = 5, dpi = 300)

# -- Figure 2: MuSiC vs BayesPrism concordance per cell type ------------------
concordance$cell_type <- tfm_normalise_types(concordance$celltype)

p_concordance <- ggplot(concordance,
                        aes(x = reorder(cell_type, spearman), y = spearman,
                            fill = cell_type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "black", linewidth = 0.4) +
  tfm_scale_fill() +
  scale_y_continuous(limits = c(min(0, min(concordance$spearman)), 1),
                     expand = c(0.02, 0)) +
  labs(
    title    = "MuSiC vs BayesPrism concordance",
    subtitle = "Per-cell-type Spearman correlation | dashed line: r = 0.70",
    x        = "Cell type", y = "Spearman r"
  ) +
  theme_tfm(show_x_text = TRUE) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "spearman_concordance.png"),
       p_concordance, width = 12, height = 5, dpi = 300)

message("Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 08 (CONSENSUS) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Methods        : MuSiC + BayesPrism (mean)"))
message(sprintf("  Common samples : %d", length(common_samples)))
message(sprintf("  Cell types     : %d", length(common_types)))
message(sprintf("  Median MuSiC-BayesPrism concordance: %.3f",
                median(concordance$spearman, na.rm = TRUE)))
message(sprintf("  Main output    : %s",
                file.path(OUT_DIR, "consensus_proportions.csv")))
message("  NOTE: validate against ESTIMATE using the UNCORRECTED consensus.")
message(strrep("=", 70))