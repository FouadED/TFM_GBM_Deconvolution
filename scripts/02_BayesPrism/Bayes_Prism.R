# ==============================================================================
# Cell-Type Deconvolution with BayesPrism
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 02_BayesPrism
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Estimate cell-type proportions in 122 TCGA oligodendroglioma bulk RNA-seq
#   samples and recover the purified malignant (Glioma) transcriptome for
#   downstream analyses (stemness index and marker identification).
#
# METHOD
#   BayesPrism (Chu et al., 2022) jointly infers, via Gibbs sampling:
#     (1) cell-type proportions (theta), refined in a final optimisation step;
#     (2) sample-specific purified expression per cell type (the Z tensor).
#   It is robust to cross-platform/cross-cohort batch effects by Bayesian
#   design, which suits a scRNA-seq reference and a bulk cohort from different
#   sources. key = "Glioma" designates the malignant compartment so its
#   purified expression can be extracted for the stemness analysis.
#
# INPUTS
#   data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds
#       Annotated scRNA-seq reference with RAW CellRanger counts (GSE182109).
#   data/bulk/oligo_counts_protein.coding.xlsx
#       Bulk RNA-seq raw counts, 122 TCGA samples (protein-coding genes).
#
# OUTPUTS
#   scripts/02_BayesPrism/results/bayesprism_proportions.csv
#   scripts/02_BayesPrism/results/bayesprism_full_result.rds
#   scripts/02_BayesPrism/results/bayesprism_expression_glioma.rds
#   scripts/02_BayesPrism/results/bayesprism_expression_glioma.csv
#   scripts/02_BayesPrism/figures/bayesprism_composition.png
#   scripts/02_BayesPrism/figures/bayesprism_boxplot.png
#
# REPRODUCIBILITY
#   Run from the project root. Deterministic given the fixed seed.
#   Computationally heavy: intended to run on the HPC cluster (>= 16 cores).
#
# REFERENCES
#   Chu et al. (2022)         BayesPrism. Nat Cancer 3:505-517.
#   Abdelfattah et al. (2022) scRNA-seq glioma atlas. Nat Commun 13:767.
#   Malta et al. (2018)       Stemness indices. Cell 173:338-354.
#   Avila Cobos et al. (2020) Deconvolution benchmarking. Nat Commun 11:5650.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SEURAT_PATH    <- "data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds"
BULK_PATH      <- "data/bulk/oligo_counts_protein.coding.xlsx"
OUT_DIR        <- "scripts/02_BayesPrism/results"
FIG_DIR        <- "scripts/02_BayesPrism/figures"

MIN_CELLS       <- 20    # minimum cells per type to build a stable profile
MAX_CELLS_TYPE  <- 500   # subsampling cap per type (memory management)
MIN_COMMON_GENE <- 1000  # abort if fewer shared genes (nomenclature mismatch)
N_CORES         <- max(1, parallel::detectCores() - 1)  # set to 16 on the cluster

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(BayesPrism)
  library(Seurat)
  library(readxl)
  library(ggplot2)
  library(reshape2)
})

# Shared figure styling: cell-type palette, theme and scale helpers.
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# 1. LOAD AND PREPARE THE scRNA-seq REFERENCE
# ==============================================================================
message("\n[1] Loading scRNA-seq reference...")
seurat_obj <- readRDS(SEURAT_PATH)

# Resolve the working cell-type label: prefer fine-grained SubAssignment,
# fall back to coarse Assignment when the subtype is not annotated.
assignment     <- as.character(seurat_obj@meta.data$Assignment)
sub_assignment <- as.character(seurat_obj@meta.data$SubAssignment)
seurat_obj$cell_label <- ifelse(
  is.na(sub_assignment) | sub_assignment %in% c("", "NA"),
  assignment, sub_assignment
)

# Retain only cell types with enough cells for a reliable profile
# (Avila Cobos et al., 2020).
type_counts <- table(seurat_obj$cell_label)
valid_types <- names(type_counts[type_counts >= MIN_CELLS])

message(sprintf("Valid cell types (>= %d cells): %d",
                MIN_CELLS, length(valid_types)))
message("Excluded (too few cells):")
print(type_counts[type_counts < MIN_CELLS])

seurat_ref <- seurat_obj[, seurat_obj$cell_label %in% valid_types]

# Subsample to MAX_CELLS_TYPE per type to keep the Gibbs sampler tractable.
set.seed(123)
cells_keep <- unlist(lapply(valid_types, function(ct) {
  idx <- which(seurat_ref$cell_label == ct)
  sample(idx, min(length(idx), MAX_CELLS_TYPE))
}))
seurat_ref <- seurat_ref[, cells_keep]
message(sprintf("Reference cells after subsampling: %d", ncol(seurat_ref)))

# Sanity check: BayesPrism expects raw counts. Confirm the reference is
# integer-valued, guarding against an accidentally log-normalised object.
ref_counts <- GetAssayData(seurat_ref, layer = "counts", assay = "RNA")
if (!all(ref_counts@x == round(ref_counts@x))) {
  stop("Reference 'counts' layer is not integer-valued. ",
       "Expected raw CellRanger counts (seurat_oligodendroglioma_RAWCOUNTS.rds).")
}
sc_matrix <- as.matrix(ref_counts)
message(sprintf("Reference matrix: %d genes x %d cells (raw integer counts)",
                nrow(sc_matrix), ncol(sc_matrix)))

# ==============================================================================
# 2. LOAD BULK RNA-seq
# ==============================================================================
message("\n[2] Loading bulk RNA-seq...")
bulk_df <- as.data.frame(read_excel(BULK_PATH))
rownames(bulk_df) <- bulk_df[[1]]
bulk_df <- bulk_df[, -1, drop = FALSE]

bulk_mat <- as.matrix(bulk_df)
storage.mode(bulk_mat) <- "numeric"
message(sprintf("Bulk matrix: %d genes x %d samples",
                nrow(bulk_mat), ncol(bulk_mat)))

# ==============================================================================
# 3. RESTRICT TO SHARED GENES
# ==============================================================================
message("\n[3] Intersecting reference and bulk gene sets...")
common_genes <- intersect(rownames(sc_matrix), rownames(bulk_mat))
message(sprintf("Shared genes: %d", length(common_genes)))

if (length(common_genes) < MIN_COMMON_GENE) {
  stop(sprintf("Only %d shared genes (< %d). ",
               length(common_genes), MIN_COMMON_GENE),
       "Check gene nomenclature (ENSEMBL vs SYMBOL) on both inputs.")
}

sc_shared   <- sc_matrix[common_genes, ]
bulk_shared <- bulk_mat[common_genes, ]

# ==============================================================================
# 4. BUILD AND RUN BAYESPRISM
# ==============================================================================
message("\n[4] Building BayesPrism object...")
cell_labels <- seurat_ref$cell_label

my_prism <- new.prism(
  reference         = t(sc_shared),    # cells x genes
  mixture           = t(bulk_shared),  # samples x genes
  input.type        = "count.matrix",
  cell.type.labels  = cell_labels,
  cell.state.labels = cell_labels,
  outlier.cut       = 0.01,
  outlier.fraction  = 0.1,
  key               = "Glioma"
)

message(sprintf("Running BayesPrism on %d cores...", N_CORES))
t_start <- Sys.time()

bp_res <- run.prism(prism = my_prism, n.cores = N_CORES)

message(sprintf("Deconvolution finished in %.1f min",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

# ==============================================================================
# 5. EXTRACT, VALIDATE, AND SAVE PROPORTIONS
# ==============================================================================
message("\n[5] Extracting proportion estimates...")

# Final theta: cell-type proportions after BayesPrism's refinement step.
theta <- get.fraction(
  bp            = bp_res,
  which.theta   = "final",
  state.or.type = "type"
)

message(sprintf("Value range: [%.4f, %.4f] | NA count: %d",
                min(theta, na.rm = TRUE), max(theta, na.rm = TRUE),
                sum(is.na(theta))))
message(sprintf("Per-sample proportion sums: [%.4f, %.4f] (expected ~1.0)",
                min(rowSums(theta)), max(rowSums(theta))))

write.csv(theta, file.path(OUT_DIR, "bayesprism_proportions.csv"))
message("Proportions written.")

# ==============================================================================
# 6. EXTRACT AND SAVE THE PURIFIED MALIGNANT TRANSCRIPTOME
# ------------------------------------------------------------------------------
# BayesPrism produces two distinct quantities from the same model:
#   - Purified expression (the Z tensor): the per-sample, per-cell-type
#     transcriptome estimated at cell-type resolution. It is exposed through
#     get.exp(state.or.type = "type"), which internally reads the cell-type
#     posterior; there is no separate "final" expression update, so this is
#     the expression output of record.
#   - Proportions (theta): refined in a subsequent step, hence "final" above.
# The two are therefore taken from their respective canonical outputs; they
# are not inconsistent.
# The purified Glioma expression (samples x genes) is used downstream to
# isolate the stemness signal to the malignant compartment.
# ==============================================================================
message("\n[6] Extracting purified Glioma expression...")

expr_glioma <- get.exp(
  bp            = bp_res,
  state.or.type = "type",
  cell.name     = "Glioma"
)
message(sprintf("Purified Glioma expression: %d samples x %d genes",
                nrow(expr_glioma), ncol(expr_glioma)))

# Persist the full result (contains purified expression for all cell types)
# and the Glioma matrix specifically, for the downstream scripts.
saveRDS(bp_res,      file.path(OUT_DIR, "bayesprism_full_result.rds"))
saveRDS(expr_glioma, file.path(OUT_DIR, "bayesprism_expression_glioma.rds"))
write.csv(expr_glioma, file.path(OUT_DIR, "bayesprism_expression_glioma.csv"))
message("Purified expression written.")

# ==============================================================================
# 7. FIGURES
# ==============================================================================
message("\n[7] Generating figures...")

theta_wide        <- as.data.frame(theta)
theta_wide$sample <- rownames(theta_wide)

theta_long <- reshape2::melt(
  theta_wide,
  id.vars       = "sample",
  variable.name = "cell_type",
  value.name    = "proportion"
)
theta_long$cell_type <- tfm_normalise_types(theta_long$cell_type)

# Order samples by malignant (Glioma) fraction for readability.
sample_order <- theta_wide[order(theta_wide$Glioma, decreasing = TRUE), "sample"]
theta_long$sample <- factor(theta_long$sample, levels = sample_order)

# -- Figure 1: stacked cell-type composition per sample -----------------------
p_composition <- ggplot(theta_long,
                        aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  tfm_scale_fill() +
  labs(
    title    = "Estimated cell-type composition (BayesPrism)",
    subtitle = "122 TCGA oligodendroglioma samples | Chu et al. 2022",
    x        = "Sample (ordered by malignant fraction)",
    y        = "Proportion",
    fill     = "Cell type"
  ) +
  theme_tfm()

ggsave(file.path(FIG_DIR, "bayesprism_composition.png"),
       p_composition, width = 14, height = 5, dpi = 300)

# -- Figure 2: per-cell-type proportion distribution --------------------------
p_box <- ggplot(theta_long,
                aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  tfm_scale_fill() +
  labs(
    title    = "Cell-type proportion distribution (BayesPrism)",
    subtitle = "122 TCGA oligodendroglioma samples",
    x        = "Cell type",
    y        = "Proportion"
  ) +
  theme_tfm(show_x_text = TRUE) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "bayesprism_boxplot.png"),
       p_box, width = 10, height = 5, dpi = 300)

message("Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 02 (BayesPrism) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Reference     : %d cells, %d cell types (LGG-04, raw counts)",
                ncol(seurat_ref), length(valid_types)))
message(sprintf("  Bulk cohort   : %d TCGA oligodendroglioma samples",
                ncol(bulk_shared)))
message(sprintf("  Shared genes  : %d", length(common_genes)))
message(sprintf("  Purified expr : %d samples x %d genes (Glioma)",
                nrow(expr_glioma), ncol(expr_glioma)))
message(sprintf("  Proportions   : %s",
                file.path(OUT_DIR, "bayesprism_proportions.csv")))
message(strrep("=", 70))