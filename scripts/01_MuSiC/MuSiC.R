# ==============================================================================
# Cell-Type Deconvolution with MuSiC
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 01_MuSiC
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Estimate cell-type proportions in 122 TCGA oligodendroglioma bulk RNA-seq
#   samples using MuSiC, with a patient-matched single-cell RNA-seq reference
#   (LGG-04; Abdelfattah et al., 2022).
#
# METHOD
#   music_prop() is used (not music2_prop): the reference and the bulk cohort
#   share the same clinical condition (IDH-mutant, 1p/19q-codeleted), so the
#   multi-condition modelling of MuSiC2 is unnecessary (Fan et al., 2022).
#   MuSiC weights genes by cross-subject consistency via W-NNLS. Here the three
#   multiregional tumour fragments of LGG-04 serve as the subject grouping;
#   the resulting cross-subject variance is therefore intra-tumoral
#   (between regions) rather than inter-patient. This is a known limitation of
#   a single-donor reference and is reported as such in the manuscript.
#
# INPUTS
#   data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds
#       Annotated scRNA-seq reference with RAW CellRanger counts (GSE182109).
#   data/bulk/oligo_counts_protein.coding.xlsx
#       Bulk RNA-seq raw counts, 122 TCGA samples (protein-coding genes).
#
# OUTPUTS
#   scripts/01_MuSiC/results/music_proportions.csv
#   scripts/01_MuSiC/figures/music_composition.png
#
# REPRODUCIBILITY
#   Run from the project root (setwd to the repository root before sourcing).
#   Deterministic given a fixed seed; see sessionInfo() for package versions.
#
# REFERENCES
#   Wang et al. (2019)        MuSiC. Nat Commun 10:380.
#   Abdelfattah et al. (2022) scRNA-seq glioma atlas. Nat Commun 13:767.
#   Fan et al. (2022)         MuSiC2. Brief Bioinform 23:bbac430.
#   Avila Cobos et al. (2020) Deconvolution benchmarking. Nat Commun 11:5650.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SEURAT_PATH <- "data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds"
BULK_PATH   <- "data/bulk/oligo_counts_protein.coding.xlsx"
OUT_DIR     <- "scripts/01_MuSiC/results"
FIG_DIR     <- "scripts/01_MuSiC/figures"

MIN_CELLS       <- 20    # minimum cells per type to build a stable profile
MIN_COMMON_GENE <- 1000  # abort if fewer shared genes (nomenclature mismatch)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(MuSiC)
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

# Resolve the working cell-type label: prefer the fine-grained SubAssignment,
# falling back to the coarse Assignment when the subtype is not annotated.
assignment     <- as.character(seurat_obj@meta.data$Assignment)
sub_assignment <- as.character(seurat_obj@meta.data$SubAssignment)
seurat_obj$cell_label <- ifelse(
  is.na(sub_assignment) | sub_assignment %in% c("", "NA"),
  assignment, sub_assignment
)

message("Cell-type counts (pre-filter):")
print(table(seurat_obj$cell_label))

# Retain only cell types with enough cells to estimate a reliable profile.
# Below this threshold the per-type mean expression is statistically unstable
# (Avila Cobos et al., 2020).
type_counts <- table(seurat_obj$cell_label)
valid_types <- names(type_counts[type_counts >= MIN_CELLS])

message(sprintf("\nValid cell types (>= %d cells): %d",
                MIN_CELLS, length(valid_types)))
message("Excluded (too few cells):")
print(type_counts[type_counts < MIN_CELLS])

seurat_ref <- seurat_obj[, seurat_obj$cell_label %in% valid_types]
message("\nCells per tumour fragment:")
print(table(seurat_ref$orig.ident))

# Sanity check: the reference must contain RAW integer counts.
# This guards against accidentally loading a log-normalised object, which
# would silently corrupt the deconvolution (MuSiC expects counts).
ref_counts <- GetAssayData(seurat_ref, layer = "counts", assay = "RNA")
if (!all(ref_counts@x == round(ref_counts@x))) {
  stop("Reference 'counts' layer is not integer-valued. ",
       "Expected raw CellRanger counts (seurat_oligodendroglioma_RAWCOUNTS.rds).")
}
message("Reference count matrix verified as raw integer counts.")

# Build the SingleCellExperiment consumed by MuSiC.
# sampleID = tumour fragment -> defines the subject grouping for W-NNLS.
sce <- SingleCellExperiment(
  assays  = list(counts = ref_counts),
  colData = DataFrame(
    cellType = seurat_ref$cell_label,
    sampleID = seurat_ref$orig.ident
  )
)
message(sprintf("SCE built: %d cells x %d genes", ncol(sce), nrow(sce)))

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
common_genes <- intersect(rownames(sce), rownames(bulk_mat))
message(sprintf("Shared genes: %d", length(common_genes)))

if (length(common_genes) < MIN_COMMON_GENE) {
  stop(sprintf("Only %d shared genes (< %d). ",
               length(common_genes), MIN_COMMON_GENE),
       "Check gene nomenclature (ENSEMBL vs SYMBOL) on both inputs.")
}

sce_shared  <- sce[common_genes, ]
bulk_shared <- bulk_mat[common_genes, ]

# ==============================================================================
# 4. DECONVOLUTION
# ==============================================================================
message("\n[4] Running MuSiC (music_prop)...")
t_start <- Sys.time()

music_result <- music_prop(
  bulk.mtx  = bulk_shared,
  sc.sce    = sce_shared,
  clusters  = "cellType",
  samples   = "sampleID",
  select.ct = valid_types,
  verbose   = FALSE
)

message(sprintf("Deconvolution finished in %.1f s",
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

# ==============================================================================
# 5. EXTRACT, VALIDATE, AND SAVE PROPORTIONS
# ==============================================================================
message("\n[5] Extracting proportion estimates...")

# music_prop returns weighted (W-NNLS) and unweighted (NNLS) estimates;
# the weighted estimate is the MuSiC output of record.
proportions <- if ("Est.prop.weighted" %in% names(music_result)) {
  music_result$Est.prop.weighted
} else if ("Est.prop" %in% names(music_result)) {
  music_result$Est.prop
} else {
  music_result[[1]]
}

message("Estimated proportions (first 3 samples):")
print(round(head(proportions, 3), 4))

# Quality checks
message(sprintf("Value range: [%.4f, %.4f] | NA count: %d",
                min(proportions, na.rm = TRUE),
                max(proportions, na.rm = TRUE),
                sum(is.na(proportions))))

row_sums <- rowSums(proportions, na.rm = TRUE)
message(sprintf("Per-sample proportion sums: [%.4f, %.4f] (expected ~1.0)",
                min(row_sums), max(row_sums)))

out_csv <- file.path(OUT_DIR, "music_proportions.csv")
write.csv(proportions, out_csv)
message(sprintf("Proportions written to: %s", out_csv))

# ==============================================================================
# 6. FIGURE: STACKED CELL-TYPE COMPOSITION PER SAMPLE
# ==============================================================================
message("\n[6] Generating composition figure...")

prop_wide        <- as.data.frame(proportions)
prop_wide$sample <- rownames(prop_wide)

prop_long <- reshape2::melt(
  prop_wide,
  id.vars       = "sample",
  variable.name = "cell_type",
  value.name    = "proportion"
)

# Normalise cell-type labels to the palette convention (dotted names),
# so the colour mapping in tfm_colors matches regardless of the separator.
prop_long$cell_type <- tfm_normalise_types(prop_long$cell_type)

# Order samples by malignant (Glioma) fraction for readability.
sample_order <- prop_wide[order(prop_wide$Glioma, decreasing = TRUE), "sample"]
prop_long$sample <- factor(prop_long$sample, levels = sample_order)

p_composition <- ggplot(prop_long,
                        aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col(width = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(
    title    = "Estimated cell-type composition (MuSiC)",
    subtitle = "122 TCGA oligodendroglioma samples | music_prop",
    x        = "Sample (ordered by malignant fraction)",
    y        = "Proportion",
    fill     = "Cell type"
  ) +
  tfm_scale_fill() +
  theme_tfm()

fig_path <- file.path(FIG_DIR, "music_composition.png")
ggsave(fig_path, p_composition, width = 14, height = 6, dpi = 300)
message(sprintf("Figure written to: %s", fig_path))

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 01 (MuSiC) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Method        : music_prop (W-NNLS)"))
message(sprintf("  Reference     : %d cells, %d cell types (LGG-04, raw counts)",
                ncol(seurat_ref), length(valid_types)))
message(sprintf("  Bulk cohort   : %d TCGA oligodendroglioma samples",
                ncol(bulk_shared)))
message(sprintf("  Shared genes  : %d", length(common_genes)))
message(sprintf("  Proportions   : %s", out_csv))
message(strrep("=", 70))

