# ==============================================================================
# Reference Coverage and Coarse Deconvolution with EPIC
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 03_EPIC
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Run EPIC as a complementary, NON-consensus method providing two outputs:
#     (1) otherCells -- the transcriptional fraction not explained by the
#         reference, used as a reference-coverage / confidence metric;
#     (2) coarse 4-category proportions, used only as supplementary
#         concordance evidence against the ROL1 consensus.
#   EPIC proportions are NOT included in the consensus (see RATIONALE).
#
# RATIONALE (why 4 categories, not 15 cell types)
#   EPIC solves constrained least squares and was validated for ~5-7 reference
#   cell types (Racle & Gfeller, 2017). With the full 15-type reference the
#   eight myeloid subtypes are strongly collinear, making the reference matrix
#   near-singular and the per-type fractions unstable. EPIC is therefore run
#   at a 4-category resolution (Tumour / Myeloid / Lymphoid / Stromal) in BOTH
#   the pipeline and the benchmark, ensuring methodological consistency. The
#   15-type consensus is built from MuSiC and BayesPrism only.
#
# NORMALIZATION
#   EPIC requires TPM (Racle & Gfeller, 2017, eq. 2). Bulk and reference are
#   both TPM-normalised using gene lengths retrieved from Ensembl via biomaRt
#   and cached locally so offline re-runs (e.g. on the HPC cluster) need no
#   internet access.
#
# INPUTS
#   data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds
#       Annotated scRNA-seq reference with RAW CellRanger counts (GSE182109).
#   data/bulk/oligo_counts_protein.coding.xlsx
#       Bulk RNA-seq raw counts, 122 TCGA samples (protein-coding genes).
#
# OUTPUTS
#   scripts/03_EPIC/results/epic_reference_coverage.csv      (otherCells; main)
#   scripts/03_EPIC/results/epic_proportions_4category.csv   (supplementary)
#   scripts/03_EPIC/results/gene_lengths_cache.csv           (reused in re-runs)
#   scripts/03_EPIC/figures/*.png
#
# REFERENCES
#   Racle & Gfeller (2017)    EPIC. eLife 6:e26476.
#   Finotello & Trajanoski (2018) Immune deconvolution review. Cancer Immunol
#                             Immunother 67:1031-1040.
#   Avila Cobos et al. (2020) Deconvolution benchmarking. Nat Commun 11:5650.
# ==============================================================================

set.seed(123)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SEURAT_PATH       <- "data/scrna/seurat_oligodendroglioma_RAWCOUNTS.rds"
BULK_PATH         <- "data/bulk/oligo_counts_protein.coding.xlsx"
OUT_DIR           <- "scripts/03_EPIC/results"
FIG_DIR           <- "scripts/03_EPIC/figures"
GENE_LENGTH_CACHE <- file.path(OUT_DIR, "gene_lengths_cache.csv")

MIN_CELLS         <- 20      # minimum cells per type to build a stable profile
MIN_COMMON_GENE   <- 1000    # abort if fewer shared genes
MIN_GENES_FOR_TPM <- 5000    # minimum genes with length info to use TPM
OTHERCELLS_THRESH <- 0.35    # annotation threshold for the coverage figure

# Broad biological categories. EPIC is run at this resolution to avoid the
# near-singular reference produced by the eight collinear myeloid subtypes.
# Names use the hyphen/space convention of the SubAssignment labels.
TYPE_GROUPS <- list(
  Tumour   = c("Glioma", "Oligo", "Proliferating"),
  Myeloid  = c("a-microglia", "AP-microglia", "h-microglia", "i-microglia",
               "MDSC", "Myeloid", "s-mac 1", "s-mac 2"),
  Lymphoid = c("TCells", "CD8 TCells"),
  Stromal  = c("Endo", "Pericytes")
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
message("Loading libraries...")
suppressPackageStartupMessages({
  library(EPIC)
  library(Seurat)
  library(readxl)
  library(ggplot2)
  library(reshape2)
  library(biomaRt)
  library(dplyr)
})

# Shared figure styling: palettes, theme and scale helpers.
source("scripts/utils/tfm_theme.R")

# ==============================================================================
# HELPER: TPM normalisation
#   TPM = (counts / length_kb) / sum(counts / length_kb) * 1e6
#   Bulk and reference must share the same normalisation (Racle & Gfeller 2017).
# ==============================================================================
compute_tpm <- function(count_matrix, gene_lengths_bp) {
  shared <- intersect(rownames(count_matrix), names(gene_lengths_bp))
  if (length(shared) < 100) {
    stop("Fewer than 100 genes with length information. Check gene names.")
  }
  mat   <- count_matrix[shared, , drop = FALSE]
  gl_kb <- gene_lengths_bp[shared] / 1000
  rpk   <- mat / gl_kb
  t(t(rpk) / colSums(rpk) * 1e6)
}

# ==============================================================================
# 1. LOAD BULK RNA-seq
# ==============================================================================
message("\n[1] Loading bulk RNA-seq...")
bulk_df <- as.data.frame(read_excel(BULK_PATH))
rownames(bulk_df) <- bulk_df[[1]]
bulk_df <- bulk_df[, -1, drop = FALSE]
bulk_mat <- as.matrix(bulk_df)
storage.mode(bulk_mat) <- "numeric"
message(sprintf("Bulk matrix: %d genes x %d samples",
                nrow(bulk_mat), ncol(bulk_mat)))

# ==============================================================================
# 2. GENE LENGTHS FROM ENSEMBL (cached for offline re-runs)
# ==============================================================================
gene_lengths <- NULL
if (file.exists(GENE_LENGTH_CACHE)) {
  message(sprintf("[2] Loading gene lengths from cache: %s", GENE_LENGTH_CACHE))
  gl_df        <- read.csv(GENE_LENGTH_CACHE)
  gene_lengths <- setNames(gl_df$median_length_bp, gl_df$hgnc_symbol)
  message(sprintf("    Cached lengths: %d genes", length(gene_lengths)))
} else {
  message("[2] Querying Ensembl for gene lengths (one-time, ~2 min)...")
  tryCatch({
    mart   <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gl_raw <- getBM(
      attributes = c("hgnc_symbol", "transcript_length"),
      filters    = "hgnc_symbol",
      values     = unique(rownames(bulk_mat)),
      mart       = mart
    )
    gl_df <- gl_raw %>%
      dplyr::filter(hgnc_symbol != "" & transcript_length > 0) %>%
      dplyr::group_by(hgnc_symbol) %>%
      dplyr::summarise(median_length_bp = median(transcript_length),
                       .groups = "drop") %>%
      as.data.frame()
    write.csv(gl_df, GENE_LENGTH_CACHE, row.names = FALSE)
    gene_lengths <- setNames(gl_df$median_length_bp, gl_df$hgnc_symbol)
    message(sprintf("    Retrieved %d gene lengths; cache written.",
                    length(gene_lengths)))
  }, error = function(e) {
    message("    WARNING: biomaRt query failed: ", conditionMessage(e))
    message("    Falling back to CPM. For TPM, obtain gene_lengths_cache.csv ",
            "on a machine with internet and place it in ", OUT_DIR)
  })
}

# ==============================================================================
# 3. NORMALISE BULK TO TPM (CPM fallback if lengths unavailable)
# ==============================================================================
message("\n[3] Normalising bulk...")
genes_with_length <- intersect(rownames(bulk_mat), names(gene_lengths))

if (!is.null(gene_lengths) && length(genes_with_length) >= MIN_GENES_FOR_TPM) {
  bulk_norm   <- compute_tpm(bulk_mat, gene_lengths)
  norm_method <- "TPM"
  message(sprintf("    TPM normalisation (%d genes with length info).",
                  length(genes_with_length)))
} else {
  bulk_norm   <- t(t(bulk_mat) / colSums(bulk_mat) * 1e6)
  norm_method <- "CPM (fallback)"
  message("    WARNING: using CPM fallback; deviates from EPIC's TPM ",
          "requirement (Racle & Gfeller 2017).")
}

# ==============================================================================
# 4. LOAD REFERENCE AND BUILD 4-CATEGORY PROFILES
# ==============================================================================
message("\n[4] Building 4-category reference...")
seurat_obj <- readRDS(SEURAT_PATH)

assignment     <- as.character(seurat_obj@meta.data$Assignment)
sub_assignment <- as.character(seurat_obj@meta.data$SubAssignment)
seurat_obj$cell_label <- ifelse(
  is.na(sub_assignment) | sub_assignment %in% c("", "NA"),
  assignment, sub_assignment
)

type_counts <- table(seurat_obj$cell_label)
valid_types <- names(type_counts[type_counts >= MIN_CELLS])
seurat_ref  <- seurat_obj[, seurat_obj$cell_label %in% valid_types]

# Sanity check: reference must contain raw integer counts.
ref_counts_cells <- GetAssayData(seurat_ref, layer = "counts", assay = "RNA")
if (!all(ref_counts_cells@x == round(ref_counts_cells@x))) {
  stop("Reference 'counts' layer is not integer-valued. ",
       "Expected raw CellRanger counts (seurat_oligodendroglioma_RAWCOUNTS.rds).")
}
sc_counts   <- as.matrix(ref_counts_cells)
cell_labels <- seurat_ref$cell_label

# Mean count profile per broad category (only types present in the reference).
ref_counts_grouped <- sapply(names(TYPE_GROUPS), function(grp) {
  types_in_grp <- intersect(TYPE_GROUPS[[grp]], valid_types)
  idx          <- which(cell_labels %in% types_in_grp)
  if (length(idx) == 0) return(rep(0, nrow(sc_counts)))
  rowMeans(sc_counts[, idx, drop = FALSE])
})
message(sprintf("    Categories: %s",
                paste(colnames(ref_counts_grouped), collapse = ", ")))

# Normalise reference profiles to the same space as the bulk.
if (norm_method == "TPM") {
  ref_norm <- compute_tpm(ref_counts_grouped, gene_lengths)
} else {
  ref_norm <- t(t(ref_counts_grouped) / colSums(ref_counts_grouped) * 1e6)
}

# ==============================================================================
# 5. RESTRICT TO SHARED GENES
# ==============================================================================
message("\n[5] Intersecting reference and bulk gene sets...")
common_genes <- intersect(rownames(ref_norm), rownames(bulk_norm))
message(sprintf("    Shared genes: %d", length(common_genes)))
if (length(common_genes) < MIN_COMMON_GENE) {
  stop(sprintf("Only %d shared genes (< %d). Check gene nomenclature.",
               length(common_genes), MIN_COMMON_GENE))
}
ref_shared  <- ref_norm[common_genes, ]
bulk_shared <- bulk_norm[common_genes, ]

# ==============================================================================
# 6. DECONVOLUTION WITH EPIC
#   withOther = TRUE  -> estimate otherCells (uncovered fraction)
#   scaleExprs = FALSE -> inputs are already on a common (TPM) scale
# ==============================================================================
message("\n[6] Running EPIC (4 categories)...")
epic_ref <- list(
  refProfiles = ref_shared,
  sigGenes    = rownames(ref_shared)
)

epic_result <- EPIC(
  bulk       = bulk_shared,
  reference  = epic_ref,
  withOther  = TRUE,
  scaleExprs = FALSE
)

# ==============================================================================
# 7. otherCells -> REFERENCE COVERAGE METRIC (MAIN OUTPUT)
# ==============================================================================
message("\n[7] Extracting otherCells coverage metric...")
epic_fractions <- as.data.frame(epic_result$cellFractions)

coverage_df <- data.frame(
  sample_id   = rownames(epic_fractions),
  otherCells  = epic_fractions$otherCells,
  norm_method = norm_method
)
message("    otherCells summary:")
print(summary(coverage_df$otherCells))
message(sprintf("    Samples > %.0f%% otherCells: %d / %d",
                100 * OTHERCELLS_THRESH,
                sum(coverage_df$otherCells > OTHERCELLS_THRESH),
                nrow(coverage_df)))

write.csv(coverage_df,
          file.path(OUT_DIR, "epic_reference_coverage.csv"),
          row.names = FALSE)
message("    Coverage metric written (main output).")

# ==============================================================================
# 8. 4-CATEGORY PROPORTIONS -> SUPPLEMENTARY (NOT consensus)
# ==============================================================================
message("\n[8] Extracting 4-category proportions (supplementary)...")
prop_cols  <- setdiff(colnames(epic_fractions), "otherCells")
prop_clean <- epic_fractions[, prop_cols, drop = FALSE]

# Renormalise across categories (excluding otherCells) for concordance plots.
prop_renorm <- t(apply(prop_clean, 1, function(r) {
  s <- sum(r, na.rm = TRUE); if (s == 0) r else r / s
}))
message(sprintf("    Per-sample sums after renormalisation: [%.4f, %.4f]",
                min(rowSums(prop_renorm)), max(rowSums(prop_renorm))))

write.csv(prop_renorm,
          file.path(OUT_DIR, "epic_proportions_4category.csv"))
message("    4-category proportions written (supplementary).")

# ==============================================================================
# 9. FIGURES
# ==============================================================================
message("\n[9] Generating figures...")

# -- Figure 1: otherCells distribution (coverage; the main EPIC result) -------
p_coverage <- ggplot(coverage_df, aes(x = otherCells)) +
  geom_histogram(bins = 30, fill = tfm_category_colors[["otherCells"]],
                 color = "white", alpha = 0.9) +
  geom_vline(xintercept = OTHERCELLS_THRESH, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  annotate("text", x = OTHERCELLS_THRESH + 0.02, y = Inf, vjust = 1.5,
           label = sprintf("%.0f%% threshold", 100 * OTHERCELLS_THRESH),
           color = "#D73027", size = 3.5) +
  labs(
    title    = "EPIC reference coverage (otherCells)",
    subtitle = paste0("Fraction not explained by the LGG-04 reference | ",
                      norm_method),
    x        = "otherCells fraction",
    y        = "Number of samples"
  ) +
  theme_tfm(show_x_text = TRUE)

ggsave(file.path(FIG_DIR, "epic_othercells_coverage.png"),
       p_coverage, width = 8, height = 5, dpi = 300)

# -- Supplementary figures: 4-category proportions ----------------------------
prop_wide        <- as.data.frame(prop_renorm)
prop_wide$sample <- rownames(prop_wide)
prop_long <- reshape2::melt(prop_wide, id.vars = "sample",
                            variable.name = "category",
                            value.name = "proportion")

if ("Tumour" %in% colnames(prop_wide)) {
  sample_order <- prop_wide[order(prop_wide$Tumour, decreasing = TRUE), "sample"]
  prop_long$sample <- factor(prop_long$sample, levels = sample_order)
}

p_composition <- ggplot(prop_long,
                        aes(x = sample, y = proportion, fill = category)) +
  geom_col(width = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_fill_manual(values = tfm_category_colors) +
  labs(
    title    = "EPIC coarse composition (4 categories, supplementary)",
    subtitle = paste0("122 TCGA oligodendroglioma samples | ", norm_method,
                      " | otherCells excluded and renormalised"),
    x        = "Sample (ordered by Tumour fraction)",
    y        = "Proportion", fill = "Category"
  ) +
  theme_tfm()

ggsave(file.path(FIG_DIR, "epic_composition_4category.png"),
       p_composition, width = 14, height = 5, dpi = 300)

p_box <- ggplot(prop_long,
                aes(x = category, y = proportion, fill = category)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  scale_fill_manual(values = tfm_category_colors) +
  labs(
    title    = "EPIC category proportion distribution (supplementary)",
    subtitle = "122 TCGA oligodendroglioma samples (otherCells excluded)",
    x        = "Category", y = "Proportion"
  ) +
  theme_tfm(show_x_text = TRUE) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "epic_boxplot_4category.png"),
       p_box, width = 8, height = 5, dpi = 300)

message("    Figures written.")

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 70))
message("  SCRIPT 03 (EPIC) COMPLETED")
message(strrep("=", 70))
message(sprintf("  Role          : reference coverage (otherCells) + ",
                "supplementary 4-category proportions"))
message(sprintf("  Normalisation : %s", norm_method))
message(sprintf("  Categories    : %d (Tumour/Myeloid/Lymphoid/Stromal)",
                ncol(ref_shared)))
message(sprintf("  Shared genes  : %d", length(common_genes)))
message(sprintf("  otherCells med: %.3f", median(coverage_df$otherCells)))
message(sprintf("  Main output   : %s",
                file.path(OUT_DIR, "epic_reference_coverage.csv")))
message(strrep("=", 70))