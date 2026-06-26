# ==============================================================================
# Standardised Colour Palette and ggplot2 Theme
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : utils/tfm_theme
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Single source of truth for figure styling across the whole pipeline.
#   Defines the cell-type colour palette, auxiliary palettes (deconvolution
#   methods, simulation strategies), a shared ggplot2 theme, and convenience
#   scale helpers. Sourcing this file guarantees that every figure in the
#   repository uses identical colours and styling.
#
# USAGE
#   source("scripts/utils/tfm_theme.R")
#   ggplot(df, aes(x, y, fill = cell_type)) +
#     geom_col() +
#     tfm_scale_fill() +
#     theme_tfm()
#
# NOTE ON CELL-TYPE NAMES
#   Palette names use dots (e.g. "s.mac.1", "CD8.TCells", "h.microglia"),
#   matching how R rewrites column names via make.names() after as.data.frame().
#   If a data frame still carries hyphen/space names (e.g. "s-mac 1"), normalise
#   them with tfm_normalise_types() before plotting so the colour mapping works.
#
# DEPENDENCIES
#   ggplot2 (functions are namespaced, so the package need not be attached).
# ==============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("tfm_theme.R requires the 'ggplot2' package to be installed.")
}

# ------------------------------------------------------------------------------
# 1. CELL-TYPE PALETTE
#    Stable colour per cell type, used in all composition and proportion plots.
# ------------------------------------------------------------------------------
tfm_colors <- c(
  "Glioma"        = "#D73027",
  "Oligo"         = "#1A9850",
  "Proliferating" = "#E7298A",
  "h.microglia"   = "#4575B4",
  "i.microglia"   = "#74ADD1",
  "a.microglia"   = "#ABD9E9",
  "AP.microglia"  = "#313695",
  "s.mac.1"       = "#542788",
  "s.mac.2"       = "#8073AC",
  "MDSC"          = "#B2ABD2",
  "Myeloid"       = "#D4B9DA",
  "TCells"        = "#F46D43",
  "CD8.TCells"    = "#FDAE61",
  "Endo"          = "#8C510A",
  "Pericytes"     = "#BF812D",
  "otherCells"    = "#BDBDBD"
)

# ------------------------------------------------------------------------------
# 2. AUXILIARY PALETTES
# ------------------------------------------------------------------------------
# Deconvolution methods (benchmark and method-comparison plots).
tfm_method_colors <- c(
  "MuSiC"      = "#2E75B6",
  "EPIC"       = "#27AE60",
  "BayesPrism" = "#E74C3C",
  "Consensus"  = "#8E44AD"
)

# Pseudobulk simulation strategies (benchmark plots).
tfm_strategy_colors <- c(
  "homo"  = "#74ADD1",
  "heter" = "#D73027"
)

# Broad biological categories (EPIC, run at 4-category resolution).
# Each category takes the representative colour of its lineage in tfm_colors.
tfm_category_colors <- c(
  "Tumour"     = "#D73027",  # Glioma lineage
  "Myeloid"    = "#4575B4",  # microglia/macrophage lineage
  "Lymphoid"   = "#F46D43",  # T-cell lineage
  "Stromal"    = "#8C510A",  # endothelial/pericyte lineage
  "otherCells" = "#BDBDBD"   # fraction not covered by the reference
)

# Molecular subtypes (Ceccarelli 2016): transcriptome classes used in the
# stemness comparison (script 12). NE = neural, PN = proneural; CL/ME included
# for completeness so plots that show all four classes also map correctly.
tfm_subtype_colors <- c(
  "NE" = "#4575B4",   # neural
  "PN" = "#D73027",   # proneural
  "CL" = "#1A9850",   # classical
  "ME" = "#FDAE61"    # mesenchymal
)

# ------------------------------------------------------------------------------
# 3. SHARED THEME
#    show_x_text = TRUE  -> rotated x-axis labels (e.g. categorical comparisons)
#    show_x_text = FALSE -> hidden x-axis labels (e.g. per-sample composition)
# ------------------------------------------------------------------------------
theme_tfm <- function(show_x_text = FALSE) {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = if (show_x_text)
        ggplot2::element_text(angle = 45, hjust = 1)
      else ggplot2::element_blank(),
      axis.ticks.x     = if (show_x_text)
        ggplot2::element_line()
      else ggplot2::element_blank(),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(size = 11, color = "grey30")
    )
}

# ------------------------------------------------------------------------------
# 4. CONVENIENCE SCALE HELPERS
#    Drop-in replacements for scale_*_manual using the project palettes.
#    drop = FALSE keeps absent cell types in the legend for cross-figure
#    consistency; set drop = TRUE to show only the types present.
# ------------------------------------------------------------------------------
tfm_scale_fill <- function(palette = tfm_colors, drop = FALSE, ...) {
  ggplot2::scale_fill_manual(values = palette, drop = drop, ...)
}

tfm_scale_color <- function(palette = tfm_colors, drop = FALSE, ...) {
  ggplot2::scale_color_manual(values = palette, drop = drop, ...)
}

# ------------------------------------------------------------------------------
# 5. NAME NORMALISATION HELPER
#    Converts cell-type labels with hyphens/spaces to the dotted convention
#    used by the palette, so colour mapping is robust to the input separator.
#    e.g. "s-mac 1" -> "s.mac.1", "CD8 TCells" -> "CD8.TCells"
# ------------------------------------------------------------------------------
tfm_normalise_types <- function(x) {
  make.names(as.character(x))
}