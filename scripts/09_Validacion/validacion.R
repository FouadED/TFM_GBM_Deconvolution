# ==============================================================================
# ROL2 Biological Validation of ROL1 Findings
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 09_Validacion
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Independent validation of the ROL1 consensus using reference-free (ROL2)
#   methods, each applied to its intended task. Tests whether the cellular
#   composition patterns replicate, and quantifies reference coverage.
#
# ANALYSES
#   A) xCell 2.0    -> NE vs PN Wilcoxon per cell category
#   B) MCP-counter  -> NE vs PN Wilcoxon per cell category
#   C) quanTIseq    -> NE vs PN Wilcoxon per cell category
#   D) ESTIMATE     -> Glioma% (ROL1) vs TumorPurity, Spearman
#   E) Bland-Altman -> agreement ROL1 vs ESTIMATE purity
#   F) EPIC otherCells -> single-patient reference coverage across the cohort
#
#   Multiple testing: Benjamini-Hochberg across all Wilcoxon tests.
#
# METHODOLOGICAL NOTES (logic, not fixed numbers; values are printed at runtime)
#   - Scales differ across methods (ssGSEA enrichment vs absolute fractions),
#     so direct numeric concordance is invalid; the replicated quantity is the
#     DIRECTION of the NE-vs-PN difference, not its magnitude (Sturm et al.,
#     2019). xCell T-cell scores are near-zero in both groups (immune-cold
#     tumour); MCP-counter is on a more interpretable scale. Both are read for
#     direction only.
#   - NK cells are the least reliable population to estimate (White et al.,
#     2024, DREAM challenge); methods disagree on direction, so NK is excluded
#     from the main interpretation.
#
# CIRCULARITY SAFEGUARD (Section E)
#   The Bland-Altman and Spearman analyses compare the UNCORRECTED consensus
#   (consensus_proportions.csv) against ESTIMATE. The purity-corrected consensus
#   must NOT be used here: it sets Glioma = ESTIMATE purity, which would make the
#   agreement trivially perfect. A negative bias (ROL1 - ESTIMATE) is expected
#   and structural: ROL1 separates malignant Glioma from normal Oligo, whereas
#   ESTIMATE groups both as "tumour". The Oligo-dominant samples drive the bias.
#
# EPIC otherCells (Section F)
#   otherCells = bulk fraction unexplained by the single-patient LGG-04
#   reference; its cohort median quantifies that reference's coverage limit.
#   The script reports, at runtime: whether the Oligo-dominant samples sit below
#   the cohort median (consistent with an OC-like malignant state, Tirosh 2016),
#   and whether NE differs from PN (LGG-04 being proneural-like is expected to
#   cover PN better).
#
# INPUTS
#   scripts/04_xCell/results/xcell2_scores.csv
#   scripts/05_MCP_Counter/results/mcpcounter_scores.csv
#   scripts/06_ESTIMATE/results/estimate_scores.csv
#   scripts/07_quanTIseq/results/quantiseq_proportions.csv
#   scripts/08_Consenso/results/consensus_proportions.csv   (UNCORRECTED)
#   scripts/03_EPIC/results/epic_reference_coverage.csv
#   scripts/12_Stemness/results/stemness_scores.csv         (clinical metadata)
#
# OUTPUTS
#   scripts/09_Validacion/figures/rol2_replication/*.png
#   scripts/09_Validacion/figures/rol2_replication/*.csv
#
# DEPENDENCY NOTE
#   Clinical subtype (NE/PN/Grade) is currently read from the stemness output,
#   so script 12 must run before this one. Consider moving clinical metadata to
#   a dedicated 00_Preprocessing output to decouple the two.
#
# REFERENCES
#   Bland & Altman (1986) Lancet; Sturm et al. (2019) Bioinformatics;
#   White et al. (2024) Nat Commun; Racle & Gfeller (2017) eLife;
#   Tirosh et al. (2016) Nature.
# ==============================================================================

set.seed(123)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(ggpubr)
})

# Shared palettes and theme (provides tfm_subtype_colors, theme_tfm, etc.).
source("scripts/utils/tfm_theme.R")

NE_COLOR <- tfm_subtype_colors[["NE"]]
PN_COLOR <- tfm_subtype_colors[["PN"]]

# ==============================================================================
# CONFIGURATION
# ==============================================================================
OUTPUT_DIR <- "scripts/09_Validacion/figures/rol2_replication"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EXPECTED_XCELL_TYPES <- 43  # BlueprintEncode panel; warn (not stop) if different
OLIGO_DOMINANT <- c("TCGA-DB-A64V", "TCGA-FG-8187", "TCGA-HT-7620",
                    "TCGA-HT-A4DV", "TCGA-HW-7486", "TCGA-QH-A6CU")

message(strrep("=", 60))
message(" ROL2 -- Independent biological validation")
message(strrep("=", 60))

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
normalize_sample_col <- function(df) {
  candidates <- c("sample_id", "sample", "Sample", "SampleID", "Sample_ID", "...1")
  sc <- intersect(candidates, names(df))[1]
  if (is.na(sc)) sc <- names(df)[1]
  if (sc != "sample_id") df <- df %>% rename(sample_id = !!sym(sc))
  df
}

truncate_tcga_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)

aggregate_category <- function(df, types) {
  present <- intersect(types, names(df))
  if (length(present) == 0) return(NULL)
  rowSums(df[, present, drop = FALSE], na.rm = TRUE)
}

# ==============================================================================
# LOAD DATA
# ==============================================================================
message("\n--- Loading data ---")

meta <- read_csv("scripts/12_Stemness/results/stemness_scores.csv",
                 show_col_types = FALSE) %>%
  rename(sample_id = sample_id_12, subtype = Transcriptome.Subtype) %>%
  filter(subtype %in% c("NE", "PN")) %>%
  dplyr::select(sample_id, subtype, Grade)
message(sprintf("  Clinical metadata: %d samples | NE=%d | PN=%d",
                nrow(meta), sum(meta$subtype == "NE"), sum(meta$subtype == "PN")))

load_scores <- function(path, label) {
  df <- read_csv(path, show_col_types = FALSE) %>%
    normalize_sample_col() %>%
    mutate(sample_id = truncate_tcga_id(sample_id))
  message(sprintf("  %s: %d samples x %d columns", label, nrow(df), ncol(df) - 1))
  df
}

xcell    <- load_scores("scripts/04_xCell/results/xcell2_scores.csv", "xCell 2.0")
if (ncol(xcell) - 1 != EXPECTED_XCELL_TYPES) {
  warning(sprintf("xcell2_scores.csv has %d cell types (expected %d). ",
                  ncol(xcell) - 1, EXPECTED_XCELL_TYPES),
          "Check the BlueprintEncode reference in 04_xCell.R.")
}
mcp       <- load_scores("scripts/05_MCP_Counter/results/mcpcounter_scores.csv", "MCP-counter")
quant     <- load_scores("scripts/07_quanTIseq/results/quantiseq_proportions.csv", "quanTIseq")
estimate  <- load_scores("scripts/06_ESTIMATE/results/estimate_scores.csv", "ESTIMATE")
consensus <- load_scores("scripts/08_Consenso/results/consensus_proportions.csv", "ROL1 consensus")

epic_cov_path <- "scripts/03_EPIC/results/epic_reference_coverage.csv"
epic_cov <- if (file.exists(epic_cov_path)) {
  load_scores(epic_cov_path, "EPIC otherCells")
} else {
  message("  WARNING: epic_reference_coverage.csv not found; Section F skipped.")
  NULL
}

# ==============================================================================
# CELL TYPE MAPPING (method-specific labels per broad category)
# ==============================================================================
MAPPING <- list(
  T_cells = list(
    xcell = c("CD4_positive_alpha_beta_T_cell", "CD8_positive_alpha_beta_T_cell",
              "regulatory_T_cell",
              "central_memory_CD4_positive_alpha_beta_T_cell",
              "effector_memory_CD4_positive_alpha_beta_T_cell",
              "central_memory_CD8_positive_alpha_beta_T_cell",
              "effector_memory_CD8_positive_alpha_beta_T_cell"),
    mcp   = c("T cells", "CD8 T cells", "Cytotoxic lymphocytes"),
    quant = c("T cell CD4+ (non-regulatory)", "T cell CD8+",
              "T cell regulatory (Tregs)")
  ),
  NK_cells = list(
    xcell = c("natural_killer_cell"),
    mcp   = c("NK cells"),
    quant = c("NK cell")
  ),
  Myeloid = list(
    xcell = c("monocyte", "macrophage", "inflammatory_macrophage",
              "alternatively_activated_macrophage", "dendritic_cell",
              "neutrophil", "granulocyte_monocyte_progenitor_cell"),
    mcp   = c("Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"),
    quant = c("Macrophage M1", "Macrophage M2", "Monocyte",
              "Neutrophil", "Myeloid dendritic cell")
  ),
  Endothelial = list(
    xcell = c("endothelial_cell", "microvascular_endothelial_cell"),
    mcp   = c("Endothelial cells"),
    quant = NULL
  )
)

# ==============================================================================
# WILCOXON HELPER (NE vs PN for one method/category)
# ==============================================================================
run_wilcoxon <- function(df_with_meta, mapping_types, method_name, cat_name) {
  na_row <- tibble(
    category = cat_name, method = method_name,
    n_NE = NA_integer_, n_PN = NA_integer_,
    median_NE = NA_real_, median_PN = NA_real_,
    delta_median = NA_real_, p_value = NA_real_, p_adj = NA_real_,
    direction = NA_character_, sig_label = "n.d.", n_types_used = 0L,
    note = "not estimated by this method")
  
  if (is.null(mapping_types) || length(mapping_types) == 0) return(na_row)
  
  v <- aggregate_category(
    df_with_meta %>% dplyr::select(-sample_id, -subtype), mapping_types)
  n_types <- length(intersect(mapping_types, names(df_with_meta)))
  if (is.null(v) || n_types == 0) {
    message(sprintf("  %-12s (%s): no matching columns", cat_name, method_name))
    return(na_row %>% mutate(note = "no matching columns found"))
  }
  
  ne  <- v[df_with_meta$subtype == "NE"]
  pn  <- v[df_with_meta$subtype == "PN"]
  wt  <- suppressWarnings(wilcox.test(ne, pn, exact = FALSE))
  dir <- ifelse(median(ne, na.rm = TRUE) > median(pn, na.rm = TRUE),
                "NE > PN", "PN > NE")
  message(sprintf("  %-12s (%s): %d types | %s | p=%.3e",
                  cat_name, method_name, n_types, dir, wt$p.value))
  
  tibble(
    category = cat_name, method = method_name,
    n_NE = length(ne), n_PN = length(pn),
    median_NE = round(median(ne, na.rm = TRUE), 6),
    median_PN = round(median(pn, na.rm = TRUE), 6),
    delta_median = round(median(ne, na.rm = TRUE) - median(pn, na.rm = TRUE), 6),
    p_value = wt$p.value, p_adj = NA_real_,
    direction = dir, sig_label = "pending",
    n_types_used = n_types, note = "")
}

# ==============================================================================
# ANALYSES A-C: NE vs PN per method
# ==============================================================================
attach_meta <- function(df) {
  df %>% inner_join(meta %>% dplyr::select(sample_id, subtype), by = "sample_id")
}

message("\n--- Analysis A: xCell 2.0 ---")
xcell_meta <- attach_meta(xcell)
xcell_results <- map_dfr(names(MAPPING),
                         ~ run_wilcoxon(xcell_meta, MAPPING[[.x]]$xcell, "xCell 2.0", .x))

message("\n--- Analysis B: MCP-counter ---")
mcp_meta <- attach_meta(mcp)
mcp_results <- map_dfr(names(MAPPING),
                       ~ run_wilcoxon(mcp_meta, MAPPING[[.x]]$mcp, "MCP-counter", .x))

message("\n--- Analysis C: quanTIseq ---")
quant_meta <- attach_meta(quant)
quant_results <- map_dfr(names(MAPPING),
                         ~ run_wilcoxon(quant_meta, MAPPING[[.x]]$quant, "quanTIseq", .x))

# ==============================================================================
# ANALYSIS D: ESTIMATE cross-validation (Spearman vs uncorrected consensus)
# ==============================================================================
message("\n--- Analysis D: ESTIMATE cross-validation ---")
est_glioma <- estimate %>%
  dplyr::select(sample_id, TumorPurity) %>%
  inner_join(consensus %>% dplyr::select(sample_id, Glioma), by = "sample_id") %>%
  dplyr::filter(!is.na(TumorPurity), !is.na(Glioma)) %>%
  left_join(meta %>% dplyr::select(sample_id, subtype), by = "sample_id")
ct <- suppressWarnings(
  cor.test(est_glioma$Glioma, est_glioma$TumorPurity, method = "spearman"))
message(sprintf("  n=%d | Spearman rho=%.3f | p=%.2e",
                nrow(est_glioma), ct$estimate, ct$p.value))

# ==============================================================================
# MULTIPLE TESTING CORRECTION (BH across all Wilcoxon tests)
# ==============================================================================
message("\n--- Multiple testing correction (BH) ---")
all_results <- bind_rows(xcell_results, mcp_results, quant_results)
valid_idx   <- !is.na(all_results$p_value)
all_results$p_adj[valid_idx] <- p.adjust(all_results$p_value[valid_idx], method = "BH")

all_results <- all_results %>%
  mutate(
    sig_label = case_when(
      sig_label == "n.d." | is.na(p_adj) ~ "n.d.",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "n.s."
    ),
    # Qualitative notes only; numeric effect sizes live in delta_median.
    note = case_when(
      category == "T_cells" & method == "xCell 2.0"   ~ "near-zero scores (immune-cold)",
      category == "T_cells" & method == "MCP-counter" ~ "interpretable scale",
      category == "NK_cells"                          ~ "discordant across methods -- excluded",
      TRUE ~ note
    )
  )

message("Corrected results (BH):")
print(all_results %>% dplyr::filter(!is.na(p_value)) %>%
        dplyr::select(method, category, direction, p_value, p_adj, sig_label) %>%
        dplyr::arrange(p_adj))
write_csv(all_results, file.path(OUTPUT_DIR, "rol2_replication_results.csv"))

# NK discordance note
nk_s <- all_results %>% dplyr::filter(category == "NK_cells", !is.na(p_value))
message("\n--- NK cells: discordance across methods (excluded from main) ---")
for (i in seq_len(nrow(nk_s))) {
  message(sprintf("  %-12s: %s | p_adj=%.4f",
                  nk_s$method[i], nk_s$direction[i], nk_s$p_adj[i]))
}

# ==============================================================================
# FIGURE 1: replication matrix
# ==============================================================================
message("\n--- Generating figures ---")
plot_data <- all_results %>%
  mutate(
    fill_val = case_when(
      sig_label %in% c("n.d.", "n.s.") ~ 0,
      direction == "NE > PN" ~  pmin(-log10(p_adj + 1e-15), 5),
      direction == "PN > NE" ~ -pmin(-log10(p_adj + 1e-15), 5),
      TRUE ~ 0
    ),
    cell_label = case_when(
      sig_label %in% c("n.d.", "n.s.") ~ sig_label,
      category == "NK_cells" ~ paste0(sig_label, "\n", gsub(" > ", ">", direction),
                                      "\n(discordant)"),
      TRUE ~ paste0(sig_label, "\n", gsub(" > ", ">", direction))
    ),
    category = factor(category, levels = rev(names(MAPPING))),
    method   = factor(method, levels = c("xCell 2.0", "MCP-counter", "quanTIseq"))
  )

p_matrix <- ggplot(plot_data, aes(x = method, y = category, fill = fill_val)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_label), size = 3.2, color = "black",
            fontface = "bold", lineheight = 0.85) +
  scale_fill_gradient2(low = PN_COLOR, mid = "#F5F5F5", high = NE_COLOR,
                       midpoint = 0, limits = c(-5, 5), na.value = "grey85",
                       name = "-log10(p_adj)\nblue = NE>PN\nred = PN>NE") +
  labs(
    title    = "Independent biological validation -- ROL2",
    subtitle = sprintf("Wilcoxon NE (n=%d) vs PN (n=%d) | BH correction (%d tests)",
                       sum(meta$subtype == "NE"), sum(meta$subtype == "PN"),
                       sum(!is.na(all_results$p_value))),
    x = NULL, y = NULL,
    caption = paste0("*** p_adj<0.001  ** p_adj<0.01  * p_adj<0.05  ",
                     "n.s. p_adj>=0.05  n.d. = not estimated\n",
                     "NK cells discordant across methods -- excluded (White et al. 2024)")
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(face = "bold", size = 11),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0),
        legend.key.height = unit(1.2, "cm"))
ggsave(file.path(OUTPUT_DIR, "fig1_ROL2_replication_matrix.png"),
       p_matrix, width = 9, height = 6, dpi = 300, bg = "white")
message("  fig1_ROL2_replication_matrix.png")

# ==============================================================================
# FIGURES 2A-2C: NE vs PN boxplots per method
#   stat_compare_means shows the RAW Wilcoxon p-value (the BH-adjusted values
#   are in fig1 / the results CSV); the caption states this explicitly.
# ==============================================================================
make_boxplots <- function(df_meta, mapping_list, method_label) {
  long_df <- map_dfr(names(mapping_list), function(cat_name) {
    types <- mapping_list[[cat_name]]
    if (is.null(types) || length(types) == 0) return(NULL)
    v <- aggregate_category(df_meta %>% dplyr::select(-sample_id, -subtype), types)
    if (is.null(v)) return(NULL)
    tibble(sample_id = df_meta$sample_id, subtype = df_meta$subtype,
           category = cat_name, score = v)
  }) %>% dplyr::filter(!is.na(score))
  if (nrow(long_df) == 0) return(invisible(NULL))
  
  ggplot(long_df, aes(x = subtype, y = score, fill = subtype)) +
    geom_boxplot(outlier.size = 1.0, alpha = 0.85, width = 0.55) +
    geom_jitter(width = 0.14, size = 0.6, alpha = 0.4) +
    scale_fill_manual(values = tfm_subtype_colors, guide = "none") +
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.y.npc = 0.92, size = 4.5) +
    facet_wrap(~ category, scales = "free_y", nrow = 1) +
    labs(title = paste0(method_label, " -- NE vs PN"), x = NULL,
         y = "Score / abundance",
         caption = "Wilcoxon rank-sum test (raw p shown; BH-adjusted values in fig1)") +
    theme_tfm(show_x_text = TRUE) +
    theme(strip.text = element_text(face = "bold"))
}

p2a <- make_boxplots(xcell_meta, lapply(MAPPING, `[[`, "xcell"),
                     "xCell 2.0 (Angel et al. 2025)")
if (!is.null(p2a)) ggsave(file.path(OUTPUT_DIR, "fig2a_xCell_NE_PN_boxplot.png"),
                          p2a, width = 13, height = 5, dpi = 300, bg = "white")
p2b <- make_boxplots(mcp_meta, lapply(MAPPING, `[[`, "mcp"),
                     "MCP-counter (Becht et al. 2016)")
if (!is.null(p2b)) ggsave(file.path(OUTPUT_DIR, "fig2b_MCP_NE_PN_boxplot.png"),
                          p2b, width = 13, height = 5, dpi = 300, bg = "white")
p2c <- make_boxplots(quant_meta, lapply(MAPPING, `[[`, "quant"),
                     "quanTIseq (Finotello et al. 2019)")
if (!is.null(p2c)) ggsave(file.path(OUTPUT_DIR, "fig2c_quanTIseq_NE_PN_boxplot.png"),
                          p2c, width = 13, height = 5, dpi = 300, bg = "white")
message("  fig2a/b/c boxplots")

# ==============================================================================
# FIGURE 3: ESTIMATE scatter (secondary; Bland-Altman is primary)
# ==============================================================================
p_scatter <- ggplot(est_glioma,
                    aes(x = Glioma * 100, y = TumorPurity * 100, color = subtype)) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "grey30",
              linewidth = 0.8, linetype = "dashed", formula = y ~ x) +
  scale_color_manual(values = tfm_subtype_colors, name = "Subtype",
                     na.value = "grey70") +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Spearman rho = %.2f\np = %.2e", ct$estimate, ct$p.value),
           fontface = "bold", size = 4.2) +
  labs(
    title    = "Cross-validation of tumour purity: ROL1 vs ESTIMATE",
    subtitle = "Glioma% consensus (ROL1) vs ESTIMATE purity (secondary figure)",
    x        = "Glioma fraction -- ROL1 consensus (%)",
    y        = "Tumour purity -- ESTIMATE (%)",
    caption  = sprintf(paste0("Spearman | n=%d | restricted purity range ",
                              "attenuates rho; see Bland-Altman (fig3b)."),
                       nrow(est_glioma))
  ) +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(OUTPUT_DIR, "fig3_ESTIMATE_purity_scatter.png"),
       p_scatter, width = 8, height = 6, dpi = 300, bg = "white")
message("  fig3_ESTIMATE_purity_scatter.png (secondary)")

# ==============================================================================
# ANALYSIS E: BLAND-ALTMAN (uncorrected consensus vs ESTIMATE) -- BRACKET
# ------------------------------------------------------------------------------
# Purity vs purity. The malignant fraction cannot be pinned to a single value
# because the deconvolution cannot separate OC-like malignant cells from normal
# oligodendrocytes (both share the oligodendroglial lineage program; Tirosh
# 2016). We therefore bracket it with two biologically motivated definitions:
#   LOWER bound  -> Glioma only        (excludes OC-like that fell into Oligo)
#   UPPER bound  -> Glioma + Oligo      (absorbs residual normal oligodendrocyte)
# Both are compared to ESTIMATE TumorPurity.
#
# Diagnostic (NOT a pipeline step; reported in the thesis discussion):
#   Wilcoxon TumorPurity oligo-dominant (n=6) vs rest (n=116): W=599, p=0.003;
#   medians 0.911 vs 0.958. The 6 Oligo-dominant samples are high-purity
#   (0.82-0.93) yet significantly less pure than the rest, consistent with an
#   Oligo fraction that is mostly OC-like malignant with a residual normal
#   oligodendrocyte contribution.
#
# Oligo-dominant samples are flagged from the explicit OLIGO_DOMINANT vector
# (Glioma < 0.22; 6 samples), consistent with Section F. The earlier
# Glioma<0.10 flag (1 sample) is replaced by this criterion.
# ==============================================================================
message("\n--- Analysis E: Bland-Altman agreement ROL1 vs ESTIMATE (bracket) ---")

# Build the base table once: both purity definitions + ESTIMATE + subtype.
ba_base <- data.frame(
  sid          = truncate_tcga_id(consensus[["sample_id"]]),
  glioma_only  = consensus[["Glioma"]],
  glioma_oligo = consensus[["Glioma"]] + consensus[["Oligo"]]) %>%
  inner_join(data.frame(sid = truncate_tcga_id(estimate[["sample_id"]]),
                        purity_estimate = estimate[["TumorPurity"]]),
             by = "sid") %>%
  left_join(meta %>% mutate(sid = truncate_tcga_id(sample_id)) %>%
              dplyr::select(sid, subtype), by = "sid") %>%
  mutate(is_oligo_dominant = sid %in% OLIGO_DOMINANT)

# Helper: Bland-Altman stats for one purity definition vs ESTIMATE.
bland_altman_stats <- function(purity_rol1, purity_estimate) {
  d    <- purity_rol1 - purity_estimate
  keep <- !is.na(d)
  d    <- d[keep]
  bias <- mean(d)
  s    <- sd(d)
  list(n = length(d), bias = bias, sd = s,
       loa_lower = bias - 1.96 * s, loa_upper = bias + 1.96 * s,
       pct_within = mean(d >= bias - 1.96 * s & d <= bias + 1.96 * s) * 100,
       mean_purity = (purity_rol1[keep] + purity_estimate[keep]) / 2,
       diff_purity = d)
}

ba_lower <- bland_altman_stats(ba_base$glioma_only,  ba_base$purity_estimate)
ba_upper <- bland_altman_stats(ba_base$glioma_oligo, ba_base$purity_estimate)

message(sprintf("  LOWER (Glioma only):  n=%d | bias=%.3f | LOA=[%.3f, %.3f] | within=%.1f%%",
                ba_lower$n, ba_lower$bias, ba_lower$loa_lower, ba_lower$loa_upper,
                ba_lower$pct_within))
message(sprintf("  UPPER (Glioma+Oligo): n=%d | bias=%.3f | LOA=[%.3f, %.3f] | within=%.1f%%",
                ba_upper$n, ba_upper$bias, ba_upper$loa_lower, ba_upper$loa_upper,
                ba_upper$pct_within))
message(sprintf("  Oligo-dominant samples flagged: %d", sum(ba_base$is_oligo_dominant)))

# ---- Long table for a faceted two-panel Bland-Altman ------------------------
ba_long <- bind_rows(
  data.frame(definition = "Lower bound: Glioma only",
             mean_purity = ba_lower$mean_purity, diff_purity = ba_lower$diff_purity,
             subtype = ba_base$subtype[!is.na(ba_base$glioma_only - ba_base$purity_estimate)],
             is_oligo_dominant = ba_base$is_oligo_dominant[!is.na(ba_base$glioma_only - ba_base$purity_estimate)]),
  data.frame(definition = "Upper bound: Glioma + Oligo",
             mean_purity = ba_upper$mean_purity, diff_purity = ba_upper$diff_purity,
             subtype = ba_base$subtype[!is.na(ba_base$glioma_oligo - ba_base$purity_estimate)],
             is_oligo_dominant = ba_base$is_oligo_dominant[!is.na(ba_base$glioma_oligo - ba_base$purity_estimate)]))
ba_long$definition <- factor(ba_long$definition,
                             levels = c("Lower bound: Glioma only",
                                        "Upper bound: Glioma + Oligo"))

# Per-panel bias / LOA lines.
ba_lines <- data.frame(
  definition = factor(c("Lower bound: Glioma only", "Upper bound: Glioma + Oligo"),
                      levels = levels(ba_long$definition)),
  bias      = c(ba_lower$bias,      ba_upper$bias),
  loa_upper = c(ba_lower$loa_upper, ba_upper$loa_upper),
  loa_lower = c(ba_lower$loa_lower, ba_upper$loa_lower))

p_ba <- ggplot(ba_long, aes(x = mean_purity, y = diff_purity)) +
  geom_point(aes(color = subtype, shape = is_oligo_dominant),
             alpha = 0.75, size = 2.2) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "Standard",
                                "TRUE" = "Oligo-dominant (Glioma<22%)"),
                     name = NULL) +
  geom_hline(data = ba_lines, aes(yintercept = bias),
             color = "#333333", linewidth = 0.9) +
  geom_hline(data = ba_lines, aes(yintercept = loa_upper),
             color = "#D73027", linetype = "dashed", linewidth = 0.7) +
  geom_hline(data = ba_lines, aes(yintercept = loa_lower),
             color = "#D73027", linetype = "dashed", linewidth = 0.7) +
  geom_text(data = ba_lines, aes(x = Inf, y = bias, label = sprintf("Bias = %.3f", bias)),
            hjust = 1.05, vjust = -0.5, size = 3.2, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = tfm_subtype_colors, name = "Subtype",
                     na.value = "grey70") +
  facet_wrap(~ definition, nrow = 1) +
  labs(
    title    = "Agreement between ROL1 and ESTIMATE tumour purity (bracket)",
    subtitle = sprintf("Bland-Altman | n=%d | lower bias=%.3f, upper bias=%.3f",
                       ba_lower$n, ba_lower$bias, ba_upper$bias),
    x = "Mean purity  (ROL1 + ESTIMATE) / 2",
    y = "Difference  ROL1 - ESTIMATE",
    caption = paste0("Solid: mean bias. Dashed: 95% limits of agreement ",
                     "(Bland & Altman 1986). Triangles: Oligo-dominant.\n",
                     "True malignant fraction lies between the two bounds: Glioma ",
                     "only excludes OC-like cells; Glioma+Oligo absorbs residual ",
                     "normal oligodendrocytes (Tirosh 2016).")
  ) +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(OUTPUT_DIR, "fig3b_bland_altman_rol1_vs_estimate.png"),
       p_ba, width = 13, height = 6, dpi = 300, bg = "white")
message("  fig3b_bland_altman_rol1_vs_estimate.png (primary, bracket)")

# ---- Stats table: both bounds + the purity diagnostic -----------------------
write.csv(
  data.frame(
    bound          = c("lower_glioma_only", "upper_glioma_oligo"),
    n              = c(ba_lower$n, ba_upper$n),
    bias           = round(c(ba_lower$bias, ba_upper$bias), 4),
    sd_diff        = round(c(ba_lower$sd, ba_upper$sd), 4),
    loa_lower      = round(c(ba_lower$loa_lower, ba_upper$loa_lower), 4),
    loa_upper      = round(c(ba_lower$loa_upper, ba_upper$loa_upper), 4),
    pct_within_loa = round(c(ba_lower$pct_within, ba_upper$pct_within), 1),
    n_oligo_dominant = sum(ba_base$is_oligo_dominant),
    spearman_rho   = round(as.numeric(ct$estimate), 3),
    spearman_p     = ct$p.value,
    reference      = "Bland & Altman 1986"),
  file.path(OUTPUT_DIR, "bland_altman_stats.csv"), row.names = FALSE)
# ==============================================================================
# ANALYSIS F: EPIC otherCells -- reference coverage
# ==============================================================================
message("\n--- Analysis F: EPIC reference coverage (otherCells) ---")
if (!is.null(epic_cov)) {
  meta_all <- read_csv("scripts/12_Stemness/results/stemness_scores.csv",
                       show_col_types = FALSE) %>%
    rename(sample_id = sample_id_12, subtype = Transcriptome.Subtype) %>%
    mutate(sample_id = truncate_tcga_id(sample_id)) %>%
    dplyr::select(sample_id, subtype)
  
  oc_df <- epic_cov %>%
    left_join(consensus %>% dplyr::select(sample_id, Glioma, Oligo), by = "sample_id") %>%
    left_join(meta_all, by = "sample_id") %>%
    mutate(coverage = 1 - otherCells, is_outlier = sample_id %in% OLIGO_DOMINANT)
  oc_median <- median(oc_df$otherCells, na.rm = TRUE)
  
  message(sprintf("  otherCells: median=%.3f | range=[%.3f, %.3f]",
                  oc_median, min(oc_df$otherCells), max(oc_df$otherCells)))
  message("  Oligo-dominant samples (expected below median):")
  print(oc_df %>% dplyr::filter(is_outlier) %>%
          dplyr::select(sample_id, Glioma, Oligo, otherCells, subtype) %>%
          dplyr::arrange(otherCells))
  
  oc_sub <- oc_df %>% dplyr::filter(subtype %in% c("NE", "PN"))
  wt_oc  <- wilcox.test(otherCells ~ subtype, data = oc_sub)
  ne_med <- median(oc_sub$otherCells[oc_sub$subtype == "NE"])
  pn_med <- median(oc_sub$otherCells[oc_sub$subtype == "PN"])
  message(sprintf("  NE median=%.3f | PN median=%.3f | Wilcoxon p=%.2e",
                  ne_med, pn_med, wt_oc$p.value))
  
  p_oc <- ggplot(oc_df, aes(x = Glioma, y = otherCells)) +
    geom_point(aes(color = subtype, shape = is_outlier), alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", color = "grey30", se = TRUE,
                linewidth = 0.8, formula = y ~ x) +
    geom_hline(yintercept = oc_median, linetype = "dashed", color = "grey50") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                       labels = c("FALSE" = "Standard", "TRUE" = "Oligo-dominant"),
                       name = NULL) +
    scale_color_manual(values = tfm_subtype_colors, name = "Subtype",
                       na.value = "grey70") +
    annotate("text", x = 0.05, y = oc_median + 0.025,
             label = sprintf("Cohort median = %.2f", oc_median),
             size = 3.2, color = "grey40", hjust = 0) +
    labs(
      title    = "EPIC reference coverage across the cohort",
      subtitle = sprintf("otherCells (unexplained fraction) | median=%.2f | LGG-04 single-patient reference",
                         oc_median),
      x = "Glioma fraction (ROL1 consensus)", y = "otherCells (EPIC)",
      caption = paste0("Dashed: cohort median. Triangles: Oligo-dominant samples. ",
                       "High overall otherCells quantifies the single-patient limit.\n",
                       "Reference: Racle & Gfeller 2017.")
    ) +
    theme_tfm(show_x_text = TRUE)
  ggsave(file.path(OUTPUT_DIR, "figF_EPIC_otherCells_coverage.png"),
         p_oc, width = 9, height = 6, dpi = 300, bg = "white")
  message("  figF_EPIC_otherCells_coverage.png")
  
  write.csv(
    oc_df %>% dplyr::select(sample_id, otherCells, coverage, Glioma, Oligo,
                            subtype, is_outlier),
    file.path(OUTPUT_DIR, "epic_otherCells_analysis.csv"), row.names = FALSE)
} else {
  message("  [skip] epic_reference_coverage.csv not found.")
}

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 60))
message(" SUMMARY -- ROL2 validation")
message(strrep("=", 60))
message(sprintf("  ESTIMATE Spearman: rho=%.3f | p=%.2e | n=%d",
                ct$estimate, ct$p.value, nrow(est_glioma)))
message(sprintf("  Bland-Altman: bias=%.3f | LOA=[%.3f, %.3f] | %.1f%% within",
                ba_bias, ba_loa_lower, ba_loa_upper, ba_pct))
if (!is.null(epic_cov)) {
  message(sprintf("  otherCells median=%.3f | NE=%.3f vs PN=%.3f",
                  oc_median, ne_med, pn_med))
}
message(sprintf("  Figures and tables in: %s", OUTPUT_DIR))
message(strrep("=", 60))