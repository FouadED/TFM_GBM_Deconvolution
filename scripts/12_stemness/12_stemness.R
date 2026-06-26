# ==============================================================================
# Glioma Stemness Analysis (with sensitivity analysis)
# ------------------------------------------------------------------------------
# Project : Reproducible cell-type deconvolution pipeline for IDH-mutant
#           1p/19q-codeleted oligodendroglioma (bulk RNA-seq)
# Script  : 12_Stemness
# Author  : Fouad Eddaoudi Lakraichi
# ------------------------------------------------------------------------------
# PURPOSE
#   Characterise stemness at cell-type resolution and test whether the proneural
#   (PN) subtype is more stem-like than the neural (NE) subtype. Stemness is
#   scored on the BayesPrism-purified Glioma expression (Z matrix), removing
#   immune/stromal contamination, rather than via the pan-cancer mRNAsi (Malta
#   et al., 2018), which needs restricted OCLR weights.
#
# STEMNESS PROXIES
#   1. Proliferating-cell proportion (consensus deconvolution; Tirosh 2016).
#   2. Curated GSC + oligodendroglioma stem-like signature score on the Z matrix
#      (primary proxy; see SIGNATURE).
#   3. Survival association (SUPPLEMENTARY ONLY -- underpowered; see caveat).
#
# SIGNATURE (curated, justified per gene; not a validated pan-cancer index)
#   GSC core regulators : SOX2, POU3F2, SALL2  (Suva et al. 2014, Cell)
#   Self-renewal        : MYC, MYCN            (Wang 2008; Tateishi 2016)
#   Oligo stem-like     : SOX4, SOX11, CCND2, CTNNB1, USP22, MSI1
#                         (Tirosh et al. 2016, Nature, Ext. Data Fig. 3b)
#   Excluded: OLIG1/OLIG2 (OC-like / differentiated pole in Tirosh 2016);
#   LIN28A/B (absent from the Z matrix).
#   Scoring: per-gene z-score across samples, then per-sample mean of the
#   signature z-scores (compute_stemness, used for both the main score and the
#   sensitivity control).
#
# PRIMARY TEST (Step 6)
#   Wilcoxon rank-sum on the stemness score, NE vs PN. Malta et al. (2018) is
#   cited only for the biological interpretation (stemness ~ immune evasion).
#
# SENSITIVITY ANALYSIS (Step 7) -- circularity guard
#   The main score uses the purified Z matrix, which comes from the
#   deconvolution. Step 7 recomputes the SAME score on the raw bulk (TPM), which
#   is deconvolution-independent. If PN > NE holds on both and the two scores
#   correlate, the result is robust. (This controls for the deconvolution; it
#   does not remove the subtler point that NE/PN are themselves transcriptionally
#   defined, which is addressed in the discussion.)
#
# SURVIVAL CAVEAT (Step 8)
#   Underpowered (events per variable ~5, below the EPV>=10 rule of thumb);
#   reported as SUPPLEMENTARY only and retained for reproducibility.
#
# CLINICAL DATA (no dependency on script 09)
#   Obtained from TCGAquery_subtype and cached locally (clinical_cache.csv), so
#   the script runs independently and offline after the first run. This breaks
#   the previous circular dependency with script 09.
#
# INPUTS
#   scripts/08_Consenso/results/consensus_proportions.csv
#   scripts/08_Consenso/results/consensus_proportions_purity_corrected.csv
#   scripts/02_BayesPrism/results/bayesprism_expression_glioma.rds
#   scripts/06_ESTIMATE/results/estimate_scores.csv
#   scripts/06_ESTIMATE/results/bulk_tpm_for_estimate.txt   (Step 7 control)
#
# OUTPUTS
#   scripts/12_Stemness/results/stemness_scores.csv
#   scripts/12_Stemness/results/clinical_cache.csv
#   scripts/12_Stemness/results/wilcoxon_stemness_NE_vs_PN.csv
#   scripts/12_Stemness/results/stemness_sensitivity_bulk_vs_purified.csv
#   scripts/12_Stemness/results/cox_stemness_SUPPLEMENTARY.csv  (if computed)
#   scripts/12_Stemness/figures/*.png
#
# REFERENCES
#   Suva et al. (2014) Cell; Tirosh et al. (2016) Nature; Wang et al. (2008)
#   PLoS ONE; Tateishi et al. (2016) Clin Cancer Res; Malta et al. (2018) Cell.
# ==============================================================================
set.seed(123)
# ==============================================================================
# DEPENDENCIES
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(TCGAbiolinks)
})
# dplyr verbs must win over any namespace masking (e.g. TCGAbiolinks).
select <- dplyr::select
filter <- dplyr::filter
source("scripts/utils/tfm_theme.R")  # tfm_colors, tfm_subtype_colors, theme_tfm
# ==============================================================================
# CONFIGURATION
# ==============================================================================
RESULTS_DIR  <- "scripts/12_Stemness/results"
FIGURES_DIR  <- "scripts/12_Stemness/figures"
CONSENSUS    <- "scripts/08_Consenso/results/consensus_proportions.csv"
CONSENSUS_PC <- "scripts/08_Consenso/results/consensus_proportions_purity_corrected.csv"
BP_EXP_PATH  <- "scripts/02_BayesPrism/results/bayesprism_expression_glioma.rds"
ESTIMATE     <- "scripts/06_ESTIMATE/results/estimate_scores.csv"
BULK_TPM     <- "scripts/06_ESTIMATE/results/bulk_tpm_for_estimate.txt"
CLIN_CACHE   <- file.path(RESULTS_DIR, "clinical_cache.csv")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
# Neutral colours for tumour grade (a distinct axis from NE/PN subtype; do NOT
# reuse subtype colours here, that would imply a G2=NE / G3=PN equivalence).
GRADE_COLORS <- c("G2" = "#9ECAE1", "G3" = "#3182BD")
# Single source of truth for the stemness signature (used by both the main
# score and the sensitivity control).
STEMNESS_SIGNATURE <- c(
  "SOX2", "POU3F2", "SALL2", "MYC", "MYCN",            # GSC core + self-renewal
  "SOX4", "SOX11", "CCND2", "CTNNB1", "USP22", "MSI1"  # Tirosh 2016 stem-like
)
# ==============================================================================
# HELPERS
# ==============================================================================
# Stemness score from an expression matrix (genes x samples): per-gene z-score
# across samples, then the per-sample mean of the signature z-scores.
compute_stemness <- function(expr_matrix, signature) {
  genes_ok <- intersect(signature, rownames(expr_matrix))
  message(sprintf("    Signature genes present: %d / %d",
                  length(genes_ok), length(signature)))
  if (length(genes_ok) < 3) stop("Too few signature genes present (<3).")
  sub <- expr_matrix[genes_ok, , drop = FALSE]
  if (max(sub, na.rm = TRUE) > 100) sub <- log2(sub + 1)  # linear -> log2
  z <- t(scale(t(sub)))   # scale standardises columns -> transpose for per-gene
  z[is.na(z)] <- 0        # zero-variance genes -> 0
  setNames(as.numeric(colMeans(z, na.rm = TRUE)), colnames(expr_matrix))
}
norm_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)
# ==============================================================================
# STEP 1: LOAD CONSENSUS PROPORTIONS AND ESTIMATE PURITY
# ==============================================================================
message("[1] Loading consensus proportions and ESTIMATE...")
read_consensus <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  colnames(df)[1] <- "sample_id"
  colnames(df) <- make.names(colnames(df))
  df %>% mutate(sample_id_12 = substr(sample_id, 1, 12))
}
consensus_df    <- read_consensus(CONSENSUS)
consensus_pc_df <- read_consensus(CONSENSUS_PC)
estimate_df <- read.csv(ESTIMATE, check.names = FALSE)
colnames(estimate_df)[1] <- "sample_id"
colnames(estimate_df) <- make.names(colnames(estimate_df))
estimate_df <- estimate_df %>%
  mutate(sample_id_12 = substr(sample_id, 1, 12)) %>%
  dplyr::select(sample_id_12, TumorPurity)
message(sprintf("  Consensus samples: %d", nrow(consensus_df)))
# ==============================================================================
# STEP 1b: CLINICAL DATA (TCGAquery_subtype, cached locally)
#   Column names from TCGAquery_subtype have shifted across package versions, so
#   the survival column is resolved defensively: if the expected name is absent,
#   OS_days is set to NA (the survival step is supplementary and will skip).
# ==============================================================================
message("[1b] Loading clinical data...")
if (file.exists(CLIN_CACHE)) {
  clin_df <- read.csv(CLIN_CACHE, stringsAsFactors = FALSE)
  message(sprintf("  Loaded clinical cache: %d samples", nrow(clin_df)))
} else {
  message("  Querying TCGAquery_subtype(tumor='lgg') [requires internet]...")
  clinical_raw <- TCGAquery_subtype(tumor = "lgg")
  
  # Resolve the survival-time column defensively (name varies by version).
  surv_col <- intersect(c("Survival..months.", "Survival", "OS.time",
                          "OS_MONTHS"),
                        colnames(clinical_raw))
  if (length(surv_col) >= 1) {
    os_days <- as.numeric(clinical_raw[[surv_col[1]]]) * 30.44
    message(sprintf("  Survival column used: '%s'", surv_col[1]))
  } else {
    os_days <- NA_real_
    message("  WARNING: no survival column found; OS_days = NA ",
            "(supplementary survival step will skip).")
  }
  
  vital_col <- intersect(c("Vital.status..1.dead.", "vital_status",
                           "Vital.Status"),
                         colnames(clinical_raw))
  vital <- if (length(vital_col) >= 1)
    suppressWarnings(as.numeric(clinical_raw[[vital_col[1]]])) else NA_real_
  
  clin_df <- data.frame(
    sample_id_12          = substr(toupper(clinical_raw$patient), 1, 12),
    Transcriptome.Subtype = clinical_raw$Transcriptome.Subtype,
    Grade                 = clinical_raw$Grade,
    OS_days               = os_days,
    Vital.status..1.dead. = vital,
    ABSOLUTE.purity       = clinical_raw$ABSOLUTE.purity,
    IDH.codel.subtype     = clinical_raw$IDH.codel.subtype,
    Age                   = clinical_raw$Age..years.at.diagnosis.,
    stringsAsFactors      = FALSE
  ) %>% distinct(sample_id_12, .keep_all = TRUE)
  write.csv(clin_df, CLIN_CACHE, row.names = FALSE)
  message(sprintf("  Retrieved and cached %d samples", nrow(clin_df)))
}
message("  Subtype counts:")
print(table(clin_df$Transcriptome.Subtype))
# ==============================================================================
# STEP 2: PROXY 1 (proliferating proportion); assemble base dataset
# ==============================================================================
message("\n[2] Assembling stemness dataset (proxy 1: proliferating)...")
stemness_df <- data.frame(
  sample_id_12            = consensus_df$sample_id_12,
  Proliferating           = consensus_df$Proliferating,
  Glioma                  = consensus_df$Glioma,
  Glioma_purity_corrected = consensus_pc_df$Glioma,
  stringsAsFactors        = FALSE
) %>%
  merge(estimate_df, by = "sample_id_12", all.x = TRUE) %>%
  merge(clin_df,     by = "sample_id_12", all.x = TRUE)
message(sprintf("  Dataset: %d samples", nrow(stemness_df)))
# ==============================================================================
# STEP 3: PROXY 2 -- curated stemness score on BayesPrism Z matrix
# ==============================================================================
message("\n[3] Computing stemness score on BayesPrism Z matrix...")
common_sig    <- character(0)
stem_purified <- NULL
if (file.exists(BP_EXP_PATH)) {
  glioma_exp <- readRDS(BP_EXP_PATH)                 # samples x genes
  Z_gxs      <- t(as.matrix(glioma_exp))             # genes x samples
  common_sig <- intersect(STEMNESS_SIGNATURE, rownames(Z_gxs))
  
  if (length(common_sig) >= 3) {
    stem_purified <- compute_stemness(Z_gxs, STEMNESS_SIGNATURE)
    glioma_stem_df <- data.frame(
      sample_id_12   = substr(names(stem_purified), 1, 12),
      stemness_score = as.numeric(stem_purified),
      stringsAsFactors = FALSE
    )
    stemness_df <- merge(stemness_df, glioma_stem_df, by = "sample_id_12", all.x = TRUE)
    message(sprintf("  Scores for %d samples | range [%.3f, %.3f]",
                    sum(!is.na(stemness_df$stemness_score)),
                    min(stemness_df$stemness_score, na.rm = TRUE),
                    max(stemness_df$stemness_score, na.rm = TRUE)))
  } else {
    message("  WARNING: < 3 signature genes in Z matrix; stemness_score = NA.")
    stemness_df$stemness_score <- NA_real_
  }
} else {
  message("  WARNING: Z matrix not found; stemness_score = NA.")
  stemness_df$stemness_score <- NA_real_
}
# ==============================================================================
# STEP 4: SAVE SCORES
# ==============================================================================
write.csv(stemness_df, file.path(RESULTS_DIR, "stemness_scores.csv"), row.names = FALSE)
message("\n[4] Saved stemness_scores.csv")
# ==============================================================================
# STEP 5: FIGURES (all PNG, consistent with the rest of the pipeline)
# ==============================================================================
message("\n[5] Generating figures...")
score_median <- median(stemness_df$stemness_score, na.rm = TRUE)
p1 <- ggplot(stemness_df %>% filter(!is.na(stemness_score)), aes(x = stemness_score)) +
  geom_histogram(bins = 25, fill = tfm_colors[["Glioma"]], color = "white", alpha = 0.85) +
  geom_vline(xintercept = score_median, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  labs(title = "Distribution of Glioma stemness score",
       subtitle = sprintf("Curated GSC + stem-like signature | n=%d | median=%.3f",
                          sum(!is.na(stemness_df$stemness_score)), score_median),
       x = "Stemness score (mean z-score)", y = "Count") +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(FIGURES_DIR, "01_stemness_distribution.png"),
       p1, width = 6, height = 4, dpi = 300, bg = "white")
grade_df <- stemness_df %>% filter(Grade %in% c("G2", "G3"), !is.na(stemness_score))
if (nrow(grade_df) > 0) {
  p2 <- ggplot(grade_df, aes(x = Grade, y = stemness_score, fill = Grade)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
    scale_fill_manual(values = GRADE_COLORS) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
    labs(title = "Glioma stemness score by tumour grade",
         subtitle = "Curated signature | Wilcoxon rank-sum test",
         x = "Grade", y = "Stemness score") +
    theme_tfm(show_x_text = TRUE) + theme(legend.position = "none")
  ggsave(file.path(FIGURES_DIR, "02_stemness_by_grade.png"),
         p2, width = 4, height = 4, dpi = 300, bg = "white")
}
subtype_df <- stemness_df %>% filter(Transcriptome.Subtype %in% c("NE", "PN"),
                                     !is.na(stemness_score))
p3 <- ggplot(subtype_df, aes(x = Transcriptome.Subtype, y = stemness_score,
                             fill = Transcriptome.Subtype)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  scale_fill_manual(values = tfm_subtype_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  labs(title = "Glioma stemness score by transcriptome subtype",
       subtitle = sprintf("NE (n=%d) vs PN (n=%d) | Wilcoxon rank-sum test",
                          sum(subtype_df$Transcriptome.Subtype == "NE"),
                          sum(subtype_df$Transcriptome.Subtype == "PN")),
       x = "Transcriptome subtype", y = "Stemness score") +
  theme_tfm(show_x_text = TRUE) + theme(legend.position = "none")
ggsave(file.path(FIGURES_DIR, "03_stemness_by_subtype.png"),
       p3, width = 4, height = 4, dpi = 300, bg = "white")
ggsave(file.path(FIGURES_DIR, "fig1_stemness_NE_vs_PN.png"),
       p3, width = 5, height = 5, dpi = 300, bg = "white")
p4 <- ggplot(stemness_df %>% filter(!is.na(Glioma), !is.na(stemness_score)),
             aes(x = Glioma, y = stemness_score)) +
  geom_point(alpha = 0.6, color = tfm_colors[["Glioma"]], size = 2) +
  geom_smooth(method = "lm", color = "grey30", se = TRUE, linewidth = 0.8, formula = y ~ x) +
  stat_cor(method = "spearman") +
  labs(title = "Stemness score vs Glioma proportion",
       subtitle = "Integration with deconvolution",
       x = "Glioma proportion (consensus)", y = "Stemness score") +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(FIGURES_DIR, "04_stemness_vs_glioma.png"),
       p4, width = 5, height = 4, dpi = 300, bg = "white")
p5 <- ggplot(stemness_df %>% filter(!is.na(TumorPurity), !is.na(stemness_score)),
             aes(x = TumorPurity, y = stemness_score)) +
  geom_point(alpha = 0.6, color = tfm_colors[["Endo"]], size = 2) +
  geom_smooth(method = "lm", color = "grey30", se = TRUE, linewidth = 0.8, formula = y ~ x) +
  stat_cor(method = "spearman") +
  labs(title = "Stemness score vs tumour purity",
       subtitle = "Sanity check: stemness should be independent of purity",
       x = "Tumour purity (ESTIMATE)", y = "Stemness score") +
  theme_tfm(show_x_text = TRUE)
ggsave(file.path(FIGURES_DIR, "05_stemness_vs_purity.png"),
       p5, width = 5, height = 4, dpi = 300, bg = "white")
message("  Figures 01-05 saved (PNG).")
# ==============================================================================
# STEP 6: PRIMARY TEST -- Wilcoxon stemness NE vs PN
# ==============================================================================
message("\n[6] Primary test: Wilcoxon stemness NE vs PN...")
ne_vec <- subtype_df$stemness_score[subtype_df$Transcriptome.Subtype == "NE"]
pn_vec <- subtype_df$stemness_score[subtype_df$Transcriptome.Subtype == "PN"]
wt <- wilcox.test(ne_vec, pn_vec, exact = FALSE)
sig <- if (wt$p.value < 0.001) "***" else if (wt$p.value < 0.01) "**" else
  if (wt$p.value < 0.05) "*" else "n.s."
message(sprintf("  NE median=%.3f (n=%d) | PN median=%.3f (n=%d)",
                median(ne_vec), length(ne_vec), median(pn_vec), length(pn_vec)))
message(sprintf("  delta=%.3f (PN-NE) | W=%.0f | p=%.4e %s",
                median(pn_vec) - median(ne_vec), wt$statistic, wt$p.value, sig))
write.csv(
  data.frame(comparison = "PN vs NE", n_NE = length(ne_vec), n_PN = length(pn_vec),
             median_NE = round(median(ne_vec), 3), median_PN = round(median(pn_vec), 3),
             delta_median = round(median(pn_vec) - median(ne_vec), 3),
             W_statistic = wt$statistic, p_value = wt$p.value, significance = sig,
             method = "Wilcoxon rank-sum (exact=FALSE)",
             signature = paste(STEMNESS_SIGNATURE, collapse = ";"),
             n_genes_used = length(common_sig)),
  file.path(RESULTS_DIR, "wilcoxon_stemness_NE_vs_PN.csv"), row.names = FALSE)
message("  Saved wilcoxon_stemness_NE_vs_PN.csv")
# ==============================================================================
# STEP 7: SENSITIVITY ANALYSIS -- purified vs bulk (circularity guard)
#   Reuses stem_purified (Step 3) and the same compute_stemness on the raw bulk.
# ==============================================================================
message("\n[7] Sensitivity: purified vs bulk...")
if (!is.null(stem_purified) && file.exists(BULK_TPM)) {
  bulk <- as.matrix(read.table(BULK_TPM, header = TRUE, sep = "\t",
                               row.names = 1, check.names = FALSE))
  message(sprintf("  Bulk TPM: %d genes x %d samples", nrow(bulk), ncol(bulk)))
  stem_bulk <- compute_stemness(bulk, STEMNESS_SIGNATURE)
  
  sens_df <- stemness_df %>%
    dplyr::select(sample_id_12, subtype = Transcriptome.Subtype,
                  stem_purified = stemness_score) %>%
    mutate(sid = norm_id(sample_id_12)) %>%
    inner_join(data.frame(sid = norm_id(names(stem_bulk)),
                          stem_bulk = as.numeric(stem_bulk)), by = "sid") %>%
    filter(subtype %in% c("NE", "PN"), !is.na(stem_purified))
  message(sprintf("  NE/PN samples with both scores: %d", nrow(sens_df)))
  
  w_pur  <- wilcox.test(stem_purified ~ subtype, data = sens_df)
  w_bulk <- wilcox.test(stem_bulk ~ subtype, data = sens_df)
  ct_sb  <- cor.test(sens_df$stem_purified, sens_df$stem_bulk, method = "spearman")
  message(sprintf("  Purified: W=%.0f, p=%.3e | Bulk: W=%.0f, p=%.3e | rho=%.3f",
                  w_pur$statistic, w_pur$p.value,
                  w_bulk$statistic, w_bulk$p.value, ct_sb$estimate))
  
  write.csv(sens_df, file.path(RESULTS_DIR, "stemness_sensitivity_bulk_vs_purified.csv"),
            row.names = FALSE)
  
  sens_long <- sens_df %>%
    pivot_longer(c(stem_purified, stem_bulk), names_to = "input", values_to = "score") %>%
    mutate(input = recode(input, stem_purified = "Purified (Z matrix)",
                          stem_bulk = "Bulk (control)"))
  p_sens <- ggplot(sens_long, aes(x = subtype, y = score, fill = subtype)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.85, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    facet_wrap(~ input) +
    scale_fill_manual(values = tfm_subtype_colors) +
    labs(title = "Stemness sensitivity: purified vs bulk",
         subtitle = sprintf("Purified p=%.1e | Bulk p=%.1e | concordance rho=%.2f",
                            w_pur$p.value, w_bulk$p.value, ct_sb$estimate),
         x = NULL, y = "Stemness score (mean z-score)") +
    theme_tfm(show_x_text = TRUE) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  ggsave(file.path(FIGURES_DIR, "stemness_sensitivity_bulk_vs_purified.png"),
         p_sens, width = 9, height = 5, dpi = 300, bg = "white")
  message("  Saved sensitivity outputs (PN>NE on both inputs => robust).")
} else {
  message("  [skip] purified score or bulk TPM unavailable.")
}
# ==============================================================================
# STEP 8: SURVIVAL (SUPPLEMENTARY, UNDERPOWERED -- not for main text)
# ==============================================================================
message("\n[8] Survival (SUPPLEMENTARY; underpowered, EPV<10)...")
surv_df <- stemness_df %>%
  filter(!is.na(OS_days), !is.na(Vital.status..1.dead.), !is.na(stemness_score),
         Grade %in% c("G2", "G3")) %>%
  mutate(stem_group = ifelse(stemness_score >= median(stemness_score, na.rm = TRUE),
                             "High stemness", "Low stemness"),
         event = as.integer(Vital.status..1.dead.))
message(sprintf("  n=%d | events=%d (EPV~%.1f for a 2-covariate Cox)",
                nrow(surv_df), sum(surv_df$event), sum(surv_df$event) / 2))
if (nrow(surv_df) >= 20 && sum(surv_df$event) >= 5) {
  message("  NOTE: supplementary only; insufficient power for inference.")
  fit_km <- survfit(Surv(OS_days, event) ~ stem_group, data = surv_df)
  p_km <- ggsurvplot(fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
                     risk.table.height = 0.28,
                     palette = c(tfm_subtype_colors[["PN"]], tfm_subtype_colors[["NE"]]),
                     title = "Overall survival by stemness (SUPPLEMENTARY, underpowered)",
                     xlab = "Time (days)", ylab = "Survival probability",
                     legend.title = "Stemness", ggtheme = theme_tfm(show_x_text = TRUE))
  png(file.path(FIGURES_DIR, "08_survival_stemness_SUPPLEMENTARY.png"),
      width = 7, height = 6, units = "in", res = 300, bg = "white")
  print(p_km); dev.off()
  cox_model <- coxph(Surv(OS_days, event) ~ stemness_score + Grade, data = surv_df)
  write.csv(broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE),
            file.path(RESULTS_DIR, "cox_stemness_SUPPLEMENTARY.csv"), row.names = FALSE)
  message("  Saved supplementary survival outputs (interpret with caution).")
} else {
  message("  Skipped: too few events even for a supplementary analysis.")
}
# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", strrep("=", 60))
message("  SCRIPT 12 (STEMNESS) COMPLETED")
message(strrep("=", 60))
message(sprintf("  Signature genes used : %d / %d",
                length(common_sig), length(STEMNESS_SIGNATURE)))
message(sprintf("  Primary result       : NE=%.3f vs PN=%.3f | p=%.4e %s",
                median(ne_vec), median(pn_vec), wt$p.value, sig))
message("  Sensitivity (Step 7) : purified vs bulk (circularity guard)")
message("  Survival   (Step 8)  : supplementary only (underpowered)")
message(sprintf("  Outputs in           : %s", RESULTS_DIR))
message(strrep("=", 60))