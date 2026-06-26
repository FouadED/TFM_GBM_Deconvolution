# ==============================================================================
# SCRIPT 10 - Benchmark LOFO
# Author: Fouad Eddaoudi Lakraichi
# TFM: Cell type deconvolution in oligodendroglioma
# ==============================================================================
# Design: Leave-One-Fragment-Out (LOFO)
#   3 iterations (LGG-04-1, LGG-04-2, LGG-04-3 as test)
#   N_SIM Dirichlet simulations per iteration
#   N_SIM_BAYES BayesPrism simulations (fewer due to computational cost)
# Methods: MuSiC, EPIC, BayesPrism, Consenso ROL1
# Metrics: RMSE, MAE, Pearson, Spearman, CCC
# Reference: Hu & Chikina 2024, Genome Biology
# ==============================================================================
set.seed(123)

# ==============================================================================
# CONFIG CICA
# ==============================================================================
N_SIM       <- 100
N_SIM_BAYES <- 100
N_CORES     <- 16

SEURAT_RDS <- "~/tfm_deconv/data/scrna/seurat_oligodendroglioma_RAW.rds"
BULK_XLSX  <- "~/tfm_deconv/data/bulk/oligo_counts_protein.coding.xlsx"
OUT_DIR    <- "~/tfm_deconv/results/benchmark"
FIG_DIR    <- "~/tfm_deconv/figures/benchmark"

MAX_CELLS_PER_TYPE <- 500
MIN_CELLS_TYPE     <- 20

# ==============================================================================
# LIBRARIES
# ==============================================================================
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(MuSiC)
  library(EPIC)
  library(BayesPrism)
  library(gtools)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(readxl)
  library(BiocParallel)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FUNCTION: Lin 1989 CCC
# ==============================================================================
ccc_lin <- function(x, y) {
  complete <- complete.cases(x, y)
  if (sum(complete) < 2) return(NA)
  x <- x[complete]; y <- y[complete]
  mx <- mean(x); my <- mean(y)
  sx <- var(x);  sy <- var(y)
  sxy <- cov(x, y)
  (2 * sxy) / (sx + sy + (mx - my)^2)
}

# ==============================================================================
# METRICS FUNCTION
# ==============================================================================
calc_m <- function(pred, true, fragmento, sim_i, metodo) {
  pred <- as.numeric(pred)
  true <- as.numeric(true)
  data.frame(
    fragmento = fragmento, sim = sim_i, metodo = metodo,
    RMSE     = sqrt(mean((pred - true)^2, na.rm = TRUE)),
    MAE      = mean(abs(pred - true), na.rm = TRUE),
    Pearson  = tryCatch(cor(pred, true, method = "pearson",
                            use = "complete.obs"), error = function(e) NA),
    Spearman = tryCatch(cor(pred, true, method = "spearman",
                            use = "complete.obs"), error = function(e) NA),
    CCC      = ccc_lin(pred, true)
  )
}

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================
cat("\nLoading Seurat object...\n")
seurat_obj <- readRDS(SEURAT_RDS)

asig <- as.character(seurat_obj@meta.data$Assignment)
sub  <- as.character(seurat_obj@meta.data$SubAssignment)
seurat_obj$cell_label <- ifelse(
  is.na(sub) | sub == "" | sub == "NA", asig, sub
)

frags <- unique(seurat_obj@meta.data$orig.ident)
cat("Fragments:", paste(frags, collapse = ", "), "\n")

tab <- table(seurat_obj$cell_label)
tipos_validos <- names(tab[tab >= MIN_CELLS_TYPE])
cat("Valid cell types (>=", MIN_CELLS_TYPE, "cells):", length(tipos_validos), "\n")
cat(paste(" -", tipos_validos), sep = "\n")

seurat_obj <- seurat_obj[, seurat_obj$cell_label %in% tipos_validos]

cat("\nReconstructing counts from log-normalized data...\n")
norm_matrix <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
lib_size    <- as.numeric(seurat_obj@meta.data$nCount_RNA)
sc_counts   <- round(sweep(expm1(norm_matrix), 2, lib_size / 10000, "*"))
all_genes   <- rownames(sc_counts)
cat("Total genes:", length(all_genes), "\n")

cat("Loading bulk...\n")
bulk_raw <- as.data.frame(read_excel(BULK_XLSX))
rownames(bulk_raw) <- bulk_raw[, 1]
bulk_raw    <- bulk_raw[, -1]
bulk_matrix <- as.matrix(bulk_raw)
cat("Bulk:", nrow(bulk_matrix), "genes x", ncol(bulk_matrix), "samples\n")

# ==============================================================================
# LOFO MAIN LOOP
# ==============================================================================
cat("\n=== START BENCHMARK LOFO ===\n")
resultados_todos <- list()

for (frag_test in frags) {
  cat("\n--- Iteration: test =", frag_test, "---\n")
  
  frags_train <- setdiff(frags, frag_test)
  
  cells_train <- colnames(seurat_obj)[
    seurat_obj@meta.data$orig.ident %in% frags_train
  ]
  cells_test <- colnames(seurat_obj)[
    seurat_obj@meta.data$orig.ident == frag_test
  ]
  
  seurat_train <- seurat_obj[, cells_train]
  seurat_test  <- seurat_obj[, cells_test]
  
  # Subsample training cells (cap per cell type)
  set.seed(123)
  cells_sub <- unlist(lapply(tipos_validos, function(ct) {
    idx <- which(seurat_train$cell_label == ct)
    if (length(idx) == 0) return(NULL)
    sample(idx, min(length(idx), MAX_CELLS_PER_TYPE))
  }))
  seurat_train_sub <- seurat_train[, cells_sub]
  sc_counts_train  <- sc_counts[, colnames(seurat_train_sub)]
  
  # --- MuSiC reference: SCE duplicated (single-subject workaround) ----------
  sce_A <- SingleCellExperiment(
    assays  = list(counts = sc_counts_train),
    colData = DataFrame(
      cellType = seurat_train_sub$cell_label,
      sampleID = "LGG-04-A"
    )
  )
  sce_B <- sce_A
  colData(sce_B)$sampleID <- "LGG-04-B"
  sce_music <- cbind(sce_A, sce_B)
  
  # --- EPIC reference: raw count means --------------------------------------
  # IMPORTANT: no CPM normalization, no explicit sigGenes.
  # EPIC handles scale internally via its mRNA renormalization step.
  # Adding CPM or explicit sigGenes caused internal dimension errors in
  # EPIC v1.1.7 -- reverting to original working setup.
  ref_epic <- lapply(tipos_validos, function(ct) {
    idx <- which(seurat_train_sub$cell_label == ct)
    if (length(idx) == 0) return(NULL)
    rowMeans(sc_counts_train[, idx, drop = FALSE])
  })
  names(ref_epic) <- tipos_validos
  ref_epic        <- Filter(Negate(is.null), ref_epic)
  ref_mat_epic    <- do.call(cbind, ref_epic)
  
  # No sigGenes: EPIC selects informative genes internally
  epic_ref <- list(refProfiles = ref_mat_epic)
  
  # --- BayesPrism reference -------------------------------------------------
  genes_comunes   <- all_genes
  bp_counts_train <- t(as.matrix(sc_counts_train))
  bp_labels_cell  <- seurat_train_sub$cell_label
  bp_labels_type  <- seurat_train_sub$cell_label
  
  cat("  Generating", N_SIM, "Dirichlet simulations...\n")
  res_iter <- data.frame()
  
  # ==========================================================================
  # MuSiC + EPIC SIMULATIONS
  # ==========================================================================
  for (sim_i in seq_len(N_SIM)) {
    alpha_true <- gtools::rdirichlet(1, rep(1, length(tipos_validos)))[1, ]
    names(alpha_true) <- tipos_validos
    
    # FIX: use colnames() not which() for correct cell indexing on sc_counts
    cells_by_type <- lapply(tipos_validos, function(ct)
      colnames(seurat_test)[seurat_test$cell_label == ct])
    names(cells_by_type) <- tipos_validos
    
    sampled_cells <- unlist(lapply(tipos_validos, function(ct) {
      n_ct     <- round(alpha_true[ct] * 500)
      cell_nms <- cells_by_type[[ct]]
      if (length(cell_nms) == 0 || n_ct == 0) return(NULL)
      sample(cell_nms, min(length(cell_nms), n_ct), replace = TRUE)
    }))
    
    pb_counts <- rowSums(sc_counts[, sampled_cells, drop = FALSE])
    
    # --- MuSiC ---------------------------------------------------------------
    pb_vec   <- pb_counts[genes_comunes]
    pb_music <- cbind(pb_vec, pb_vec)
    colnames(pb_music) <- c("pb", "pb_dup")
    
    music_res <- tryCatch({
      music_prop(
        bulk.mtx  = pb_music,
        sc.sce    = sce_music,
        clusters  = "cellType",
        samples   = "sampleID",
        select.ct = tipos_validos
      )$Est.prop.weighted[1, ]
    }, error = function(e) {
      cat("    MuSiC error sim", sim_i, ":", conditionMessage(e), "\n")
      rep(NA, length(tipos_validos))
    })
    if (!is.null(names(music_res)))
      music_res <- music_res[tipos_validos]
    
    # --- EPIC ----------------------------------------------------------------
    # Raw count pseudo-bulk with dimnames set in matrix() -- original setup
    pb_dense <- as.numeric(pb_counts[genes_comunes])
    pb_dense[is.na(pb_dense)] <- 0
    pb_bulk_epic <- matrix(pb_dense, ncol = 1,
                           dimnames = list(genes_comunes, "pb"))
    
    epic_res <- tryCatch({
      out   <- EPIC(bulk       = pb_bulk_epic,
                    reference  = epic_ref,
                    withOther  = TRUE,
                    scaleExprs = FALSE)
      cf    <- out$cellFractions
      props <- if (is.matrix(cf)) cf[1, ] else as.numeric(cf)
      names(props) <- if (is.matrix(cf)) colnames(cf) else names(cf)
      props <- props[setdiff(names(props), "otherCells")]
      s <- sum(props, na.rm = TRUE)
      if (s > 0) props / s else rep(NA, length(tipos_validos))
    }, error = function(e) {
      cat("    EPIC error sim", sim_i, ":", conditionMessage(e), "\n")
      rep(NA, length(tipos_validos))
    })
    if (!is.null(names(epic_res)))
      epic_res <- epic_res[tipos_validos]
    
    # --- Consensus MuSiC + EPIC ----------------------------------------------
    prop_mat     <- rbind(music_res, epic_res)
    consenso_res <- apply(prop_mat, 2, median, na.rm = TRUE)
    s <- sum(consenso_res, na.rm = TRUE)
    if (!is.na(s) && s > 0) consenso_res <- consenso_res / s
    
    # --- Metrics -------------------------------------------------------------
    true_props <- alpha_true[tipos_validos]
    res_iter   <- rbind(res_iter,
                        calc_m(music_res,    true_props, frag_test, sim_i, "MuSiC"),
                        calc_m(epic_res,     true_props, frag_test, sim_i, "EPIC"),
                        calc_m(consenso_res, true_props, frag_test, sim_i, "Consenso")
    )
  }
  
  # ==========================================================================
  # BayesPrism SIMULATIONS
  # ==========================================================================
  cat("  BayesPrism:", N_SIM_BAYES, "simulations...\n")
  
  for (sim_i in seq_len(N_SIM_BAYES)) {
    alpha_true <- gtools::rdirichlet(1, rep(1, length(tipos_validos)))[1, ]
    names(alpha_true) <- tipos_validos
    
    cells_by_type <- lapply(tipos_validos, function(ct)
      colnames(seurat_test)[seurat_test$cell_label == ct])
    names(cells_by_type) <- tipos_validos
    
    sampled_cells <- unlist(lapply(tipos_validos, function(ct) {
      n_ct     <- round(alpha_true[ct] * 500)
      cell_nms <- cells_by_type[[ct]]
      if (length(cell_nms) == 0 || n_ct == 0) return(NULL)
      sample(cell_nms, min(length(cell_nms), n_ct), replace = TRUE)
    }))
    
    pb_counts <- rowSums(sc_counts[, sampled_cells, drop = FALSE])
    pb_vec    <- as.numeric(pb_counts[genes_comunes])
    
    bp_bulk <- matrix(pb_vec, nrow = 1)
    rownames(bp_bulk) <- "pb"
    colnames(bp_bulk) <- genes_comunes
    
    bp_res <- tryCatch({
      myPrism <- new.prism(
        reference         = bp_counts_train,
        mixture           = bp_bulk,
        input.type        = "count.matrix",
        cell.type.labels  = bp_labels_type,
        cell.state.labels = bp_labels_cell,
        key               = "Glioma",
        outlier.cut       = 0.01,
        outlier.fraction  = 0.1
      )
      bp_out <- run.prism(prism = myPrism, n.cores = N_CORES)
      theta  <- get.fraction(bp            = bp_out,
                             which.theta   = "final",
                             state.or.type = "type")
      theta[1, tipos_validos]
    }, error = function(e) {
      cat("    BayesPrism error sim", sim_i, ":", conditionMessage(e), "\n")
      rep(NA, length(tipos_validos))
    })
    
    true_props <- alpha_true[tipos_validos]
    res_iter   <- rbind(res_iter,
                        calc_m(as.numeric(bp_res), true_props, frag_test, sim_i, "BayesPrism")
    )
  }
  
  resultados_todos[[frag_test]] <- res_iter
  cat("  Iteration", frag_test, "completed:", nrow(res_iter), "rows\n")
}

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================
cat("\nCombining results...\n")
resultados <- do.call(rbind, resultados_todos)

write.csv(resultados,
          file.path(OUT_DIR, "benchmark_resultados.csv"),
          row.names = FALSE)
cat("Results:", nrow(resultados), "rows saved\n")

resumen <- resultados %>%
  group_by(fragmento, metodo) %>%
  summarise(
    n           = n(),
    RMSE_med    = round(median(RMSE,    na.rm = TRUE), 4),
    CCC_med     = round(median(CCC,     na.rm = TRUE), 4),
    Pearson_med = round(median(Pearson, na.rm = TRUE), 4),
    .groups     = "drop"
  )
cat("\nSummary:\n"); print(resumen)

write.csv(resumen,
          file.path(OUT_DIR, "benchmark_resumen.csv"),
          row.names = FALSE)

# ==============================================================================
# FIGURES
# ==============================================================================
cat("\nGenerating figures...\n")
colores <- c("MuSiC"      = "#2E75B6",
             "EPIC"       = "#27AE60",
             "BayesPrism" = "#E74C3C",
             "Consenso"   = "#8E44AD")

p_rmse <- ggplot(resultados, aes(x = metodo, y = RMSE, fill = metodo)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  facet_wrap(~fragmento, labeller = label_both) +
  scale_fill_manual(values = colores) +
  labs(title    = "Benchmark LOFO - RMSE per method",
       subtitle = paste0("Leave-One-Fragment-Out | ", N_SIM,
                         " Dirichlet sim (", N_SIM_BAYES, " BayesPrism)"),
       x = NULL, y = "RMSE",
       caption = "Lower RMSE indicates higher precision") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "benchmark_rmse.png"),
       p_rmse, width = 12, height = 5, dpi = 300)

p_ccc <- ggplot(resultados, aes(x = metodo, y = CCC, fill = metodo)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "black") +
  facet_wrap(~fragmento, labeller = label_both) +
  scale_fill_manual(values = colores) +
  labs(title    = "Benchmark LOFO - CCC per method",
       subtitle = "Concordance Correlation Coefficient (Lin 1989)",
       x = NULL, y = "CCC",
       caption = "Higher CCC = better agreement | Line: CCC = 0.70") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "benchmark_ccc.png"),
       p_ccc, width = 12, height = 5, dpi = 300)

p_pear <- ggplot(resultados, aes(x = metodo, y = Pearson, fill = metodo)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "black") +
  facet_wrap(~fragmento, labeller = label_both) +
  scale_fill_manual(values = colores) +
  labs(title    = "Benchmark LOFO - Pearson per method",
       subtitle = "Pearson correlation between true and estimated proportions",
       x = NULL, y = "Pearson r",
       caption = "Higher Pearson indicates better linear correlation") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "benchmark_pearson.png"),
       p_pear, width = 12, height = 5, dpi = 300)

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n====================================================\n")
cat("  BENCHMARK LOFO COMPLETED\n")
cat("====================================================\n")
cat("MuSiC/EPIC/Consenso simulations:", N_SIM, "\n")
cat("BayesPrism simulations:         ", N_SIM_BAYES, "\n")
cat("Total rows:                     ", nrow(resultados), "\n\n")
print(resumen)
cat("\nScript 10 completed OK\n")