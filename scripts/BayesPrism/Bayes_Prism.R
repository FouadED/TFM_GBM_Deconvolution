###############################################################################
# DECONVOLUCIÓN BayesPrism - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

#devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

library(BayesPrism)
library(Seurat)
library(readxl)
library(ggplot2)
library(reshape2)

# =============================================================================
# 1. CARGAR Y PREPARAR REFERENCIA SINGLE-CELL
# =============================================================================

seurat_obj <- readRDS("../data/seurat_oligodendroglioma_RAW.rds")

#  quitar BCells por n=3
seurat_filt <- subset(seurat_obj, subset = Assignment != "BCells")

# Extraer counts crudos — BayesPrism necesita counts, NO data normalizada
sc_matrix <- GetAssayData(
  object = seurat_filt,
  layer  = "counts",
  assay  = "RNA"
)

sc_matrix <- as.matrix(sc_matrix)  # genes x células


# Extraer etiquetas 
cell_type_vector <- as.character(seurat_filt@meta.data$Assignment)

cat("\nDistribución de tipos celulares:\n")
print(table(cell_type_vector))

# =============================================================================
# 2. CARGAR Y PREPARAR BULK RNA-SEQ
# =============================================================================

bulk_raw <- read_excel("oligo_counts_protein.coding.xlsx")

gene_col  <- colnames(bulk_raw)[1]
bulk_mat  <- as.matrix(bulk_raw[, -1])
rownames(bulk_mat) <- bulk_raw[[gene_col]]
storage.mode(bulk_mat) <- "numeric"

cat("\nDimensiones bulk:", nrow(bulk_mat), "genes x", ncol(bulk_mat), "muestras\n")


# =============================================================================
# 3. ALINEAR GENES
# =============================================================================

genes_comunes <- intersect(rownames(sc_matrix), rownames(bulk_mat))
cat("\nGenes en referencia:", nrow(sc_matrix))
cat("\nGenes en bulk:", nrow(bulk_mat))
cat("\nGenes compartidos:", length(genes_comunes), "\n")

if (length(genes_comunes) < 1000) {
  warning("Menos de 1000 genes comunes — revisa nomenclatura ENSEMBL vs SYMBOL")
}

sc_matrix_filt  <- sc_matrix[genes_comunes, ]
bulk_mat_filt   <- bulk_mat[genes_comunes, ]

# =============================================================================
# 4. SUBMUESTREAR REFERENCIA
# Máximo 300 células por tipo — evita problemas de memoria
# =============================================================================

set.seed(042)
tipos    <- unique(cell_type_vector)
idx_keep <- c()

for (tipo in tipos) {
  idx_tipo <- which(cell_type_vector == tipo)
  n_max    <- min(300, length(idx_tipo))
  idx_keep <- c(idx_keep, sample(idx_tipo, n_max))
}

sc_matrix_sub  <- sc_matrix_filt[, idx_keep]
cell_type_sub  <- cell_type_vector[idx_keep]

cat("\nCélulas tras submuestreo:\n")
print(table(cell_type_sub))

# =============================================================================
# 5. CONSTRUIR OBJETO PRISM
# BayesPrism quiere: muestras x genes (transpuesto respecto a lo habitual)
# tener cuidado con transformaciones log
# =============================================================================

myPrism <- new.prism(
  reference         = t(sc_matrix_sub),
  mixture           = t(bulk_mat_filt),
  input.type        = "count.matrix",
  cell.type.labels  = cell_type_sub,
  cell.state.labels = cell_type_sub,
  outlier.cut       = 0.01,
  outlier.fraction  = 0.1,
  key               = "Glioma"   # tipo celular tumoral de tu referencia
)

# =============================================================================
# 6. EJECUTAR DECONVOLUCIÓN
# AVISO: puede tardar 10-30 min
# =============================================================================

bp.res <- run.prism(
  prism   = myPrism,
  n.cores = 8   # ajustar cores  — comprobar con parallel::detectCores()
)

# =============================================================================
# 7. EXTRAER PROPORCIONES
# =============================================================================

theta_type <- get.fraction(
  bp            = bp.res,
  which.theta   = "final",
  state.or.type = "type"
)

# theta_type es muestras x tipos_celulares
cat("\nProporciones estimadas BayesPrism:\n")
print(round(theta_type, 3))

# =============================================================================
# 8. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

saveRDS(bp.res,     "results/bayesprism_full_result.rds")
write.csv(theta_type, "results/bayesprism_proportions.csv")

# =============================================================================
# 9. COMPARACIÓN DIRECTA MuSiC vs BayesPrism
# =============================================================================

music_prop <- read.csv("proporciones_estimadas_music.csv", row.names = 1)

# Alinear muestras — mismo orden
muestras_comunes <- intersect(rownames(music_prop), rownames(theta_type))
music_alin  <- music_prop[muestras_comunes, ]
bayes_alin  <- as.data.frame(theta_type[muestras_comunes, ])

cat("\nMuestras en MuSiC:", nrow(music_alin))
cat("\nMuestras en BayesPrism:", nrow(bayes_alin), "\n")

# Alinear tipos celulares comunes
tipos_comunes <- intersect(colnames(music_alin), colnames(bayes_alin))
cat("Tipos celulares comunes:", tipos_comunes, "\n")

# =============================================================================
# 10. GRÁFICOS
# =============================================================================

# --- 10A. Barras apiladas BayesPrism ---
theta_df       <- as.data.frame(theta_type)
theta_df$sample <- rownames(theta_df)

theta_long <- melt(theta_df,
                   id.vars      = "sample",
                   variable.name = "cell_type",
                   value.name   = "proportion")

p1 <- ggplot(theta_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(
    title = "Composición celular estimada (BayesPrism) — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Proporción",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/bayesprism_proporciones_barras.png", p1, width = 10, height = 6, dpi = 300)

# --- 10B. Correlación MuSiC vs BayesPrism por tipo celular ---
# Solo para tipos celulares presentes en ambos
for (tipo in tipos_comunes) {
  
  df_cor <- data.frame(
    MuSiC      = music_alin[, tipo],
    BayesPrism = bayes_alin[, tipo],
    sample     = muestras_comunes
  )
  
  cor_val <- round(cor(df_cor$MuSiC, df_cor$BayesPrism, 
                       method = "spearman"), 3)
  
  p <- ggplot(df_cor, aes(x = MuSiC, y = BayesPrism, label = sample)) +
    geom_point(size = 3, color = "#2C7BB6") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    geom_text(vjust = -0.5, size = 2.5) +
    labs(
      title    = paste0("MuSiC vs BayesPrism — ", tipo),
      subtitle = paste0("Spearman r = ", cor_val),
      x        = "MuSiC proporción",
      y        = "BayesPrism proporción"
    ) +
    theme_bw()
  
  ggsave(paste0("figures/correlacion_", tipo, ".png"),
         p, width = 6, height = 5, dpi = 300)
}

cat("\n✓ Script completado. Revisa carpetas results/ y figures/\n")