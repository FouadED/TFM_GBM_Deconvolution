###############################################################################
# DECONVOLUCIÓN BayesPrism - DATOS REALES OLIGODENDROGLIOMA
# Optimizado con filtrado de genes variables
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

library(Matrix)

seurat_obj  <- readRDS("../seurat_oligodendroglioma_RAW.rds")
seurat_filt <- subset(seurat_obj, subset = Assignment != "BCells")

cat("Distribución de tipos celulares:\n")
print(table(seurat_filt@meta.data$Assignment))

# Extraer como sparse y reconstruir counts crudos
sc_sparse    <- GetAssayData(seurat_filt, layer = "counts", assay = "RNA")
lib_size     <- as.numeric(seurat_filt@meta.data$nCount_RNA)
scale_factors <- lib_size / 10000

sc_sparse@x  <- expm1(sc_sparse@x)
sc_sparse_raw <- sc_sparse %*% Diagonal(x = scale_factors)
sc_sparse_raw@x <- round(sc_sparse_raw@x)

# Convertir a matriz densa SOLO después de filtrar genes variables
# para no cargar 3.4GB innecesariamente
cat("\nMax counts:", max(sc_sparse_raw@x), "\n")
cat("Error relativo máximo:", max(abs(colSums(sc_sparse_raw) - lib_size) / lib_size), "\n")

cell_type_vector <- as.character(seurat_filt@meta.data$Assignment)

# =============================================================================
# 2. CARGAR BULK RNA-SEQ
# =============================================================================

bulk_raw <- read_excel("../oligo_counts_protein.coding.xlsx")

gene_col <- colnames(bulk_raw)[1]
bulk_mat <- as.matrix(bulk_raw[, -1])
rownames(bulk_mat) <- bulk_raw[[gene_col]]
storage.mode(bulk_mat) <- "numeric"

cat("Dimensiones bulk:", nrow(bulk_mat), "genes x", ncol(bulk_mat), "muestras\n")
cat("¿Bulk son counts enteros?", all(bulk_mat == floor(bulk_mat)), "\n")

# =============================================================================
# 3. FILTRAR TOP 3000 GENES MÁS VARIABLES
# Convertir a densa SOLO con los genes filtrados — evita el error de memoria
# =============================================================================

seurat_filt     <- FindVariableFeatures(seurat_filt, nfeatures = 3000)
genes_variables <- VariableFeatures(seurat_filt)

genes_comunes <- intersect(genes_variables, rownames(bulk_mat))
cat("Genes comunes tras filtrado:", length(genes_comunes), "\n")

# Ahora sí convertir a densa — solo 3000 genes, manejable
sc_matrix_filt <- as.matrix(sc_sparse_raw[genes_comunes, ])
bulk_mat_filt  <- bulk_mat[genes_comunes, ]

cat("Dimensiones referencia filtrada:", 
    nrow(sc_matrix_filt), "genes x", ncol(sc_matrix_filt), "células\n")

# =============================================================================
# 4. SUBMUESTREAR REFERENCIA
# Máximo 300 células por tipo — evita problemas de memoria
# =============================================================================

set.seed(123)
tipos    <- unique(cell_type_vector)
idx_keep <- c()

for (tipo in tipos) {
  idx_tipo <- which(cell_type_vector == tipo)
  n_max    <- min(300, length(idx_tipo))
  idx_keep <- c(idx_keep, sample(idx_tipo, n_max))
}

sc_matrix_sub <- sc_matrix_filt[, idx_keep]
cell_type_sub <- cell_type_vector[idx_keep]

cat("\nCélulas tras submuestreo:\n")
print(table(cell_type_sub))

# =============================================================================
# 5. CONSTRUIR OBJETO PRISM
# =============================================================================

myPrism <- new.prism(
  reference         = t(sc_matrix_sub),
  mixture           = t(bulk_mat_filt),
  input.type        = "count.matrix",
  cell.type.labels  = cell_type_sub,
  cell.state.labels = cell_type_sub,
  outlier.cut       = 0.01,
  outlier.fraction  = 0.1,
  key               = "Glioma"
)

# =============================================================================
# 6. EJECUTAR DECONVOLUCIÓN
# Con 3000 genes debería tardar 15-20 minutos
# =============================================================================

cat("\nIniciando deconvolución BayesPrism...\n")
cat("Hora de inicio:", format(Sys.time(), "%H:%M:%S"), "\n")

bp.res <- run.prism(
  prism   = myPrism,
  n.cores = parallel::detectCores() - 1
)

cat("Hora de fin:", format(Sys.time(), "%H:%M:%S"), "\n")

# =============================================================================
# 7. EXTRAER PROPORCIONES
# =============================================================================

theta_type <- get.fraction(
  bp            = bp.res,
  which.theta   = "final",
  state.or.type = "type"
)

cat("\nProporciones estimadas BayesPrism:\n")
print(round(theta_type, 3))

# =============================================================================
# 8. EXTRAER EXPRESIÓN PURIFICADA POR TIPO CELULAR
# Este es el output más valioso de BayesPrism para el índice de stemness
# Permite identificar genes específicos de células stem-like
# =============================================================================

# Expresión purificada de células de Glioma (stem-like)
expr_purificada <- get.exp(
  bp         = bp.res,
  state.or.type = "type",
  cell.name  = "Glioma"
)

cat("\nDimensiones expresión purificada Glioma:",
    nrow(expr_purificada), "muestras x", ncol(expr_purificada), "genes\n")

# =============================================================================
# 9. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

saveRDS(bp.res,         "results/bayesprism_full_result.rds")
saveRDS(expr_purificada,"results/bayesprism_expresion_glioma.rds")

write.csv(theta_type,    "results/bayesprism_proportions.csv")
write.csv(expr_purificada, "results/bayesprism_expresion_glioma.csv")

# =============================================================================
# 10. GRÁFICOS
# =============================================================================

theta_df        <- as.data.frame(theta_type)
theta_df$sample <- rownames(theta_df)

theta_long <- melt(theta_df,
                   id.vars       = "sample",
                   variable.name = "cell_type",
                   value.name    = "proportion")

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
ggsave("figures/bayesprism_proporciones_barras.png", p1, width = 14, height = 6, dpi = 300)

# --- Comparación MuSiC vs BayesPrism ---
music_prop <- read.csv("../MuSiC/results/proporciones_estimadas_music.csv",
                       row.names = 1)

tipos_comunes    <- intersect(colnames(theta_type), colnames(music_prop))
muestras_comunes <- intersect(rownames(theta_type), rownames(music_prop))

cat("\nTipos comunes MuSiC vs BayesPrism:", tipos_comunes, "\n")

for (tipo in tipos_comunes) {
  df_cor <- data.frame(
    MuSiC      = music_prop[muestras_comunes, tipo],
    BayesPrism = theta_type[muestras_comunes, tipo],
    sample     = muestras_comunes
  )
  
  cor_val <- round(cor(df_cor$MuSiC, df_cor$BayesPrism, method = "spearman"), 3)
  
  p <- ggplot(df_cor, aes(x = MuSiC, y = BayesPrism, label = sample)) +
    geom_point(size = 3, color = "#2C7BB6") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    geom_text(vjus