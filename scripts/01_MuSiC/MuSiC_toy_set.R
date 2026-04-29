###############################################################################
# DECONVOLUCIÓN MuSiC - TOY DATASET
###############################################################################

library(Biobase)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)
library(reshape2)

# =============================================================================
# 1. CARGAR DATOS
# =============================================================================

sce       <- readRDS("sc_ref_SingleCellExperiment.rds")
bulk_eset <- readRDS("bulk_data_ExpressionSet.rds")
true_prop <- read.table("true_proportions.txt", header = TRUE, sep = "\t", row.names = 1)

# =============================================================================
# 2. FILTRAR GENES COMUNES ENTRE SCE Y BULK
# =============================================================================

genes_common <- intersect(rownames(sce), rownames(bulk_eset))

sce_filt  <- sce[genes_common, ]
bulk_filt <- bulk_eset[genes_common, ]

# =============================================================================
# 3. DECONVOLUCIÓN CON MuSiC
# =============================================================================

result <- music_prop(
  bulk.mtx = exprs(bulk_filt),
  sc.sce   = sce_filt,
  clusters = "cellType",
  samples  = "sampleID",
  verbose  = TRUE
)

# =============================================================================
# 4. EXTRAER Y ORDENAR PROPORCIONES ESTIMADAS
# =============================================================================

prop_est       <- result$Est.prop.weighted
prop_df        <- as.data.frame(prop_est)
prop_df$sample <- rownames(prop_df)
prop_df        <- prop_df[, c("sample", setdiff(colnames(prop_df), "sample"))]

# =============================================================================
# 5. ALINEAR true_prop CON prop_est (mismas muestras y tipos celulares)
# =============================================================================

# Asignamos los nombres de muestra de prop_est a true_prop (Bulk1, Bulk2, ...)
# ya que ambos tienen el mismo orden pero distintos rownames (1,2,3 vs Bulk1,Bulk2,...)
rownames(true_prop) <- rownames(prop_est)

# Reordenamos columnas de true_prop para que coincidan con prop_est
celltypes_common <- intersect(colnames(prop_est), colnames(true_prop))
est_alin  <- prop_est[, celltypes_common]
true_alin <- as.matrix(true_prop)[, celltypes_common]

# =============================================================================
# 6. MÉTRICAS DE COMPARACIÓN (estimado vs. verdadero)
# =============================================================================

rmse <- function(a, b) sqrt(mean((a - b)^2))

metrics <- data.frame(
  cell_type   = celltypes_common,
  RMSE        = sapply(celltypes_common, function(ct) rmse(est_alin[, ct], true_alin[, ct])),
  Pearson_cor = sapply(celltypes_common, function(ct) cor(est_alin[, ct],  true_alin[, ct]))
)

print(metrics)
cat("RMSE global:", rmse(as.vector(est_alin), as.vector(true_alin)), "\n")

# =============================================================================
# 7. GRÁFICOS
# =============================================================================


# --- 7A. Barras apiladas: proporciones estimadas por muestra ---
prop_long <- melt(prop_df, id.vars = "sample",
                  variable.name = "cell_type", value.name = "proportion")

p1 <- ggplot(prop_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(title = "Composición celular estimada (MuSiC)",
       x = "Muestra", y = "Proporción", fill = "Tipo celular") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("proporciones_estimadas.png", p1, width = 10, height = 6, dpi = 300)



# --- 7B. Scatter: estimado vs. verdadero por tipo celular ---
comp_df <- data.frame(
  estimated = as.vector(est_alin),
  true      = as.vector(true_alin),
  cell_type = rep(celltypes_common, each = nrow(est_alin))
)

p2 <- ggplot(comp_df, aes(x = true, y = estimated, color = cell_type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ cell_type, scales = "free") +
  labs(title = "Proporciones estimadas vs. verdaderas",
       x = "Proporción verdadera", y = "Proporción estimada (MuSiC)") +
  theme_bw() +
  theme(legend.position = "none")

print(p2)
ggsave("estimado_vs_verdadero.png", p2, width = 12, height = 8, dpi = 300)



# --- 7C. Panel comparativo: verdadero vs estimado lado a lado ---
true_df        <- as.data.frame(true_alin)
true_df$sample <- rownames(true_df)
true_long      <- melt(true_df, id.vars = "sample",
                       variable.name = "cell_type", value.name = "proportion")
true_long$tipo <- "Verdadero"

est_df        <- as.data.frame(est_alin)
est_df$sample <- rownames(est_df)
est_long      <- melt(est_df, id.vars = "sample",
                      variable.name = "cell_type", value.name = "proportion")
est_long$tipo <- "Estimado (MuSiC)"

comp_long <- rbind(true_long, est_long)

p3 <- ggplot(comp_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ tipo) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(title = "Verdadero vs. Estimado — proporciones por muestra",
       x = "Muestra", y = "Proporción", fill = "Tipo celular") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

print(p3)
ggsave("verdadero_vs_estimado_barras.png", p3, width = 14, height = 6, dpi = 300)

# =============================================================================
# 8. GUARDAR RESULTADOS
# =============================================================================

write.csv(prop_df, "proporciones_estimadas.csv", row.names = FALSE)
write.csv(metrics, "metricas_comparacion.csv",    row.names = FALSE)
