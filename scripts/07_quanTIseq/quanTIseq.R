###############################################################################
# DECONVOLUCIÓN quanTIseq - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

# devtools::install_github("omnideconv/immunedeconv")

library(immunedeconv)
library(readxl)
library(ggplot2)
library(reshape2)

# =============================================================================
# 1. CARGAR BULK RNA-SEQ
# =============================================================================

bulk_raw <- read_excel("../oligo_counts_protein.coding.xlsx")

gene_col <- colnames(bulk_raw)[1]
bulk_mat <- as.matrix(bulk_raw[, -1])
rownames(bulk_mat) <- bulk_raw[[gene_col]]
storage.mode(bulk_mat) <- "numeric"

cat("Dimensiones bulk:", nrow(bulk_mat), "genes x", ncol(bulk_mat), "muestras\n")

# =============================================================================
# 2. NORMALIZAR A TPM
# quanTIseq necesita TPM
# =============================================================================

counts_to_tpm <- function(counts) {
  rpk <- counts / 1000
  tpm <- t(t(rpk) / colSums(rpk) * 1e6)
  return(tpm)
}

bulk_tpm <- counts_to_tpm(bulk_mat)

cat("Suma por muestra (debe ser ~1e6):\n")
print(round(colSums(bulk_tpm)))

# =============================================================================
# 3. DECONVOLUCIÓN CON quanTIseq
# =============================================================================

result_quantiseq <- deconvolute(
  gene_expression = bulk_tpm,
  method          = "quantiseq",
  tumor           = TRUE,   # indica que son muestras tumorales
  scale_mrna      = TRUE
)

# Resultado es tipos_celulares x muestras
cat("\nTipos celulares estimados:\n")
print(result_quantiseq$cell_type)

# =============================================================================
# 4. PROCESAR RESULTADOS
# =============================================================================

# Convertir a matriz muestras x tipos
quant_df        <- as.data.frame(t(result_quantiseq[, -1]))
colnames(quant_df) <- result_quantiseq$cell_type
quant_df$sample <- rownames(quant_df)
quant_df        <- quant_df[, c("sample", setdiff(colnames(quant_df), "sample"))]

cat("\nProporciones estimadas quanTIseq:\n")
print(round(quant_df[, -1], 3))

# =============================================================================
# 5. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

write.csv(quant_df, "results/quantiseq_proportions.csv", row.names = FALSE)

# =============================================================================
# 6. GRÁFICOS
# =============================================================================

# --- 6A. Barras apiladas por muestra ---
quant_long <- melt(quant_df,
                   id.vars       = "sample",
                   variable.name = "cell_type",
                   value.name    = "proportion")

p1 <- ggplot(quant_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(
    title = "Composición inmune estimada (quanTIseq) — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Proporción",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/quantiseq_proporciones_barras.png", p1, width = 14, height = 6, dpi = 300)

# --- 6B. Boxplot por tipo celular ---
p2 <- ggplot(quant_long, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(
    title = "Distribución de proporciones (quanTIseq) — Oligodendroglioma",
    x     = "Tipo celular",
    y     = "Proporción"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

print(p2)
ggsave("figures/quantiseq_boxplot.png", p2, width = 10, height = 6, dpi = 300)

cat("\n✓ Script quanTIseq completado\n")
