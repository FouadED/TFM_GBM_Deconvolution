###############################################################################
# DECONVOLUCIÓN xCell - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

# devtools::install_github("dviraran/xCell")

library(xCell)
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
# xCell necesita TPM o FPKM, no raw counts
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
# 3. DECONVOLUCIÓN CON xCell
# No necesita referencia scRNA-seq propia
# Devuelve enrichment scores, NO proporciones
# =============================================================================

xcell_scores <- xCellAnalysis(bulk_tpm)

# xcell_scores es tipos_celulares x muestras
cat("\nDimensiones resultado:", nrow(xcell_scores), "tipos x", ncol(xcell_scores), "muestras\n")
cat("\nTipos celulares estimados:\n")
print(rownames(xcell_scores))

# =============================================================================
# 4. FILTRAR TIPOS CELULARES RELEVANTES PARA GBM
# =============================================================================

tipos_relevantes <- c(
  "CD4+ T-cells", "CD8+ T-cells", "Tregs",
  "Macrophages", "Monocytes", "DC",
  "NK cells", "Astrocytes", "Neurons",
  "Oligodendrocytes", "Endothelial cells"
)

# Filtrar solo los que están en el resultado
tipos_disponibles <- intersect(tipos_relevantes, rownames(xcell_scores))
cat("\nTipos relevantes disponibles:", tipos_disponibles, "\n")

xcell_filt <- xcell_scores[tipos_disponibles, ]

# =============================================================================
# 5. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

xcell_df        <- as.data.frame(t(xcell_filt))
xcell_df$sample <- rownames(xcell_df)
xcell_df        <- xcell_df[, c("sample", setdiff(colnames(xcell_df), "sample"))]

write.csv(xcell_df, "results/xcell_scores.csv", row.names = FALSE)

# =============================================================================
# 6. GRÁFICOS
# =============================================================================

# --- 6. Barras por tipo celular ---
p2 <- ggplot(xcell_melt, aes(x = sample, y = score, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Enrichment scores xCell — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Score",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave("figures/xcell_barras_scores.png", p2, width = 14, height = 6, dpi = 300)

cat("\n✓ Script xCell completado\n")
