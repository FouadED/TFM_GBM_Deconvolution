###############################################################################
# ESTIMATE - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos = rforge, dependencies = TRUE)

library(estimate)
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
# =============================================================================

counts_to_tpm <- function(counts) {
  rpk <- counts / 1000
  tpm <- t(t(rpk) / colSums(rpk) * 1e6)
  return(tpm)
}

bulk_tpm <- counts_to_tpm(bulk_mat)

# =============================================================================
# 3. ESTIMATE NECESITA UN ARCHIVO DE TEXTO INTERMEDIO
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Guardar TPM como archivo de texto
write.table(
  bulk_tpm,
  file      = "results/bulk_tpm_for_estimate.txt",
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

# =============================================================================
# 4. FILTRAR GENES COMUNES CON FIRMAS DE ESTIMATE
# =============================================================================

filterCommonGenes(
  input.f  = "results/bulk_tpm_for_estimate.txt",
  output.f = "results/bulk_estimate_filtered.gct",
  id       = "GeneSymbol"
)

# =============================================================================
# 5. CALCULAR SCORES
# =============================================================================

estimateScore(
  input.ds   = "results/bulk_estimate_filtered.gct",
  output.ds  = "results/estimate_scores.gct",
  platform   = "illumina"  # ajusta si usas affymetrix
)

# =============================================================================
# 6. CARGAR Y PROCESAR RESULTADOS
# =============================================================================

scores_raw <- read.table(
  "results/estimate_scores.gct",
  skip      = 2,
  header    = TRUE,
  row.names = 1,
  sep       = "\t"
)

# Eliminar columna de descripción
scores_raw <- scores_raw[, -1]

# Transponer — muestras x scores
scores_df        <- as.data.frame(t(scores_raw))
scores_df$sample <- rownames(scores_df)
scores_df        <- scores_df[, c("sample", setdiff(colnames(scores_df), "sample"))]

cat("\nScores ESTIMATE:\n")
print(round(scores_df[, -1], 3))


# =============================================================================
# 7. GUARDAR RESULTADOS
# =============================================================================

write.csv(scores_df, "results/estimate_scores.csv", row.names = FALSE)

# =============================================================================
# 8. GRÁFICOS
# =============================================================================

# --- 8A. ESTIMATEScore por muestra ---
p1 <- ggplot(scores_df, aes(x = sample, y = ESTIMATEScore)) +
  geom_bar(stat = "identity", fill = "#2C7BB6", width = 0.7) +
  labs(
    title = "ESTIMATE Score — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "ESTIMATE Score"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/estimate_score.png", p1, width = 14, height = 6, dpi = 300)

# --- 8B. ImmuneScore y StromalScore por muestra ---
scores_long <- melt(
  scores_df[, c("sample", "ImmuneScore", "StromalScore")],
  id.vars       = "sample",
  variable.name = "score_type",
  value.name    = "score"
)

p2 <- ggplot(scores_long, aes(x = sample, y = score, fill = score_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "ImmuneScore y StromalScore (ESTIMATE) — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Score",
    fill  = "Tipo"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave("figures/estimate_immune_stroma_scores.png", p2, width = 14, height = 6, dpi = 300)
cat("\n✓ Script ESTIMATE completado\n")
