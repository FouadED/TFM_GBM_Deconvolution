###############################################################################
# DECONVOLUCIÓN MCP-counter - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

devtools::install_github("ebecht/MCPcounter", 
                         ref = "master", 
                         subdir = "Source")

library(MCPcounter)
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
# MCP-counter necesita TPM, no raw counts
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
# 3. DECONVOLUCIÓN CON MCP-counter
# No necesita referencia scRNA-seq propia
# Devuelve scores de abundancia, NO proporciones
# =============================================================================

mcp_scores <- MCPcounter.estimate(
  expression    = bulk_tpm,
  featuresType  = "HUGO_symbols"  # nombres de gen en formato SYMBOL
)

# mcp_scores es tipos_celulares x muestras
cat("\nTipos celulares estimados:\n")
print(rownames(mcp_scores))
cat("\nScores estimados:\n")
print(round(mcp_scores, 3))

# =============================================================================
# 4. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

mcp_df        <- as.data.frame(t(mcp_scores))
mcp_df$sample <- rownames(mcp_df)
mcp_df        <- mcp_df[, c("sample", setdiff(colnames(mcp_df), "sample"))]

write.csv(mcp_df, "results/mcpcounter_scores.csv", row.names = FALSE)

# =============================================================================
# 5. GRÁFICOS
# =============================================================================

# --- 5A. Barras apiladas por muestra ---
mcp_long <- melt(mcp_df,
                 id.vars       = "sample",
                 variable.name = "cell_type",
                 value.name    = "score")

p1 <- ggplot(mcp_long, aes(x = sample, y = score, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Abundancia celular estimada (MCP-counter) — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Score de abundancia",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/mcpcounter_barras.png", p1, width = 14, height = 6, dpi = 300)

# --- 5B. Boxplot por tipo celular ---
p2 <- ggplot(mcp_long, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(
    title = "Distribución de scores (MCP-counter) — Oligodendroglioma",
    x     = "Tipo celular",
    y     = "Score de abundancia"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

print(p2)
ggsave("figures/mcpcounter_boxplot.png", p2, width = 10, height = 6, dpi = 300)

cat("\n✓ Script MCP-counter completado\n")
