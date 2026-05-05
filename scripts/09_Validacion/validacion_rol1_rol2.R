###############################################################################
# VALIDACIÓN AUTOMÁTICA ROL1 vs ROL2
# Compara el ranking de muestras entre métodos con y sin referencia
###############################################################################

library(ggplot2)
library(reshape2)

# =============================================================================
# 1. CARGAR RESULTADOS
# =============================================================================

# ROL 1 — consenso
consenso <- read.csv("scripts/08_Consenso/results/consenso_proporciones.csv",
                     row.names = 1)

# ROL 2 — métodos sin referencia
xcell    <- read.csv("scripts/04_xCell/results/xcell_scores.csv",
                     row.names = 1)
mcp      <- read.csv("scripts/05_MCP_Counter/results/mcpcounter_scores.csv",
                     row.names = 1)
quant    <- read.csv("scripts/07_quanTIseq/results/quantiseq_proportions.csv",
                     row.names = 1)

cat("Muestras consenso:", nrow(consenso), "\n")
cat("Muestras xCell:", nrow(xcell), "\n")
cat("Muestras MCP-counter:", nrow(mcp), "\n")
cat("Muestras quanTIseq:", nrow(quant), "\n")

# =============================================================================
# 2. MAPA DE EQUIVALENCIAS ROL1 → ROL2
# Cada tipo celular del ROL1 tiene equivalentes en ROL2
# =============================================================================

equivalencias <- list(
  Myeloid = list(
    xCell = c("Macrophages", "Monocytes", "DC"),
    mcp   = c("Monocytic.lineage"),
    quant = c("Macrophages.M1", "Macrophages.M2", "Monocytes")
  ),
  TCells = list(
    xCell = c("CD4..T.cells", "CD8..T.cells", "Tregs"),
    mcp   = c("T.cells"),
    quant = c("T.cells.CD4", "T.cells.CD8", "Tregs")
  ),
  Endo = list(
    xCell = c("Endothelial.cells"),
    mcp   = c("Endothelial.cells"),
    quant = NULL
  )
)

# =============================================================================
# 3. VALIDACIÓN AUTOMÁTICA
# Correlación de Spearman entre ROL1 y ROL2 para cada tipo celular
# =============================================================================

dir.create("scripts/09_Validacion/results", showWarnings = FALSE, recursive = TRUE)
dir.create("scripts/09_Validacion/figures", showWarnings = FALSE, recursive = TRUE)

# Muestras comunes entre todos
muestras_comunes <- Reduce(intersect, list(
  rownames(consenso),
  rownames(xcell),
  rownames(mcp),
  rownames(quant)
))
cat("\nMuestras comunes para validación:", length(muestras_comunes), "\n")

resultados_val <- data.frame()

for (tipo in names(equivalencias)) {
  
  prop_rol1 <- consenso[muestras_comunes, tipo]
  
  # xCell
  if (!is.null(equivalencias[[tipo]]$xCell)) {
    cols_xcell <- intersect(equivalencias[[tipo]]$xCell, colnames(xcell))
    if (length(cols_xcell) > 0) {
      score_xcell <- rowSums(xcell[muestras_comunes, cols_xcell, drop = FALSE])
      r_xcell <- cor(prop_rol1, score_xcell, method = "spearman")
      resultados_val <- rbind(resultados_val, data.frame(
        tipo     = tipo,
        metodo   = "xCell",
        spearman = round(r_xcell, 3),
        validado = ifelse(r_xcell > 0.5, "✓ OK", "⚠ REVISAR")
      ))
    }
  }
  
  # MCP-counter
  if (!is.null(equivalencias[[tipo]]$mcp)) {
    cols_mcp <- intersect(equivalencias[[tipo]]$mcp, colnames(mcp))
    if (length(cols_mcp) > 0) {
      score_mcp <- rowSums(mcp[muestras_comunes, cols_mcp, drop = FALSE])
      r_mcp <- cor(prop_rol1, score_mcp, method = "spearman")
      resultados_val <- rbind(resultados_val, data.frame(
        tipo     = tipo,
        metodo   = "MCP-counter",
        spearman = round(r_mcp, 3),
        validado = ifelse(r_mcp > 0.5, "✓ OK", "⚠ REVISAR")
      ))
    }
  }
  
  # quanTIseq
  if (!is.null(equivalencias[[tipo]]$quant)) {
    cols_quant <- intersect(equivalencias[[tipo]]$quant, colnames(quant))
    if (length(cols_quant) > 0) {
      score_quant <- rowSums(quant[muestras_comunes, cols_quant, drop = FALSE])
      r_quant <- cor(prop_rol1, score_quant, method = "spearman")
      resultados_val <- rbind(resultados_val, data.frame(
        tipo     = tipo,
        metodo   = "quanTIseq",
        spearman = round(r_quant, 3),
        validado = ifelse(r_quant > 0.5, "✓ OK", "⚠ REVISAR")
      ))
    }
  }
}

cat("\n═══════════════════════════════════════\n")
cat("  RESULTADOS DE VALIDACIÓN ROL1 vs ROL2\n")
cat("═══════════════════════════════════════\n")
print(resultados_val)

# =============================================================================
# 5. INTERPRETACIÓN AUTOMÁTICA
# =============================================================================

cat("\n═══════════════════════════════════════\n")
cat("  INTERPRETACIÓN\n")
cat("═══════════════════════════════════════\n")

for (tipo in unique(resultados_val$tipo)) {
  
  vals <- resultados_val$spearman[resultados_val$tipo == tipo]
  media_r <- mean(vals)
  
  if (media_r > 0.5) {
    estado <- "✓ VALIDADO — concordancia robusta entre ROL1 y ROL2"
  } else if (media_r > 0.3) {
    estado <- "~ PARCIAL — concordancia moderada, posible tipo celular raro"
  } else {
    estado <- "⚠ BAJO — tipo celular casi ausente o señal débil"
  }
  
  cat(tipo, "— media r:", round(media_r, 3), "—", estado, "\n")
}

# =============================================================================
# 4. GUARDAR Y VISUALIZAR
# =============================================================================

write.csv(resultados_val, 
          "scripts/09_Validacion/results/validacion_rol1_rol2.csv",
          row.names = FALSE)

# Gráfico de validación
p <- ggplot(resultados_val, aes(x = tipo, y = spearman, fill = validado)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~metodo) +
  scale_fill_manual(values = c("✓ OK" = "#27AE60", "⚠ REVISAR" = "#E74C3C")) +
  labs(
    title    = "Validación ROL1 vs ROL2 — Concordancia de rankings",
    subtitle = "Línea roja: umbral mínimo Spearman r = 0.5",
    x        = "Tipo celular",
    y        = "Spearman r",
    fill     = "Estado"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave("scripts/09_Validacion/figures/validacion_rol1_rol2.png",
       p, width = 12, height = 6, dpi = 300)

cat("\n✓ Validación completada\n")
