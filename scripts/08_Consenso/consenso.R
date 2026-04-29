library(ggplot2)
library(reshape2)

# =============================================================================
# CONSENSO — Mediana de MuSiC2 + EPIC + BayesPrism
# =============================================================================

# Filtrar y alinear
music_filt <- as.matrix(music[muestras_comunes, tipos_comunes])
epic_filt  <- as.matrix(epic[muestras_comunes,  tipos_comunes])
bayes_filt <- as.matrix(bayes[muestras_comunes, tipos_comunes])

# Calcular mediana por muestra y tipo celular
consenso <- matrix(
  NA,
  nrow     = length(muestras_comunes),
  ncol     = length(tipos_comunes),
  dimnames = list(muestras_comunes, tipos_comunes)
)

for (muestra in muestras_comunes) {
  for (tipo in tipos_comunes) {
    consenso[muestra, tipo] <- median(c(
      music_filt[muestra, tipo],
      epic_filt[muestra,  tipo],
      bayes_filt[muestra, tipo]
    ))
  }
}

# Renormalizar a 1
consenso <- t(apply(consenso, 1, function(x) x / sum(x)))

cat("¿Filas suman 1?", all(abs(rowSums(consenso) - 1) < 1e-10), "\n")
cat("\nPrimeras proporciones consenso:\n")
print(round(head(consenso), 3))



# Correlación de Spearman entre algoritmos para Glioma
cor_music_bayes <- cor(music_filt[, "Glioma"], 
                       bayes_filt[, "Glioma"], 
                       method = "spearman")
cor_music_epic  <- cor(music_filt[, "Glioma"], 
                       epic_filt[, "Glioma"], 
                       method = "spearman")
cor_epic_bayes  <- cor(epic_filt[, "Glioma"], 
                       bayes_filt[, "Glioma"], 
                       method = "spearman")

cat("Spearman MuSiC vs BayesPrism:", round(cor_music_bayes, 3), "\n")
cat("Spearman MuSiC vs EPIC:", round(cor_music_epic, 3), "\n")
cat("Spearman EPIC vs BayesPrism:", round(cor_epic_bayes, 3), "\n")

# Guardar
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

write.csv(consenso, "results/consenso_proporciones.csv")
cat("\nConsenso guardado\n")

# =============================================================================
# GRÁFICOS
# =============================================================================

consenso_df        <- as.data.frame(consenso)
consenso_df$sample <- rownames(consenso_df)

consenso_long <- melt(consenso_df,
                      id.vars       = "sample",
                      variable.name = "cell_type",
                      value.name    = "proportion")

p1 <- ggplot(consenso_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(
    title = "Composición celular — CONSENSO (MuSiC2 + EPIC + BayesPrism)",
    x     = "Muestra bulk",
    y     = "Proporción",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

print(p1)
ggsave("figures/consenso_proporciones.png", p1, width = 16, height = 6, dpi = 300)

cat("\n Script de consenso completado\n")

