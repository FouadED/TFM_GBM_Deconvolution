###############################################################################
# DECONVOLUCIÓN EPIC - DATOS REALES OLIGODENDROGLIOMA
# Con referencia scRNA-seq propia
###############################################################################

#devtools::install_github("GfellerLab/EPIC")

library(EPIC)
library(Seurat)
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
# EPIC necesita TPM, no raw counts
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
# 3. CONSTRUIR REFERENCIA SEURAT
# =============================================================================

seurat_obj  <- readRDS("../seurat_oligodendroglioma_RAW.rds")
seurat_filt <- subset(seurat_obj, subset = Assignment != "BCells")

cat("\nDistribución de tipos celulares:\n")
print(table(seurat_filt@meta.data$Assignment))

# Calcular perfil medio de expresión por tipo celular
sc_counts  <- GetAssayData(seurat_filt, layer = "counts", assay = "RNA")
cell_types <- seurat_filt@meta.data$Assignment
tipos      <- unique(cell_types)

ref_profiles <- sapply(tipos, function(tipo) {
  idx <- which(cell_types == tipo)
  rowMeans(sc_counts[, idx])
})

# ref_profiles: genes x tipos celulares
cat("\nDimensiones referencia:", nrow(ref_profiles), "genes x", ncol(ref_profiles), "tipos\n")
cat("Tipos celulares en referencia:", colnames(ref_profiles), "\n")

# =============================================================================
# 4. ALINEAR GENES ENTRE REFERENCIA Y BULK
# =============================================================================

genes_comunes <- intersect(rownames(ref_profiles), rownames(bulk_tpm))
cat("\nGenes compartidos referencia-bulk:", length(genes_comunes), "\n")

if (length(genes_comunes) < 1000) {
  warning("Menos de 1000 genes comunes — revisa nomenclatura ENSEMBL vs SYMBOL")
}

ref_profiles_filt <- ref_profiles[genes_comunes, ]
bulk_tpm_filt     <- bulk_tpm[genes_comunes, ]

# =============================================================================
# 5. DECONVOLUCIÓN CON EPIC
# =============================================================================

mi_referencia <- list(
  refProfiles = ref_profiles_filt,
  sigGenes    = rownames(ref_profiles_filt)
)

result_epic <- EPIC(
  bulk           = bulk_tpm_filt,
  reference      = mi_referencia,
  withOtherCells = TRUE
)

# =============================================================================
# 6. EXTRAER PROPORCIONES
# =============================================================================

prop_epic        <- as.data.frame(result_epic$cellFractions)
prop_epic$sample <- rownames(prop_epic)
prop_epic        <- prop_epic[, c("sample", setdiff(colnames(prop_epic), "sample"))]

cat("\nProporciones estimadas EPIC:\n")
print(round(result_epic$cellFractions, 3))

# =============================================================================
# 7. GUARDAR RESULTADOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

write.csv(prop_epic, "results/epic_proportions.csv", row.names = FALSE)

# =============================================================================
# 8. GRÁFICOS
# =============================================================================

# --- 8A. Barras apiladas ---
prop_long <- melt(prop_epic,
                  id.vars       = "sample",
                  variable.name = "cell_type",
                  value.name    = "proportion")

p1 <- ggplot(prop_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(
    title = "Composición celular estimada (EPIC) — Oligodendroglioma",
    x     = "Muestra bulk",
    y     = "Proporción",
    fill  = "Tipo celular"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/epic_proporciones_barras.png", p1, width = 10, height = 6, dpi = 300)

# --- 8B. Comparación con MuSiC ---
music_prop <- read.csv("../MuSiC/results/proporciones_estimadas_music.csv", row.names = 1)

tipos_comunes <- intersect(
  colnames(result_epic$cellFractions),
  colnames(music_prop)
)
cat("\nTipos celulares comunes MuSiC vs EPIC:", tipos_comunes, "\n")

muestras_comunes <- intersect(
  rownames(result_epic$cellFractions),
  rownames(music_prop)
)

for (tipo in tipos_comunes) {
  
  df_cor <- data.frame(
    MuSiC  = music_prop[muestras_comunes, tipo],
    EPIC   = result_epic$cellFractions[muestras_comunes, tipo],
    sample = muestras_comunes
  )
  
  cor_val <- round(cor(df_cor$MuSiC, df_cor$EPIC, method = "spearman"), 3)
  
  p <- ggplot(df_cor, aes(x = MuSiC, y = EPIC, label = sample)) +
    geom_point(size = 3, color = "#D7191C") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    geom_text(vjust = -0.5, size = 2.5) +
    labs(
      title    = paste0("MuSiC vs EPIC — ", tipo),
      subtitle = paste0("Spearman r = ", cor_val),
      x        = "MuSiC proporción",
      y        = "EPIC proporción"
    ) +
    theme_bw()
  
  ggsave(paste0("figures/correlacion_EPIC_", tipo, ".png"),
         p, width = 6, height = 5, dpi = 300)
}

cat("\n✓ Script EPIC completado\n")

