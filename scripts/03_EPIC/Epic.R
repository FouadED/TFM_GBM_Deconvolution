###############################################################################
# 03_EPIC — Deconvolución celular con EPIC
# Input:  seurat_processed.rds + bulk RNA-seq
# Output: proporciones celulares por muestra
# Genérico para cualquier Seurat procesado con cell_label
###############################################################################

library(EPIC)
library(Seurat)
library(Matrix)
library(readxl)
library(ggplot2)
library(reshape2)

# =============================================================================
# CONFIGURACIÓN — Solo toca este bloque
# =============================================================================

config <- list(
  seurat_path  = "data/seurat_processed.rds",
  bulk_path    = "data/oligo_counts_protein.coding.xlsx",
  cell_type_col = "cell_label",
  music_csv    = "scripts/01_MuSiC/results/proporciones_estimadas_music.csv",
  output_dir   = "scripts/03_EPIC/results",
  figures_dir  = "scripts/03_EPIC/figures"
)

# =============================================================================
# 1. CARGAR BULK RNA-SEQ
# =============================================================================

cat("Cargando bulk RNA-seq...\n")
bulk_raw <- read_excel(config$bulk_path)

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
cat("Suma por muestra (debe ser ~1e6):\n")
print(round(colSums(bulk_tpm)[1:5]))

# =============================================================================
# 3. CONSTRUIR REFERENCIA DESDE SEURAT PROCESADO
# Usa cell_label — incluye subtipos
# =============================================================================

cat("\nCargando Seurat procesado...\n")
seurat_filt <- readRDS(config$seurat_path)

cat("Tipos celulares en referencia:\n")
print(table(seurat_filt@meta.data[[config$cell_type_col]]))

# Extraer counts como sparse y reconstruir counts crudos
lib_size  <- as.numeric(seurat_filt@meta.data$nCount_RNA)
sc_sparse <- GetAssayData(seurat_filt, layer = "counts", assay = "RNA")
sc_sparse@x <- expm1(sc_sparse@x)
sc_sparse_raw <- sc_sparse %*% Matrix::Diagonal(x = lib_size / 10000)
sc_sparse_raw@x <- round(sc_sparse_raw@x)

# Calcular perfil medio por tipo celular
cell_types   <- seurat_filt@meta.data[[config$cell_type_col]]
tipos        <- unique(cell_types)

ref_profiles <- sapply(tipos, function(tipo) {
  idx <- which(cell_types == tipo)
  Matrix::rowMeans(sc_sparse_raw[, idx])
})

cat("\nDimensiones referencia:", nrow(ref_profiles), "genes x", ncol(ref_profiles), "tipos\n")

# =============================================================================
# 4. ALINEAR GENES
# =============================================================================

genes_comunes <- intersect(rownames(ref_profiles), rownames(bulk_tpm))
cat("Genes compartidos referencia-bulk:", length(genes_comunes), "\n")

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

cat("\nEjecutando EPIC...\n")
cat("Hora de inicio:", format(Sys.time(), "%H:%M:%S"), "\n")

result_epic <- EPIC(
  bulk           = bulk_tpm_filt,
  reference      = mi_referencia,
  withOtherCells = TRUE
)

cat("Hora de fin:", format(Sys.time(), "%H:%M:%S"), "\n")

# =============================================================================
# 6. EXTRAER PROPORCIONES
# =============================================================================

prop_epic        <- as.data.frame(result_epic$cellFractions)
prop_epic$sample <- rownames(prop_epic)
prop_epic        <- prop_epic[, c("sample", setdiff(colnames(prop_epic), "sample"))]

cat("\nProporciones estimadas EPIC:\n")
print(round(head(result_epic$cellFractions), 3))

# =============================================================================
# 7. GUARDAR RESULTADOS
# =============================================================================

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$figures_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(prop_epic,
          file.path(config$output_dir, "epic_proportions.csv"),
          row.names = FALSE)

# =============================================================================
# 8. GRÁFICOS
# =============================================================================

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

print(p1)
ggsave(file.path(config$figures_dir, "epic_proporciones_barras.png"),
       p1, width = 14, height = 6, dpi = 300)

# --- Comparación con MuSiC ---
if (file.exists(config$music_csv)) {
  
  music_prop <- read.csv(config$music_csv, row.names = 1)
  
  tipos_comunes    <- intersect(colnames(result_epic$cellFractions), colnames(music_prop))
  muestras_comunes <- intersect(rownames(result_epic$cellFractions), rownames(music_prop))
  
  cat("\nTipos comunes MuSiC vs EPIC:", tipos_comunes, "\n")
  
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
    
    ggsave(file.path(config$figures_dir, paste0("correlacion_EPIC_", tipo, ".png")),
           p, width = 6, height = 5, dpi = 300)
  }
}

cat("\n✓ Script EPIC completado\n")

