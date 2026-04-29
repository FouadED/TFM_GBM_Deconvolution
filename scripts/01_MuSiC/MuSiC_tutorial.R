###############################################################################
# DECONVOLUCIÓN MuSiC - DATOS DEL TUTORIAL 
###############################################################################

# Instalar librerías 
# install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("MuSiC")

# Cargar librerías
library(Biobase)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)
library(reshape2)


# =============================================================================
# 1. DESCARGAR DATOS DEL TUTORIAL (Kidney de Mouse)
# =============================================================================
# Estos son los datos que usa el tutorial oficial de MuSiC


# Descargar scRNA-seq del tutorial
url_sce <- "https://xuranw.github.io/MuSiC/data/Mousesub_sce.rds"
url_bulk <- "https://xuranw.github.io/MuSiC/data/Mousebulkeset.rds"

# Descargar scRNA
cat("Descargando scRNA-seq del tutorial...\n")
tryCatch({
  download.file(url_sce, "mouse_sce.rds", mode = "wb", quiet = TRUE)
  sce <- readRDS("mouse_sce.rds")
  cat("✓ scRNA-seq descargado\n\n")
}, error = function(e) {
  cat("Error en descarga de scRNA.\n\n")
})

# Descargar bulk
cat("Descargando bulk RNA-seq del tutorial...\n")
tryCatch({
  download.file(url_bulk, "mouse_bulk.rds", mode = "wb", quiet = TRUE)
  bulk_eset <- readRDS("mouse_bulk.rds")
  cat("✓ Bulk RNA-seq descargado\n\n")
}, error = function(e) {
  cat("Error en descarga de bulk.\n\n")
})

# =============================================================================
# 2. INFORMACIÓN DE LOS DATOS
# =============================================================================

cat("Información del dataset :\n")
cat("  scRNA-seq: ", nrow(sce), "genes ×", ncol(sce), "células\n")
cat("  Bulk RNA-seq: ", nrow(bulk_eset), "genes ×", ncol(bulk_eset), "muestras\n")
cat("  Tipos celulares: ", paste(unique(colData(sce)$cellType), collapse = ", "), "\n\n")

# =============================================================================
# 3. FILTRAR GENES COMUNES
# =============================================================================


genes_common <- intersect(rownames(sce), rownames(bulk_eset))
cat("Genes comunes:", length(genes_common), "\n\n")

sce_filt <- sce[genes_common, ]
bulk_filt <- bulk_eset[genes_common, ]

# =============================================================================
# 4. MuSiC
# =============================================================================

cat("Ejecutando music_prop() - EXACTAMENTE como en el tutorial...\n\n")

# Esta es la FUNCIÓN PRINCIPAL del tutorial
result <- music_prop(
  bulk.mtx = exprs(bulk_filt),
  sc.sce = sce_filt,
  clusters = 'cellType',        # Columna que define tipos celulares
  samples = 'sampleID',         # Columna que define muestras
  verbose = TRUE
)

# obtenemos el objeto results. Este contiene: Est.prop.weighted, Est.prop.allgene, Weight.gene, r.squared.full y Var.prop




# =============================================================================
# 5. EXTRAER PROPORCIONES estimadas
# =============================================================================

prop_est <- result$Est.prop.weighted

# Crear data.frame ordenado
prop_df <- as.data.frame(prop_est)
prop_df$sample <- rownames(prop_df)

# Reordenar columnas
cols_order <- c("sample", setdiff(colnames(prop_df), "sample"))
prop_df <- prop_df[, cols_order]

print(head(prop_df, 10))

cat("\n")


# =============================================================================
# 7. GRÁFICO DE LAS PROPORCIONES CELULARES
# =============================================================================

# Transformar a formato largo
prop_long <- melt(prop_df,  
                  id.vars = "sample",
                  variable.name = "cell_type",
                  value.name = "proportion")

# Gráfico de barras apiladas
ggplot(prop_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(title = "Composición celular por muestra",
       x = "Muestra", y = "Proporción", fill = "Tipo celular") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar
ggsave("proporciones_celulares.png", width = 10, height = 6, dpi = 300)

# =============================================================================
# 8. INFORMACIÓN TÉCNICA (DEL TUTORIAL)
# =============================================================================

cat("\n=== INFORMACIÓN TÉCNICA (del resultado de music_prop) ===\n\n")

cat("R² (bondad de ajuste):\n")
cat("  Media:", mean(result$r.squared.full, na.rm = TRUE), "\n")
cat("  Rango: [", min(result$r.squared.full, na.rm = TRUE), 
    ",", max(result$r.squared.full, na.rm = TRUE), "]\n\n")

# =============================================================================
# 9. GUARDAR RESULTADOS
# =============================================================================
write.csv(prop_df, "proporciones_tutorial_music.csv", row.names = FALSE)


