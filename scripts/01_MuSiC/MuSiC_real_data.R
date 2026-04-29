###############################################################################
# DECONVOLUCIÓN MuSiC - DATOS REALES OLIGODENDROGLIOMA
###############################################################################

library(Seurat)
library(Biobase)
library(SingleCellExperiment)
library(MuSiC)
library(readxl)
library(ggplot2)
library(reshape2)

# =============================================================================
# 1. CARGAR Y PREPARAR REFERENCIA SINGLE-CELL
# =============================================================================

seurat_obj <- readRDS("../seurat_oligodendroglioma_RAW.rds")

# Filtrar BCells
# Dada la cantidad de células (3) no hay suficiente señal estadística para construir un perfil de expresión fiable para ese tipo
seurat_filt <- subset(seurat_obj, subset = Assignment != "BCells")

# Convertir Seurat → SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_filt, assay = "RNA")

# Verificar que las columnas clave están en colData
stopifnot("Assignment" %in% colnames(colData(sce)))
stopifnot("Fragment"   %in% colnames(colData(sce)))

cat("Células por tipo celular tras filtrado:\n")
print(table(colData(sce)$Assignment))

cat("\nCélulas por fragmento:\n")
print(table(colData(sce)$Fragment))

# =============================================================================
# 2. CARGAR Y PREPARAR BULK RNA-SEQ
# =============================================================================

bulk_raw <- read_excel("../oligo_counts_protein.coding.xlsx")

gene_col <- colnames(bulk_raw)[1]
cat("\nPrimera columna (genes):", gene_col, "\n")
cat("Muestras bulk:", colnames(bulk_raw)[-1], "\n")

bulk_mat <- as.matrix(bulk_raw[, -1])
rownames(bulk_mat) <- bulk_raw[[gene_col]]
storage.mode(bulk_mat) <- "numeric"

bulk_eset <- ExpressionSet(assayData = bulk_mat)

cat("\nDimensiones bulk:", nrow(bulk_eset), "genes x", ncol(bulk_eset), "muestras\n")

# =============================================================================
# 3. FILTRAR GENES COMUNES
# =============================================================================

genes_common <- intersect(rownames(sce), rownames(bulk_eset))
cat("\nGenes comunes SCE-Bulk:", length(genes_common), "\n")

if (length(genes_common) < 1000) {
  warning("Menos de 1000 genes comunes — revisa nomenclatura (ENSEMBL vs SYMBOL)")
}

sce_filt  <- sce[genes_common, ]
bulk_filt <- bulk_eset[genes_common, ]

# =============================================================================
# 4. DECONVOLUCIÓN CON MuSiC
# =============================================================================

result <- music_prop(
  bulk.mtx = exprs(bulk_filt),
  sc.sce   = sce_filt,
  clusters = "Assignment",
  samples  = "Fragment",
  verbose  = TRUE
)

# =============================================================================
# 5. EXTRAER PROPORCIONES ESTIMADAS
# =============================================================================

prop_est       <- result$Est.prop.weighted
prop_df        <- as.data.frame(prop_est)
prop_df$sample <- rownames(prop_df)
prop_df        <- prop_df[, c("sample", setdiff(colnames(prop_df), "sample"))]

cat("\nProporciones estimadas:\n")
head(round(prop_est, 3))

# =============================================================================
# 6. GRÁFICOS
# =============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- 6A. Barras apiladas por muestra ---
prop_long <- melt(prop_df, id.vars = "sample",
                  variable.name = "cell_type", value.name = "proportion")

p1 <- ggplot(prop_long, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  labs(title = "Composición celular estimada (MuSiC) — Oligodendroglioma",
       x = "Muestra bulk", y = "Proporción", fill = "Tipo celular") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("figures/proporciones_estimadas_reales.png", p1, width = 10, height = 6, dpi = 300)



# =============================================================================
# 7. GUARDAR RESULTADOS
# =============================================================================

write.csv(prop_df, "results/proporciones_estimadas_music.csv", row.names = FALSE)

