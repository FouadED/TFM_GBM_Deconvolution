# TFM_GBM_Deconvolution
# Pipeline de Deconvolución Celular - TFM GBM

## Datos

Los archivos de datos no están incluidos en el repositorio por su tamaño.

### Archivos necesarios

| Archivo | Descripción | Fuente |
|---------|-------------|--------|
| seurat_oligodendroglioma_RAW_for_MuSiC.rds | Referencia scRNA-seq anotada | Proporcionado por el tutor |
| oligo_counts_protein.coding.xlsx | Bulk RNA-seq raw counts | Proporcionado por el tutor |

## Algoritmos implementados

### MuSiC2
**Tipo:** Con referencia scRNA-seq

**Input:**
- `seurat_oligodendroglioma_RAW_for_MuSiC.rds` — referencia scRNA-seq
- `oligo_counts_protein.coding.xlsx` — bulk RNA-seq (raw counts, genes x muestras)

**Output:**
- `proporciones_estimadas_music.csv` — proporciones celulares por muestra
- `proporciones_estimadas_reales.png` — barras apiladas
- `heatmap_proporciones.png` — heatmap de proporciones

---

### BayesPrism
**Tipo:** Con referencia scRNA-seq

**Input:**
- `seurat_oligodendroglioma_RAW_for_MuSiC.rds` — referencia scRNA-seq
- `oligo_counts_protein.coding.xlsx` — bulk RNA-seq (raw counts, genes x muestras)

**Output:**
- `results/bayesprism_full_result.rds` — objeto completo (necesario para pasos posteriores)
- `results/bayesprism_proportions.csv` — proporciones celulares por muestra
- `figures/bayesprism_proporciones_barras.png` — barras apiladas

---

### EPIC
**Tipo:** Sin referencia scRNA-seq (firmas internas TRef para tumor)

**Input:**
- `oligo_counts_protein.coding.xlsx` — bulk RNA-seq (raw counts, genes x muestras)

**Output:**
- `results/epic_proportions.csv` — proporciones celulares por muestra
- `figures/epic_proporciones_barras.png` — barras apiladas
- `figures/correlacion_EPIC_[tipo_celular].png` — correlación MuSiC vs EPIC por tipo celular

---

## Estructura del repositorio

