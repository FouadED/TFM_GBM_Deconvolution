# Cell-type Deconvolution and Molecular Characterisation of IDH-mutant 1p/19q-codeleted Oligodendroglioma

Reproducible, open-source R pipeline for cell-type deconvolution and molecular characterisation of IDH-mutant, 1p/19q-codeleted oligodendroglioma from bulk RNA-seq data, oriented towards the generation of parameters for tissue-level tumour digital twins.

Analysis code accompanying the Master's thesis (MADOBIS, Universidad de Sevilla / Universidad Internacional de Andalucía).

## Overview

The pipeline estimates the cellular composition of a TCGA oligodendroglioma cohort (122 samples) from bulk RNA-seq, using a single-cell reference (patient LGG-04, GSE182109). It integrates two reference-based deconvolution methods (MuSiC and BayesPrism) into a consensus, validates the result against an independent tumour-purity estimate, and characterises the malignant compartment through stemness scoring, cell-type marker identification and protein-protein interaction network analysis.

## Data

- **Bulk RNA-seq:** TCGA oligodendroglioma (IDH-mutant, 1p/19q-codeleted), retrieved with TCGAbiolinks. Publicly available from the Genomic Data Commons (https://portal.gdc.cancer.gov/).
- **Single-cell reference:** patient LGG-04 from GSE182109 (Abdelfattah et al., 2022), available from the Gene Expression Omnibus.

Raw count matrices and large intermediate objects are excluded from the repository (see `.gitignore`) and can be regenerated from the scripts.

## Execution order

Scripts are organised in numbered directories under `scripts/`. The numbering does not strictly reflect execution order; the pipeline follows these dependencies:

1. **Data preparation** — bulk RNA-seq download and TPM conversion; single-cell reference preparation.
2. **Deconvolution and purity estimation** — once the data are prepared, the deconvolution methods (MuSiC, BayesPrism, EPIC), the reference-free methods (xCell 2.0, MCP-counter, quanTIseq) and ESTIMATE can be run independently.
3. **Consensus** — requires MuSiC and BayesPrism outputs.
4. **Validation** — requires the consensus and ESTIMATE purity.
5. **Benchmark** — requires the single-cell reference (simulation-based).
6. **Marker identification** — requires the BayesPrism purified expression.
7. **Stemness index** — requires the BayesPrism purified expression.
8. **Network analysis** — requires the identified markers.

(Replace this list with the exact script filenames and the precise run order if a master script is provided.)

## Requirements

Analyses were run under R 4.3.3 (CICA Hercules HPC cluster, Universidad de Sevilla) and R 4.5.2 (local). Key packages:

- **Deconvolution:** MuSiC 1.0.0, BayesPrism 2.2.2, EPIC 1.1.7, xCell2 1.2.3, MCP-counter 1.2.0, quanTIseq via immunedeconv 2.1.0, ESTIMATE 1.0.13
- **Single-cell:** Seurat 5.4.0, SingleCellExperiment 1.32.0
- **Networks / enrichment:** STRINGdb 2.22.0, igraph 2.2.1, clusterProfiler 4.18.4
- **Data handling:** dplyr 1.1.4, biomaRt

Computationally intensive steps (BayesPrism deconvolution, benchmark) were executed on an HPC cluster via SLURM.

## Citation

Eddaoudi Lakraichi F., del Campo Alcoba C., Nepomuceno Chamorro I., Nepomuceno Chamorro J.A., Hajji N. (2026). *Bulk transcriptomes resolve the cellular and molecular parameters of oligodendroglioma for tumour digital twins*. Master's thesis, MADOBIS, Universidad de Sevilla / Universidad Internacional de Andalucía.


