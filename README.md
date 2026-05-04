# Single-Cell Proteomics Benchmarking Pipeline
## Imputation Strategies and Batch Effect Correction Methods

[![R version](https://img.shields.io/badge/R-v4.5.1-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.21-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the complete bioinformatic pipeline developed for the Master's Thesis:

> **"Development and evaluation of a bioinformatic workflow for Single-Cell Proteomics: Detection probability modelling and batch effect correction"**
> Gabriel Herrera Prieto | VIU — Universidad Internacional de Valencia | 2025
> Supervisor: Dr. Henry Fabian Tobar Tosse

The pipeline benchmarks six downstream analysis workflows combining three missing value strategies (KNN, QRILC and the probabilistic **limpa** method) with multiple batch effect correction approaches (removeBatchEffect and ComBat) applied to the publicly available [Leduc et al. (2022)](https://www.nature.com/articles/s41592-022-01616-x) pSCoPE dataset.

---

## Dataset

- **Source:** `leduc2022_pSCoPE()` from the [scpdata](https://bioconductor.org/packages/scpdata/) Bioconductor package
- **Cells:** 2,412 TMT channels across 134 acquisition batches
- **Cell types:** Melanoma (WM989) and Monocyte (U-937)
- **Technology:** pSCoPE — DDA-TMT 18-plex with carrier proteome

---

## Pipeline Structure

The pipeline is organized in 29 sequential R scripts:

| Script | Description |
|--------|-------------|
| `00_setup_env_R.R` | Environment setup and package installation |
| `01_carga_y_exploracion.R` | Dataset download and exploratory analysis |
| `02_2_preparacion_benchmark.R` | Benchmark design — separation of raw data and author solutions |
| `02_validacion_reproducibilidad.R` | Reproducibility validation vs Leduc et al. (r = 0.88) |
| `03_filtrado_psm.R` | PSM-level filtering (PIF > 0.6, dart_qval < 0.01) |
| `04_filtrado_scr.R` | Sample-to-Carrier Ratio filtering (SCR < 0.05) |
| `05_normalizacion_y_agregacion.R` | Reference normalization and PSM-to-peptide aggregation |
| `06_union_de_datos.R` | joinAssays — global peptide matrix (22,159 peptides) |
| `07_filtrado_celulas.R` | Cell QC via Median CV (threshold = 0.42) |
| `08_normalizacion_peptidos.R` | Double centering normalization |
| `09_refinamiento_peptidos.R` | NA filtering (>99%) — 12,295 peptides × 1,569 cells |
| `10_agregacion_y_diagnostico.R` | Peptide-to-protein aggregation (2,847 proteins, 65.9% NAs) |
| `11_imputacion_knn.R` | KNN imputation (k=3) — MAR reference method |
| `11_1_imputacion_qrilc.R` | QRILC imputation — MNAR classical method |
| `12_correccion_lote.R` | removeBatchEffect — KNN baseline |
| `12_1_correccion_lote_QRILC.R` | removeBatchEffect — QRILC pipeline |
| `12_2_DE_combat.R` | ComBat batch correction + differential analysis (lcbatch & Set) |
| `12_3_DE_knn.R` | Differential analysis on KNN matrix |
| `13_visualizacion_pca.R` | PCA visualization — baseline SCoPE2 |
| `13_1_pca_y_DE_QRILC.R` | PCA and differential analysis — QRILC pipeline |
| `14_analisis_probabilistico_LIMPA.R` | LIMPA pipeline — base model (no batch correction) |
| `14_1_analisis_probabilistico_LIMPA_lcbatch.R` | LIMPA pipeline — lcbatch covariate (2 levels) |
| `14_2_analisis_probabilistico_LIMPA_set.R` | LIMPA pipeline — Set covariate (133 plates) |
| `15_benchmark_final.R` | Visual benchmarking — 6 PCA panels |
| `15_1_metricas_cuantitativas.R` | Quantitative benchmarking — ASW, ARI, N_sig, Overlap |
| `16_Volcano_Plot_LIMPA.R` | Volcano plot — LIMPA base |
| `16_1_Volcano_Plot_LIMPA_Set.R` | Volcano plot — LIMPA Set + comparative panel |
| `16_2_Volcano_Plot_QRILC.R` | Volcano plot — QRILC vs LIMPA Set panel |
| `16_3_Volcano_Plot_LIMPA_lcbatch_comparativa.R` | Volcano plot — 3-variant LIMPA comparison panel |


---

## Benchmarking Scenarios

| Scenario | Imputation | Batch Correction |
|----------|-----------|-----------------|
| 1 — Baseline | KNN (k=3) | removeBatchEffect (lcbatch + Channel) |
| 2 — ComBat lc | KNN (k=3) | ComBat (lcbatch, 2 levels) |
| 3 — LIMPA lc | DPC-Quant | lcbatch in linear model |
| 4 — QRILC | QRILC | removeBatchEffect |
| 5 — ComBat Set | KNN (k=3) | ComBat (Set, 133 plates) |
| 6 — LIMPA Set | DPC-Quant | Set in linear model (133 levels) |

---

## Requirements

### R packages

```r
# Bioconductor
BiocManager::install(c(
  "scp", "scpdata", "QFeatures", "SingleCellExperiment",
  "limma", "sva", "impute", "MsCoreUtils", "scater"
))

# CRAN
install.packages(c(
  "ggplot2", "patchwork", "cluster", "mclust",
  "EnhancedVolcano", "matrixStats", "dplyr", "scales"
))

# limpa (GitHub)
devtools::install_github("limpa-repo/limpa")
```

### System
- R v4.5.1
- Bioconductor 3.21

---

## Key Results

| Scenario | ASW | ARI | Sig. proteins (FDR < 5%) | Overlap top 50 (%) |
|----------|-----|-----|--------------------------|-------------------|
| 1 — Baseline | 0.067 | 0.797 | 1,007 | 38% |
| 2 — ComBat lc | 0.064 | 0.785 | 984 | 38% |
| 3 — LIMPA lc | 0.059 | 0.790 | 421 | 24% |
| 4 — QRILC | 0.033 | 0.000 | 762 | 38% |
| 5 — ComBat Set | 0.062 | 0.792 | 1,001 | 38% |
| 6 — LIMPA Set | 0.059 | 0.790 | 557 | 28% |

---

## Repository Structure

```
scp-imputation-batchcorrection-benchmark/
├── scripts/          # R scripts (numbered by execution order)
├── README.md
└── LICENSE
```

---

## Citation

If you use this pipeline, please cite:

> Herrera Prieto, G. (2025). *Development and evaluation of a bioinformatic workflow for Single-Cell Proteomics: Detection probability modelling and batch effect correction*. Master's Thesis, VIU — Universidad Internacional de Valencia.

And the original dataset:

> Leduc, A., Huffman, R.G., Cantlon, J., Khan, S., & Slavov, N. (2022). Exploring functional protein covariation across single cells using nPOP. *Genome Biology*, 23(1), 261.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
