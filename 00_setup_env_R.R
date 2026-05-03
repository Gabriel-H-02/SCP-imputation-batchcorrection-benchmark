### ============================================================
### SCRIPT 00: CONFIGURACIÓN DEL ENTORNO - TFM SINGLE-CELL PROTEOMICS
### ============================================================
### Ejecutar UNA SOLA VEZ por entorno de trabajo.
### ============================================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

# --- Paquetes de Bioconductor ---
paquetes_bioc <- c(
  "scp",                  # Gestión de datos y flujos de trabajo SCP
  "scpdata",              # Repositorio de datasets de SCP curados
  "QFeatures",            # Infraestructura para jerarquías de datos
  "SingleCellExperiment", # Formato estándar para datos de célula única
  "sva",                  # Corrección de efectos de lote (ComBat)
  "limma",                # Modelos lineales y análisis diferencial
  "Seurat",               # Clustering y visualización UMAP
  "ggplot2",              # Gráficos de alta calidad
  "patchwork",            # Composición de figuras
  "scater",               # Control de calidad y visualización SC
  "ComplexHeatmap",       # Heatmaps profesionales
  "MsCoreUtils",          # Utilidades para imputación (KNN, QRILC)
  "dplyr",                # Manipulación de tablas y metadatos
  "impute",               # Imputación KNN (Bioconductor)
  "EnhancedVolcano",      # Volcano plots
  "tidyverse",            # Ecosistema de análisis de datos
  "cluster",              # Cálculo de métricas de clustering (ASW)
  "mclust",               # Métricas de clustering (ARI)
  "lobstr"                # Monitorización de memoria RAM
)


BiocManager::install(paquetes_bioc, update = FALSE, ask = FALSE)

# --- limpa: paquete central del proyecto ---

message("Instalando limpa...")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limpa")

# --- Crear estructura de carpetas del proyecto ---
dirs <- c("data", "results", "results/plots", "results/tables", 
          "results/metrics")
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Carpeta creada: ", d)
  }
}

# --- Verificación final ---
message("\n--- COMPROBANDO INSTALACIÓN ---")
todos <- c(paquetes_bioc, "limpa")
check <- sapply(todos, requireNamespace, quietly = TRUE)

if (all(check)) {
  message("Todo instalado correctamente. Entorno listo.")
} else {
  faltantes <- names(check)[!check]
  message("Atención: faltan los siguientes paquetes: ",
          paste(faltantes, collapse = ", "))
}
