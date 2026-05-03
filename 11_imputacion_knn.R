### ============================================================
### SCRIPT 11: IMPUTACIÓN KNN (MÉTODO DE REFERENCIA MAR)
### ============================================================
### Implementa la imputación por K vecinos más cercanos (KNN,
### k=3) como método de control clásico MAR. KNN estima los
### valores faltantes a partir de las k células más similares
### en el espacio proteómico.
### Posición en el benchmarking:
###   KNN (MAR): ignora el mecanismo de ausencia, imputa a
###              partir de similitud entre células.
###   QRILC (MNAR): asume ausencia por baja abundancia, imputa
###                 desde la cola izquierda de la distribución.
###   LIMPA (probabilístico): modela la DPC y propaga la
###                           incertidumbre al análisis diferencial.
###
### La imputación se realiza en el ESPACIO DE CÉLULAS
### (matriz transpuesta), buscando vecinos entre perfiles
### celulares en lugar de entre proteínas, tal como especifica
### el protocolo SCoPE2 original.
###
### INPUT:  data/leduc_proteins_final.rds
### OUTPUT: data/leduc_imputed_knn.rds
###         results/metrics/tiempo_imputacion_knn.rds
### ============================================================

library(scp)
library(SingleCellExperiment)
library(impute)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_proteins_final.rds")
message("Proteínas disponibles: ", nrow(leduc[["proteins_norm2"]]))
message("Células disponibles:   ", ncol(leduc[["proteins_norm2"]]))

mat_original  <- assay(leduc[["proteins_norm2"]])
pct_na_global <- round(100 * mean(is.na(mat_original)), 1)
message("Porcentaje global de NAs: ", pct_na_global, "%")

# --- 2. FUNCIÓN DE IMPUTACIÓN KNN ---
imputeKnnSCoPE2 <- function(object, i, name, k = 3) {
  mat <- assay(object[[i]])
  
  # Filtrar proteínas con >80% NAs — no imputables de forma fiable
  pNA_rows  <- rowMeans(is.na(mat))
  keep_rows <- pNA_rows <= 0.80
  
  message("Proteínas originales:                  ", nrow(mat))
  message("Proteínas eliminadas por >80% NAs:     ", sum(!keep_rows))
  message("Proteínas que entran en la imputación: ", sum(keep_rows))
  
  mat_filtered <- mat[keep_rows, ]
  
  # Transponer: KNN busca k células similares en el espacio proteómico
  # (protocolo SCoPE2 — vecinos entre células, no entre proteínas)
  mat_t     <- t(mat_filtered)
  imputed_t <- impute::impute.knn(mat_t, k = k)$data
  imputed   <- t(imputed_t)
  
  object <- addAssay(
    object,
    SingleCellExperiment(assays = list(imputed = imputed)),
    name = name
  )
  object <- addAssayLink(object, from = i, to = name)
  return(object)
}

# --- 3. EJECUCIÓN Y MEDICIÓN DE TIEMPO ---
message("\nIniciando imputación KNN (k=3)...")
t_knn <- system.time({
  leduc <- imputeKnnSCoPE2(leduc,
                           i    = "proteins_norm2",
                           name = "proteins_impd",
                           k    = 3)
})

message("\nImputación KNN completada.")
message("Tiempo de ejecución: ", round(t_knn["elapsed"], 1), " segundos")

if (!dir.exists("results/metrics")) dir.create("results/metrics",
                                               recursive = TRUE)
saveRDS(t_knn, "results/metrics/tiempo_imputacion_knn.rds")

# --- 4. VERIFICACIÓN ---
mat_imp <- assay(leduc[["proteins_impd"]])
n_nas   <- sum(is.na(mat_imp))
message("NAs restantes tras imputación: ", n_nas,
        ifelse(n_nas == 0, " ✓ imputación completa", ""))

# --- 5. RESUMEN FINAL ---
message("\n--- RESUMEN IMPUTACIÓN KNN ---")
message("Proteínas en matriz final: ", nrow(mat_imp))
message("Células en matriz final:   ", ncol(mat_imp))
message("NAs restantes:             ", n_nas)
message("Tiempo:                    ",
        round(t_knn["elapsed"], 1), " segundos")

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_imputed_knn.rds")
message("\nScript 11_1 finalizado.")
message("  data/leduc_imputed_knn.rds               → objeto con imputación KNN")
message("  results/metrics/tiempo_imputacion_knn.rds → tiempo registrado")