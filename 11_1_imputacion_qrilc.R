### ============================================================
### SCRIPT 11_1: IMPUTACIÓN QRILC (MÉTODO CLÁSICO MNAR)
### ============================================================
### Implementa la imputación QRILC (Quantile Regression Imputation
### of Left-Censored data) como segundo método de control clásico
### junto a KNN. QRILC es específico para datos MNAR en proteómica
### y es el contrapunto natural de LIMPA en el benchmarking.
###
### Diferencia clave vs KNN:
###   KNN (MAR): imputa a partir de células similares, ignorando
###              el mecanismo de ausencia.
###   QRILC (MNAR): modela los NAs como datos censurados por la
###                 izquierda (baja abundancia) e imputa desde la
###                 cola izquierda de la distribución estimada.
###
### Diferencia clave vs LIMPA:
###   QRILC asigna valores fijos sin propagar la incertidumbre.
###   LIMPA modela probabilísticamente la ausencia y propaga la
###   incertidumbre al análisis diferencial mediante vooma.
###
### INPUT:  data/leduc_proteins_final.rds
### OUTPUT: data/leduc_imputed_qrilc.rds
###         results/metrics/tiempo_imputacion_qrilc.rds
### ============================================================

library(scp)
library(SingleCellExperiment)
library(MsCoreUtils)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_proteins_final.rds")
message("Proteínas disponibles: ", nrow(leduc[["proteins_norm2"]]))
message("Células disponibles:   ", ncol(leduc[["proteins_norm2"]]))

# --- 2. FUNCIÓN DE IMPUTACIÓN QRILC ---
imputeQRILC <- function(object, i, name) {
  mat <- assay(object[[i]])
  
  # Filtrar proteínas con >80% NAs — no imputables de forma fiable
  pNA_rows  <- rowMeans(is.na(mat))
  keep_rows <- pNA_rows <= 0.80
  
  message("Proteínas originales:                  ", nrow(mat))
  message("Proteínas eliminadas por >80% NAs:     ", sum(!keep_rows))
  message("Proteínas que entran en la imputación: ", sum(keep_rows))
  
  mat_filtered <- mat[keep_rows, ]
  
  # QRILC estima los valores ausentes a partir de la distribución
  # truncada de cada proteína, modelando los NAs como observaciones
  # por debajo del límite de detección.
  mat_imputed <- MsCoreUtils::impute_matrix(mat_filtered, method = "QRILC")
  
  object <- addAssay(
    object,
    SingleCellExperiment(assays = list(imputed = mat_imputed)),
    name = name
  )
  object <- addAssayLink(object, from = i, to = name)
  return(object)
}

# --- 3. EJECUCIÓN Y MEDICIÓN DE TIEMPO ---
message("\nIniciando imputación QRILC...")
t_qrilc <- system.time({
  leduc <- imputeQRILC(leduc,
                       i    = "proteins_norm2",
                       name = "proteins_impd_qrilc")
})

message("\nImputación QRILC completada.")
message("Tiempo de ejecución: ", round(t_qrilc["elapsed"], 1), " segundos")

if (!dir.exists("results/metrics")) dir.create("results/metrics",
                                               recursive = TRUE)
saveRDS(t_qrilc, "results/metrics/tiempo_imputacion_qrilc.rds")

# --- 4. VERIFICACIÓN ---
mat_imp         <- assay(leduc[["proteins_impd_qrilc"]])
n_nas_restantes <- sum(is.na(mat_imp))
message("NAs restantes tras imputación: ", n_nas_restantes,
        ifelse(n_nas_restantes == 0, " ✓ imputación completa", ""))

# --- 5. RESUMEN FINAL ---
message("\n--- RESUMEN IMPUTACIÓN QRILC ---")
message("Proteínas en matriz final: ", nrow(mat_imp))
message("Células en matriz final:   ", ncol(mat_imp))
message("NAs restantes:             ", n_nas_restantes)
message("Tiempo:                    ",
        round(t_qrilc["elapsed"], 1), " segundos")

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_imputed_qrilc.rds")
message("\nScript 11_2 finalizado.")
message("  data/leduc_imputed_qrilc.rds                → objeto con imputación QRILC")
message("  results/metrics/tiempo_imputacion_qrilc.rds  → tiempo registrado")