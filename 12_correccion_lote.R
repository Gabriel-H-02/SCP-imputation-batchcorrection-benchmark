### ============================================================
### SCRIPT 12: CORRECCIÓN DE EFECTO LOTE (removeBatchEffect)
### ============================================================
### Aplica la corrección de lote del pipeline original de Leduc
### et al. usando removeBatchEffect de limma. Este es el BASELINE
### de referencia del benchmarking.
### Se corrigen dos fuentes de variabilidad técnica:
###   - lcbatch: las 2 sesiones de adquisición del instrumento
###   - Channel: el canal TMT de cada célula dentro de su set
### La variabilidad biológica (SampleType) se protege mediante
### el modelo de diseño antes de la corrección.
### INPUT:  data/leduc_imputed_knn.rds
### OUTPUT: data/leduc_final_processed.rds
###         results/metrics/tiempo_batch_correction_knn.rds
### ============================================================

library(scp)
library(limma)
library(matrixStats)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_imputed_knn.rds")
message("Objeto KNN cargado.")

# --- 2. EXTRACCIÓN ---
sce <- getWithColData(leduc, "proteins_impd")
message("Dimensiones: ", nrow(sce), " proteínas x ", ncol(sce), " células")

# --- 3. CORRECCIÓN DE LOTE Y MEDICIÓN DE TIEMPO ---
message("\nAplicando removeBatchEffect...")
t_batch <- system.time({
  
  model <- model.matrix(~ SampleType, data = colData(sce))
  
  assay(sce) <- removeBatchEffect(
    x      = assay(sce),
    batch  = sce$lcbatch,
    batch2 = sce$Channel,
    design = model
  )
  
})
message("Corrección removeBatchEffect completada.")
message("Tiempo de ejecución: ", round(t_batch["elapsed"], 1), " segundos")

if (!dir.exists("results/metrics")) dir.create("results/metrics",
                                               recursive = TRUE)
saveRDS(t_batch, "results/metrics/tiempo_batch_correction_knn.rds")

# --- 4. DEVOLVER AL OBJETO QFEATURES ---
leduc <- addAssay(leduc, y = sce, name = "proteins_batchC")
leduc <- addAssayLinkOneToOne(leduc,
                              from = "proteins_impd",
                              to   = "proteins_batchC")

# --- 5. NORMALIZACIÓN FINAL ---
# Centrado por columna (célula)
col_meds <- colMedians(assay(leduc[["proteins_batchC"]]), na.rm = TRUE)
leduc <- sweep(leduc,
               i      = "proteins_batchC",
               name   = "proteins_batchC_norm1",
               MARGIN = 2,
               FUN    = "-",
               STATS  = col_meds)

# Centrado por fila (proteína)
leduc <- sweep(leduc,
               i      = "proteins_batchC_norm1",
               name   = "proteins_processed",
               MARGIN = 1,
               FUN    = "-",
               STATS  = rowMedians(assay(leduc[["proteins_batchC_norm1"]]),
                                   na.rm = TRUE))
message("Normalización final completada.")

# --- 6. OPTIMIZACIÓN DE MEMORIA ---
leduc <- removeAssay(leduc, "proteins_batchC")
leduc <- removeAssay(leduc, "proteins_batchC_norm1")
gc()
message("Assays intermedios eliminados de memoria.")

# --- 7. RESUMEN ---
message("\n--- RESUMEN CORRECCIÓN DE LOTE (KNN baseline) ---")
message("Método:   removeBatchEffect (lcbatch + Channel)")
message("Tiempo:   ", round(t_batch["elapsed"], 1), " segundos")
message("Proteínas procesadas: ",
        nrow(assay(leduc[["proteins_processed"]])))

# --- 8. GUARDAR ---
saveRDS(leduc, file = "data/leduc_final_processed.rds")
message("\nScript 12 finalizado. Pipeline baseline (Leduc et al.) completado.")
message("  data/leduc_final_processed.rds               → objeto baseline")
message("  results/metrics/tiempo_batch_correction_knn.rds → tiempo registrado")