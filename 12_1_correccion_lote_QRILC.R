### ============================================================
### SCRIPT 12_1: CORRECCIÓN DE LOTE SOBRE MATRIZ QRILC
### ============================================================
### Aplica removeBatchEffect sobre la matriz imputada por QRILC,
### siguiendo el mismo protocolo que el pipeline baseline KNN.
### Se corrigen las mismas fuentes de variabilidad técnica:
###   - lcbatch: las 2 sesiones de adquisición del instrumento
###   - Channel: el canal TMT de cada célula dentro de su set
### INPUT:  data/leduc_imputed_qrilc.rds
### OUTPUT: data/leduc_qrilc_processed.rds
###         results/metrics/tiempo_batch_correction_qrilc.rds
### ============================================================

library(scp)
library(limma)
library(matrixStats)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_imputed_qrilc.rds")
message("Objeto QRILC cargado.")

# --- 2. EXTRACCIÓN ---
sce <- getWithColData(leduc, "proteins_impd_qrilc")
message("Dimensiones: ", nrow(sce), " proteínas x ", ncol(sce), " células")

# --- 3. CORRECCIÓN DE LOTE Y MEDICIÓN DE TIEMPO ---
message("\nAplicando removeBatchEffect sobre QRILC...")
t_batch <- system.time({
  
  model <- model.matrix(~ SampleType, data = colData(sce))
  
  assay(sce) <- removeBatchEffect(
    x      = assay(sce),
    batch  = sce$lcbatch,
    batch2 = sce$Channel,
    design = model
  )
  
})
message("Corrección removeBatchEffect sobre QRILC completada.")
message("Tiempo de ejecución: ", round(t_batch["elapsed"], 1), " segundos")

if (!dir.exists("results/metrics")) dir.create("results/metrics",
                                               recursive = TRUE)
saveRDS(t_batch, "results/metrics/tiempo_batch_correction_qrilc.rds")

# --- 4. DEVOLVER AL OBJETO QFEATURES ---
leduc <- addAssay(leduc, y = sce, name = "proteins_qrilc_batchC")
leduc <- addAssayLinkOneToOne(leduc,
                              from = "proteins_impd_qrilc",
                              to   = "proteins_qrilc_batchC")

# --- 5. NORMALIZACIÓN FINAL ---
col_meds <- colMedians(assay(leduc[["proteins_qrilc_batchC"]]), na.rm = TRUE)
leduc <- sweep(leduc,
               i      = "proteins_qrilc_batchC",
               name   = "proteins_qrilc_norm1",
               MARGIN = 2,
               FUN    = "-",
               STATS  = col_meds)

leduc <- sweep(leduc,
               i      = "proteins_qrilc_norm1",
               name   = "proteins_qrilc_processed",
               MARGIN = 1,
               FUN    = "-",
               STATS  = rowMedians(assay(leduc[["proteins_qrilc_norm1"]]),
                                   na.rm = TRUE))
message("Normalización final completada.")

# --- 6. LIMPIEZA DE MEMORIA ---
leduc <- removeAssay(leduc, "proteins_qrilc_batchC")
leduc <- removeAssay(leduc, "proteins_qrilc_norm1")
gc()

# --- 7. RESUMEN ---
message("\n--- RESUMEN CORRECCIÓN DE LOTE (QRILC) ---")
message("Método:   removeBatchEffect (lcbatch + Channel)")
message("Tiempo:   ", round(t_batch["elapsed"], 1), " segundos")
message("Proteínas procesadas: ",
        nrow(assay(leduc[["proteins_qrilc_processed"]])))

# --- 8. GUARDAR ---
saveRDS(leduc, file = "data/leduc_qrilc_processed.rds")
message("\nScript 12_1 finalizado.")
message("  data/leduc_qrilc_processed.rds                  → objeto QRILC procesado")
message("  results/metrics/tiempo_batch_correction_qrilc.rds → tiempo registrado")