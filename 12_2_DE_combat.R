### ============================================================
### SCRIPT 12_2: CORRECCIÓN COMBAT Y ANÁLISIS DIFERENCIAL
### ============================================================
### Aplica ComBat para la corrección de efectos de lote y el
### análisis de expresión diferencial (limma + eBayes) sobre
### las dos variantes:
###   - ComBat lcbatch (2 niveles de lote grueso)
###   - ComBat Set (133 niveles de placa individual)
### Los resultados alimentan el Script 15_2 para completar la
### tabla de métricas del OE2 (N_sig y Overlap para ComBat).
### INPUT:  data/leduc_combat_lcbatch.rds
###         data/leduc_combat_SET.rds
### OUTPUT: results/biomarcadores_combat_lcbatch.csv
###         results/biomarcadores_combat_SET.csv
###         results/tables/TOP_20_BIOMARCADORES_COMBAT_LC.csv
###         results/tables/TOP_20_BIOMARCADORES_COMBAT_SET.csv
###         results/metrics/tiempo_combat_lcbatch.rds
###         results/metrics/tiempo_combat_SET.rds
### ============================================================

library(scp)
library(limma)

if (!dir.exists("results"))          dir.create("results", recursive = TRUE)
if (!dir.exists("results/tables"))   dir.create("results/tables",
                                                recursive = TRUE)
if (!dir.exists("results/metrics"))  dir.create("results/metrics",
                                                recursive = TRUE)

# --- FUNCIÓN AUXILIAR DE ANÁLISIS DIFERENCIAL ---
analisis_diferencial <- function(leduc, assay_name, label) {
  sce    <- getWithColData(leduc, assay_name)
  mat    <- assay(sce)
  meta   <- as.data.frame(colData(sce))
  design <- model.matrix(~ SampleType, data = meta)
  fit    <- lmFit(mat, design)
  fit    <- eBayes(fit)
  tabla  <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  tabla$is_significant <- tabla$adj.P.Val < 0.05
  n_sig  <- sum(tabla$is_significant, na.rm = TRUE)
  message(label, " — Proteínas significativas (FDR < 0.05): ", n_sig)
  return(tabla)
}

# ── 1. COMBAT LCBATCH ──────────────────────────────────────────────────────
message("Cargando ComBat lcbatch...")
leduc_lc <- readRDS("data/leduc_combat_lcbatch.rds")

assay_lc <- grep("combat", names(leduc_lc), value = TRUE,
                 ignore.case = TRUE)[1]
message("Assay usado: ", assay_lc)

# Medición de tiempo del análisis diferencial
message("Ejecutando análisis diferencial ComBat lcbatch...")
t_lc <- system.time({
  tabla_lc <- analisis_diferencial(leduc_lc, assay_lc, "ComBat lcbatch")
})
message("Tiempo ComBat lcbatch DE: ", round(t_lc["elapsed"], 1), " segundos")
saveRDS(t_lc, "results/metrics/tiempo_combat_lcbatch.rds")

write.csv(tabla_lc, "results/biomarcadores_combat_lcbatch.csv")
top20_lc <- head(tabla_lc[order(tabla_lc$adj.P.Val), ], 20)
write.csv(top20_lc, "results/tables/TOP_20_BIOMARCADORES_COMBAT_LC.csv")

# ── 2. COMBAT SET ──────────────────────────────────────────────────────────
message("\nCargando ComBat Set...")
leduc_set <- readRDS("data/leduc_combat_SET.rds")

assay_set <- grep("combat", names(leduc_set), value = TRUE,
                  ignore.case = TRUE)[1]
message("Assay usado: ", assay_set)

message("Ejecutando análisis diferencial ComBat Set...")
t_set <- system.time({
  tabla_set <- analisis_diferencial(leduc_set, assay_set, "ComBat Set")
})
message("Tiempo ComBat Set DE: ", round(t_set["elapsed"], 1), " segundos")
saveRDS(t_set, "results/metrics/tiempo_combat_SET.rds")

write.csv(tabla_set, "results/biomarcadores_combat_SET.csv")
top20_set <- head(tabla_set[order(tabla_set$adj.P.Val), ], 20)
write.csv(top20_set, "results/tables/TOP_20_BIOMARCADORES_COMBAT_SET.csv")

# ── 3. RESUMEN ─────────────────────────────────────────────────────────────
message("\n--- RESUMEN ANÁLISIS DIFERENCIAL COMBAT ---")
message("ComBat lcbatch — proteínas FDR < 0.05: ",
        sum(tabla_lc$is_significant,  na.rm = TRUE),
        " | Tiempo: ", round(t_lc["elapsed"], 1), "s")
message("ComBat Set     — proteínas FDR < 0.05: ",
        sum(tabla_set$is_significant, na.rm = TRUE),
        " | Tiempo: ", round(t_set["elapsed"], 1), "s")

message("\nScript 12_2 finalizado.")
message("  results/biomarcadores_combat_lcbatch.csv → DE lcbatch")
message("  results/biomarcadores_combat_SET.csv     → DE Set")