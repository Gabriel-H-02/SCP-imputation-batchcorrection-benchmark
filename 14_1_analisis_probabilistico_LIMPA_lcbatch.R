### ============================================================
### SCRIPT 14_1: PIPELINE LIMPA (CORRECCIÓN POR LCBATCH)
### ============================================================
### Variante de LIMPA que incluye lcbatch en el modelo lineal
### del análisis diferencial, corrigiendo el efecto de las 2
### sesiones de adquisición del instrumento.
### Diferencia respecto al Script 14:
###   Script 14:   design = ~ SampleType
###   Script 14_1: design = ~ SampleType + lcbatch
### La DPC y el DPC-Quant son idénticos — solo cambia el modelo
### del análisis diferencial al incluir lcbatch como covariable.
### INPUT:  data/leduc_final_processed.rds
### OUTPUT: data/matriz_final_limpa_lcbatch.rds
###         results/biomarcadores_limpa_vs_leduc_lcbatch.csv
###         results/plots/DPC_limpa_lcbatch.png
###         results/metrics/tiempo_limpa_lcbatch.rds
### ============================================================

library(limpa)
library(limma)
library(scp)

if (!dir.exists("results/plots"))   dir.create("results/plots",
                                               recursive = TRUE)
if (!dir.exists("results/metrics")) dir.create("results/metrics",
                                               recursive = TRUE)

# --- 1. CARGA ---
leduc    <- readRDS("data/leduc_final_processed.rds")
mat_nas  <- assay(leduc, "proteins_norm2")
prot_ids <- rownames(mat_nas)
message("Proteínas disponibles: ", length(prot_ids))
message("Porcentaje global de NAs: ",
        round(100 * mean(is.na(mat_nas)), 1), "%")

# --- 2. ESTIMACIÓN DE LA DPC ---
message("Estimando DPC para variante lcbatch...")
dpc_est <- dpc(mat_nas)

png("results/plots/DPC_limpa_lcbatch.png", width = 700, height = 500, res = 100)
plotDPC(dpc_est)
dev.off()
message("DPC estimada y guardada.")

# --- 3. DPC-QUANT + ANÁLISIS DIFERENCIAL Y MEDICIÓN DE TIEMPO ---
message("Ejecutando DPC-Quant (operación más costosa del pipeline)...")
t_limpa_lc <- system.time({
  
  y.protein <- dpcQuant(mat_nas, prot_ids, dpc = dpc_est)
  
  # --- 4. ANÁLISIS DIFERENCIAL CON LCBATCH ---
  # lcbatch entra como covariable técnica en el modelo lineal
  # para controlar el efecto de las 2 sesiones de adquisición
  metadata       <- as.data.frame(colData(leduc))
  metadata       <- metadata[colnames(mat_nas), ]
  design_lcbatch <- model.matrix(~ SampleType + lcbatch, data = metadata)
  
  fit_lcbatch <- dpcDE(y.protein, design_lcbatch)
  fit_lcbatch <- eBayes(fit_lcbatch)
  
})

message("Pipeline LIMPA lcbatch completado.")
message("Tiempo total (DPC-Quant + dpcDE): ",
        round(t_limpa_lc["elapsed"] / 3600, 2), " horas (",
        round(t_limpa_lc["elapsed"], 0), " segundos)")
saveRDS(t_limpa_lc, "results/metrics/tiempo_limpa_lcbatch.rds")

# --- 5. EXTRACCIÓN DE BIOMARCADORES ---
tabla_resultados_lcbatch <- topTable(fit_lcbatch, coef = 2,
                                     number = Inf, sort.by = "P")
tabla_resultados_lcbatch$is_significant <-
  tabla_resultados_lcbatch$adj.P.Val < 0.05

n_sig <- sum(tabla_resultados_lcbatch$is_significant, na.rm = TRUE)
message("Proteínas significativas (FDR < 0.05): ", n_sig)

write.csv(tabla_resultados_lcbatch,
          "results/biomarcadores_limpa_vs_leduc_lcbatch.csv")

# --- 6. EXPORTAR MATRIZ PARA PCA ---
mat_limpa_lcbatch <- y.protein$E
saveRDS(mat_limpa_lcbatch, "data/matriz_final_limpa_lcbatch.rds")

# --- 7. RESUMEN ---
message("\n--- RESUMEN PIPELINE LIMPA LCBATCH ---")
message("Proteínas analizadas:              ", length(prot_ids))
message("Proteínas significativas FDR < 5%: ", n_sig)
message("Tiempo total:                      ",
        round(t_limpa_lc["elapsed"] / 3600, 2), " horas")

message("\nScript 14_1 finalizado.")
message("  data/matriz_final_limpa_lcbatch.rds              → matriz para PCA")
message("  results/biomarcadores_limpa_vs_leduc_lcbatch.csv → tabla DE completa")
message("  results/metrics/tiempo_limpa_lcbatch.rds         → tiempo registrado")