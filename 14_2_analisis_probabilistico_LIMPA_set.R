### ============================================================
### SCRIPT 14_2: PIPELINE LIMPA (CORRECCIÓN POR SET)
### ============================================================
### Variante de LIMPA que incluye Set en el modelo lineal del
### análisis diferencial, corrigiendo el efecto de las 133
### placas individuales de adquisición. Es la variante con
### mayor resolución de corrección de lote y la que produce
### el mayor número de proteínas significativas del benchmark.
### Diferencia respecto a los Scripts 14 y 14_1:
###   Script 14:   design = ~ SampleType
###   Script 14_1: design = ~ SampleType + lcbatch  (2 niveles)
###   Script 14_2: design = ~ SampleType + Set      (133 niveles)
### INPUT:  data/leduc_final_processed.rds
### OUTPUT: data/matriz_final_limpa_SET.rds
###         results/biomarcadores_limpa_RESOLUCION_SET.csv
###         results/plots/DPC_limpa_set.png
###         results/metrics/tiempo_limpa_SET.rds
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
message("Estimando DPC para variante Set...")
dpc_est <- dpc(mat_nas)

png("results/plots/DPC_limpa_set.png", width = 700, height = 500, res = 100)
plotDPC(dpc_est)
dev.off()
message("DPC estimada y guardada.")

# --- 3. DPC-QUANT + ANÁLISIS DIFERENCIAL Y MEDICIÓN DE TIEMPO ---
message("Ejecutando DPC-Quant (operación más costosa del pipeline)...")
t_limpa_set <- system.time({
  
  y.protein <- dpcQuant(mat_nas, prot_ids, dpc = dpc_est)
  
  # --- 4. ANÁLISIS DIFERENCIAL CON SET ---
  # Set entra como covariable técnica en el modelo lineal
  # para controlar el efecto de las 133 placas individuales.
  # Es el modelo más paramétrico de los tres — mayor resolución
  # pero también mayor coste computacional en el dpcDE.
  metadata   <- as.data.frame(colData(leduc))
  metadata   <- metadata[colnames(mat_nas), ]
  design_set <- model.matrix(~ SampleType + Set, data = metadata)
  
  message("Ejecutando dpcDE con corrección por ",
          length(unique(metadata$Set)), " placas (Set)...")
  
  fit_set <- dpcDE(y.protein, design_set)
  fit_set <- eBayes(fit_set)
  
})

message("Pipeline LIMPA Set completado.")
message("Tiempo total (DPC-Quant + dpcDE): ",
        round(t_limpa_set["elapsed"] / 3600, 2), " horas (",
        round(t_limpa_set["elapsed"], 0), " segundos)")
saveRDS(t_limpa_set, "results/metrics/tiempo_limpa_SET.rds")

# --- 5. EXTRACCIÓN DE BIOMARCADORES ---
tabla_resultados_set <- topTable(fit_set, coef = 2,
                                 number = Inf, sort.by = "P")
tabla_resultados_set$is_significant <-
  tabla_resultados_set$adj.P.Val < 0.05

n_sig <- sum(tabla_resultados_set$is_significant, na.rm = TRUE)
message("Proteínas significativas con resolución SET (FDR < 0.05): ", n_sig)

write.csv(tabla_resultados_set,
          "results/biomarcadores_limpa_RESOLUCION_SET.csv")

# --- 6. EXPORTAR MATRIZ PARA PCA ---
mat_limpa_set <- y.protein$E
saveRDS(mat_limpa_set, "data/matriz_final_limpa_SET.rds")

# --- 7. RESUMEN ---
message("\n--- RESUMEN PIPELINE LIMPA SET ---")
message("Proteínas analizadas:              ", length(prot_ids))
message("Proteínas significativas FDR < 5%: ", n_sig)
message("Niveles de Set (placas):           ",
        length(unique(metadata$Set)))
message("Tiempo total:                      ",
        round(t_limpa_set["elapsed"] / 3600, 2), " horas")

message("\nScript 14_2 finalizado.")
message("  data/matriz_final_limpa_SET.rds                  → matriz para PCA")
message("  results/biomarcadores_limpa_RESOLUCION_SET.csv   → tabla DE completa")
message("  results/metrics/tiempo_limpa_SET.rds             → tiempo registrado")