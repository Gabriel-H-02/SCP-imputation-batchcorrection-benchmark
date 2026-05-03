### ============================================================
### SCRIPT 14: PIPELINE LIMPA (SIN CORRECCIÓN DE LOTE)
### ============================================================
### Implementa el pipeline probabilístico limpa en su variante
### base, sin corrección de efectos de lote en el modelo lineal.
### Tres pasos secuenciales:
###   (1) Estimación de la DPC: ajusta un modelo logístico que
###       relaciona la probabilidad de detección con la intensidad.
###   (2) DPC-Quant: cuantificación probabilística que integra
###       observados y no observados, propagando la incertidumbre.
###   (3) dpcDE + eBayes: análisis diferencial con pesos de
###       precisión derivados de DPC-Quant (pipeline vooma).
### INPUT:  data/leduc_final_processed.rds
### OUTPUT: data/matriz_final_limpa.rds
###         results/biomarcadores_limpa_vs_leduc.csv
###         results/plots/DPC_limpa_base.png
###         results/metrics/tiempo_limpa_base.rds
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
message("Estimando la Curva de Probabilidad de Detección (DPC)...")
dpc_est <- dpc(mat_nas)

# plotDPC() devuelve un objeto gráfico base de R, no ggplot
# por eso usamos png() en lugar de ggsave()
png("results/plots/DPC_limpa_base.png", width = 700, height = 500, res = 100)
plotDPC(dpc_est)
dev.off()
message("DPC estimada y guardada.")

# --- 3. CUANTIFICACIÓN PROBABILÍSTICA (DPC-Quant) Y MEDICIÓN DE TIEMPO ---
# Este paso es el más costoso computacionalmente (~2.5h en M1 16GB)
message("Ejecutando DPC-Quant (operación más costosa del pipeline)...")
t_limpa <- system.time({
  
  y.protein <- dpcQuant(mat_nas, prot_ids, dpc = dpc_est)
  
  # --- 4. ANÁLISIS DIFERENCIAL ---
  metadata <- as.data.frame(colData(leduc))
  metadata <- metadata[colnames(mat_nas), ]
  design   <- model.matrix(~ SampleType, data = metadata)
  
  fit <- dpcDE(y.protein, design)
  fit <- eBayes(fit)
  
})

message("Pipeline LIMPA base completado.")
message("Tiempo total (DPC-Quant + dpcDE): ",
        round(t_limpa["elapsed"] / 3600, 2), " horas (",
        round(t_limpa["elapsed"], 0), " segundos)")
saveRDS(t_limpa, "results/metrics/tiempo_limpa_base.rds")

# --- 5. EXTRACCIÓN DE BIOMARCADORES ---
tabla_resultados <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
tabla_resultados$is_significant <- tabla_resultados$adj.P.Val < 0.05

n_sig <- sum(tabla_resultados$is_significant, na.rm = TRUE)
message("Proteínas significativas (FDR < 0.05): ", n_sig)

write.csv(tabla_resultados, "results/biomarcadores_limpa_vs_leduc.csv")

# --- 6. EXPORTAR MATRIZ PARA PCA ---
mat_limpa <- y.protein$E
saveRDS(mat_limpa, "data/matriz_final_limpa.rds")

# --- 7. RESUMEN ---
message("\n--- RESUMEN PIPELINE LIMPA BASE ---")
message("Proteínas analizadas:              ", length(prot_ids))
message("Proteínas significativas FDR < 5%: ", n_sig)
message("Tiempo total:                      ",
        round(t_limpa["elapsed"] / 3600, 2), " horas")

message("\nScript 14 finalizado.")
message("  data/matriz_final_limpa.rds                → matriz para PCA")
message("  results/biomarcadores_limpa_vs_leduc.csv   → tabla DE completa")
message("  results/metrics/tiempo_limpa_base.rds      → tiempo registrado")