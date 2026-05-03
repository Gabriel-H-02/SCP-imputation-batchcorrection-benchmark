### ============================================================
### SCRIPT 04: FILTRADO POR SAMPLE-TO-CARRIER RATIO (SCR)
### ============================================================
### Calcula el Sample-to-Carrier Ratio (SCR) por PSM y elimina
### aquellos espectros donde la señal de la célula individual
### supera el 5% del carrier proteome, indicando posible
### contaminación cruzada entre canales TMT o co-aislamiento.
###
### Fundamento: En experimentos pSCoPE con carrier de 100
### equivalentes celulares, el SCR esperado es ~1/100 = 0.01.
### Un SCR > 0.05 indica que la señal está dominada por el
### carrier y no refleja fielmente la biología celular.
###
### INPUT:  data/leduc_psm_filtered.rds
### OUTPUT: data/leduc_scr_filtered.rds
###         results/plots/01_distribucion_scr.png
### ============================================================

library(scp)
library(ggplot2)
library(dplyr)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_psm_filtered.rds")
message("Objeto cargado tras filtrado PSM.")

# --- 2. CÁLCULO DEL SCR ---
# computeSCR calcula el ratio medio entre la señal de la célula
# individual y la señal del carrier proteome por cada PSM.
leduc <- computeSCR(leduc,
                    i            = names(leduc),
                    colvar       = "SampleType",
                    samplePattern  = "Mel|Mono",
                    carrierPattern = "Carrier",
                    sampleFUN    = "mean",
                    rowDataName  = "MeanSCR")

message("SCR calculado para todos los PSMs.")

# --- 3. VISUALIZACIÓN ---
rd_scr <- rbindRowData(leduc, i = names(leduc)) %>%
  data.frame() %>%
  filter(!is.na(MeanSCR) & !is.infinite(MeanSCR) & MeanSCR > 0)

n_total   <- nrow(rd_scr)
n_eliminar <- sum(rd_scr$MeanSCR >= 0.05)

p_scr <- ggplot(rd_scr, aes(x = MeanSCR)) +
  geom_histogram(bins = 50, fill = "purple", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.05, linetype = "dashed",
             color = "red", linewidth = 0.9) +
  annotate("text", x = 0.08, y = Inf, vjust = 2,
           label = paste0("Corte: SCR < 0.05\nEliminados: ",
                          n_eliminar, " PSMs (",
                          round(100 * n_eliminar / n_total, 1), "%)"),
           color = "red", size = 3.2, hjust = 0) +
  scale_x_log10() +
  labs(
    title    = "Distribución de SCR (Sample-to-Carrier Ratio)",
    subtitle = "Corte estricto en 0.05 (5%) — escala logarítmica",
    x        = "Mean SCR",
    y        = "Frecuencia"
  ) +
  theme_minimal()

print(p_scr)

if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
ggsave("results/plots/01_distribucion_scr.png",
       p_scr, width = 8, height = 6, dpi = 300)

# --- 4. APLICAR FILTRO ---
n_antes <- sum(sapply(names(leduc), function(i) nrow(assay(leduc[[i]]))))

leduc <- filterFeatures(leduc,
                        ~ !is.na(MeanSCR) &
                          !is.infinite(MeanSCR) &
                          MeanSCR < 0.05)

# --- 5. RESUMEN ---
message("\n--- RESUMEN DEL FILTRADO SCR ---")
message("PSMs con SCR calculado: ", n_total)
message("PSMs eliminados (SCR >= 0.05): ", n_eliminar,
        " (", round(100 * n_eliminar / n_total, 1), "%)")
message("PSMs retenidos: ", n_total - n_eliminar)

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_scr_filtered.rds")
message("\nScript 04 finalizado. Objeto guardado en data/leduc_scr_filtered.rds")
message("Figura guardada en results/plots/01_distribucion_scr.png")