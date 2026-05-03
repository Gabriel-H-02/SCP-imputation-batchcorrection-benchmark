### ============================================================
### SCRIPT 03: FILTRADO DE ESPECTROS (PSM LEVEL)
### ============================================================
### Aplica el primer nivel de control de calidad sobre los
### espectros de masa (PSMs). Elimina contaminantes, péptidos
### reversos (decoys), espectros de baja pureza iónica (PIF) y
### espectros con baja confianza de identificación (dart_qval).
### INPUT:  data/leduc_raw.rds
### OUTPUT: data/leduc_psm_filtered.rds
###         results/plots/01_qc_pep_pif.png
### ============================================================

library(scp)
library(ggplot2)
library(patchwork)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_raw.rds")
message("PSMs antes del filtrado: cargando metadatos...")

# Extraemos los metadatos de filas para visualizar antes de filtrar
rd <- data.frame(rbindRowData(leduc, i = names(leduc)))
n_psm_antes <- nrow(rd)
message("Total de PSMs antes del filtrado: ", n_psm_antes)

# --- 2. VISUALIZACIÓN PRE-FILTRADO ---

# Fig A: Distribución de PEP (Posterior Error Probability de DART-ID)
# Un PEP bajo indica alta confianza de identificación.
# Corte en dart_qval < 0.01 → FDR controlado al 1%
p1 <- ggplot(rd, aes(x = dart_PEP)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.01, linetype = "dashed",
             color = "red", linewidth = 0.9) +
  annotate("text", x = 0.05, y = Inf, vjust = 2,
           label = "Corte: dart_qval < 0.01", color = "red", size = 3.5) +
  labs(
    title    = "Distribución de PEP (DART-ID)",
    subtitle = "Probabilidad de error posterior por espectro",
    x        = "PEP (Posterior Error Probability)",
    y        = "Frecuencia"
  ) +
  theme_minimal()

# Fig B: Distribución de PIF (Precursor Ion Fraction)
# El PIF mide la pureza del espectro: fracción de la señal del precursor
# que proviene del péptido de interés. PIF bajo indica co-aislamiento.
p2 <- ggplot(rd[!is.na(rd$PIF), ], aes(x = PIF)) +
  geom_histogram(bins = 50, fill = "darkorange", color = "white") +
  geom_vline(xintercept = 0.6, linetype = "dashed",
             color = "red", linewidth = 0.9) +
  annotate("text", x = 0.75, y = Inf, vjust = 2,
           label = "Corte: PIF > 0.6", color = "red", size = 3.5) +
  labs(
    title    = "Distribución de PIF (Precursor Ion Fraction)",
    subtitle = "Pureza del espectro por co-aislamiento de iones",
    x        = "PIF",
    y        = "Frecuencia"
  ) +
  theme_minimal()

panel_qc <- p1 + p2 +
  plot_annotation(
    title = "Control de Calidad Pre-filtrado — Nivel PSM",
    theme = theme(plot.title = element_text(size = 13, face = "bold",
                                            hjust = 0.5))
  )

print(panel_qc)

if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
ggsave("results/plots/01_qc_pep_pif.png",
       panel_qc, width = 12, height = 5, dpi = 300)

# --- 3. FILTRADO ---
# Criterios aplicados (reproduciendo el protocolo de Leduc et al.):
# - Eliminación de contaminantes conocidos (CON) y secuencias reversas (REV)
# - PIF > 0.6: garantiza que >60% de la señal proviene del péptido de interés
# - dart_qval < 0.01: controla el FDR de identificación al 1%
leduc <- filterFeatures(leduc,
                        ~ Potential.contaminant != "+" &
                          !grepl("CON", Proteins) &
                          Reverse != "+" &
                          !grepl("REV", Leading.razor.protein) &
                          (is.na(PIF) | PIF > 0.6) &
                          dart_qval < 0.01)

# --- 4. RESUMEN DEL FILTRADO ---
rd_post <- data.frame(rbindRowData(leduc, i = names(leduc)))
n_psm_despues <- nrow(rd_post)

message("\n--- RESUMEN DEL FILTRADO PSM ---")
message("PSMs antes:  ", n_psm_antes)
message("PSMs después: ", n_psm_despues)
message("PSMs eliminados: ", n_psm_antes - n_psm_despues,
        " (", round(100 * (n_psm_antes - n_psm_despues) / n_psm_antes, 1),
        "% del total)")

# --- 5. GUARDAR ---
saveRDS(leduc, file = "data/leduc_psm_filtered.rds")
message("\nScript 03 finalizado. Objeto guardado en data/leduc_psm_filtered.rds")
message("Figura guardada en results/plots/01_qc_pep_pif.png")