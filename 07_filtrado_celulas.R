### ============================================================
### SCRIPT 07: FILTRADO DE CÉLULAS POR MEDIAN CV
### ============================================================
### Evalúa la calidad de cada célula mediante el coeficiente de
### variación mediano (Median CV) entre péptidos de la misma
### proteína. Las células cuyo CV es indistinguible del ruido de
### fondo (controles negativos) se eliminan del análisis.
### El umbral CV < 0.42 se establece en el punto de mínimo
### solapamiento entre células reales y canales vacíos.
### INPUT:  data/leduc_peptides_joined.rds
### OUTPUT: data/leduc_cells_filtered.rds
###         results/plots/02_separacion_cv_celulas.png
### ============================================================

library(scp)
library(dplyr)
library(ggplot2)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_peptides_joined.rds")
peptideAssays <- grep("^peptides_", names(leduc), value = TRUE)
message("Células totales antes del filtrado: ", ncol(leduc[["peptides"]]))

# --- 2. CÁLCULO DEL MEDIAN CV ---
# Para cada célula: agrupa péptidos por proteína y calcula la
# variabilidad interna. La normalización SCoPE2 es necesaria
# para que el CV sea comparable entre células de distinto tamaño.
# nobs = 3: solo se usan proteínas con al menos 3 péptidos detectados.
leduc <- medianCVperCell(leduc,
                         i           = peptideAssays,
                         groupBy     = "Leading.razor.protein.symbol",
                         nobs        = 3,
                         na.rm       = TRUE,
                         colDataName = "MedianCV",
                         norm        = "SCoPE2")
message("Paso 1: Cálculo de Median CV finalizado.")

# --- 3. VISUALIZACIÓN ---
p_cv <- colData(leduc) %>%
  data.frame() %>%
  filter(grepl("Mono|Mel|Neg", SampleType)) %>%
  mutate(control = ifelse(grepl("Neg", SampleType), "ctl", "sc")) %>%
  ggplot(aes(x = MedianCV, fill = control)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.42, linetype = "dashed",
             color = "red", linewidth = 0.9) +
  annotate("text", x = 0.44, y = Inf, vjust = 2, hjust = 0,
           label = "Umbral: CV < 0.42", color = "red", size = 3.5) +
  scale_fill_manual(
    values = c("ctl" = "black", "sc" = "purple2"),
    labels = c("ctl" = "Controles negativos (canales vacíos)",
               "sc"  = "Células individuales (Melanoma + Monocito)")
  ) +
  labs(
    title    = "QC Celular: Separación de células vs ruido de fondo",
    subtitle = "El umbral de CV < 0.42 elimina células indistinguibles del ruido",
    x        = "Median CV (Variabilidad interna por proteína)",
    y        = "Densidad",
    fill     = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_cv)

if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
ggsave("results/plots/02_separacion_cv_celulas.png",
       p_cv, width = 8, height = 6, dpi = 300)

# --- 4. FILTRADO ---
n_antes <- ncol(leduc[["peptides"]])

leduc <- subsetByColData(leduc,
                         !is.na(leduc$MedianCV) &
                           leduc$MedianCV < 0.42 &
                           grepl("Mono|Mel", leduc$SampleType))

n_despues <- ncol(leduc[["peptides"]])

# --- 5. RESUMEN ---
message("\n--- RESUMEN DEL FILTRADO CELULAR ---")
message("Células antes del filtrado:   ", n_antes)
message("Células tras el filtrado:     ", n_despues)
message("Células eliminadas:           ", n_antes - n_despues,
        " (", round(100 * (n_antes - n_despues) / n_antes, 1), "%)")
message("  (incluye NegControls, Carriers y células de baja calidad)")

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_cells_filtered.rds")
message("\nScript 07 finalizado.")
message("  data/leduc_cells_filtered.rds → ", n_despues, " células retenidas")
message("  results/plots/02_separacion_cv_celulas.png → figura de QC celular")