### ============================================================
### SCRIPT 02: PREPARACIÓN DEL BENCHMARK
### ============================================================
### Separa el dataset original en dos componentes:
### (A) Los datos crudos de PSM sobre los que se aplicará el
###     pipeline propio desde cero.
### (B) Las matrices ya procesadas por Leduc et al. que servirán
###     como referencia de validación al final del TFM.
### Esta separación garantiza que el pipeline propio no usa
### ningún resultado intermedio de los autores originales.
### INPUT:  data/leduc_original.rds
### OUTPUT: data/leduc_raw.rds
###         data/leduc_author_solutions.RData
### ============================================================

library(scp)
library(ggplot2)
library(dplyr)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_original.rds")
message("Objeto original cargado: ", length(names(leduc)), " assays.")

# --- 2. EXTRAER SOLUCIONES DE LOS AUTORES ---
# Guardamos las 4 capas ya procesadas por Leduc et al.:
# - peptides: matriz de péptidos agregada
# - peptides_log: matriz de péptidos en escala log2
# - proteins_norm2: proteínas normalizadas (antes del batch correction)
# - proteins_processed: proteínas procesadas finales (referencia gold standard)
proteins_norm_leduc      <- leduc[["proteins_norm2"]]
proteins_processed_leduc <- leduc[["proteins_processed"]]

message("Soluciones de los autores extraídas:")
message("  - proteins_norm2:      ", nrow(proteins_norm_leduc),
        " proteínas x ", ncol(proteins_norm_leduc), " células")
message("  - proteins_processed:  ", nrow(proteins_processed_leduc),
        " proteínas x ", ncol(proteins_processed_leduc), " células")

# --- 3. RESET: QUEDARSE SOLO CON LOS PSM CRUDOS ---
# Eliminamos los 4 assays procesados (posiciones 135-138)
# El objeto vuelve a tener solo los 134 lotes de PSM originales
n_antes  <- length(names(leduc))
leduc    <- leduc[, , -(135:138)]
n_despues <- length(names(leduc))

message("\nObjeto reseteado a datos crudos:")
message("  Assays antes:  ", n_antes)
message("  Assays después: ", n_despues, " (134 lotes PSM crudos)")

# --- 4. FIGURA DIAGNÓSTICA ---
# Visualizamos el esquema de separación para el manuscrito
diseño_df <- data.frame(
  Componente = c("Pipeline propio\n(este TFM)",
                 "Referencia autores\n(Leduc et al.)"),
  Tipo       = c("Datos crudos PSM", "Matrices procesadas"),
  N_assays   = c(n_despues, 4),
  Color      = c("#377EB8", "#E41A1C")
)

p_esquema <- ggplot(diseño_df,
                    aes(x = Componente, y = N_assays, fill = Componente)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = paste0(N_assays, " assays")),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  scale_y_continuous(limits = c(0, 145)) +
  labs(
    title    = "Separación del dataset para el diseño del benchmarking",
    subtitle = "Los datos crudos y las soluciones de los autores se almacenan por separado",
    x        = NULL,
    y        = "Número de assays"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11))

print(p_esquema)

if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
ggsave("results/plots/00_diseno_benchmark.png",
       p_esquema, width = 7, height = 5, dpi = 300)

# --- 5. GUARDAR ---
# A) Datos crudos para el pipeline propio
saveRDS(leduc, file = "data/leduc_raw.rds")

# B) Soluciones de los autores para la validación final
save(proteins_norm_leduc, proteins_processed_leduc,
     file = "data/leduc_author_solutions.RData")

message("\nScript 02 finalizado.")
message("  data/leduc_raw.rds                  → pipeline propio")
message("  data/leduc_author_solutions.RData   → referencia autores")
message("  results/plots/00_diseno_benchmark.png → figura de diseño")