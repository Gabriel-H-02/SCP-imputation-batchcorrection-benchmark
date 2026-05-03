### ============================================================
### SCRIPT 01: CARGA Y EXPLORACIÓN DEL DATASET
### ============================================================
### Descarga el dataset de Leduc et al. (2022) desde el paquete
### scpdata, explora su estructura jerárquica y genera una figura
### diagnóstica de la distribución de tipos de muestra.
### INPUT:  Descarga automática vía scpdata::leduc2022_pSCoPE()
### OUTPUT: data/leduc_original.rds
###         results/plots/00_exploracion_dataset.png
### ============================================================

library(scp)
library(scpdata)
library(dplyr)
library(ggplot2)

if (!dir.exists("data"))           dir.create("data")
if (!dir.exists("results/plots"))  dir.create("results/plots", recursive = TRUE)

# --- 1. DESCARGA DEL DATASET ---
message("Descargando dataset leduc2022_pSCoPE desde scpdata...")
leduc_data <- leduc2022_pSCoPE()
message("Tamaño del objeto en memoria:")
print(object.size(leduc_data), units = "Mb")

# --- 2. EXPLORACIÓN DE LA ESTRUCTURA ---
message("\nEl objeto contiene ", length(names(leduc_data)), " capas (assays).")
message("Nombres de los assays:")
print(names(leduc_data))

# --- 3. EXPLORACIÓN DE METADATOS ---
metadatos <- as.data.frame(colData(leduc_data))

# Distribución de tipos de muestra
resumen_muestras <- as.data.frame(table(metadatos$SampleType))
colnames(resumen_muestras) <- c("Tipo", "N")
message("\nDistribución de muestras original:")
print(resumen_muestras)

# Número de lotes
n_sets <- length(unique(metadatos$Set))
message("Datos distribuidos en ", n_sets, " lotes (Sets).")

# --- 4. FIGURA DIAGNÓSTICA DE DISTRIBUCIÓN ---
p_dist <- ggplot(resumen_muestras, aes(x = reorder(Tipo, -N), y = N,
                                       fill = Tipo)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = N), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Melanoma"  = "#377EB8",
    "Monocyte"  = "#E41A1C",
    "Carrier"   = "#FF7F00",
    "Reference" = "#4DAF4A",
    "Negative"  = "#999999"
  )) +
  labs(
    title    = "Distribución de tipos de muestra — Dataset Leduc et al. (2022)",
    subtitle = paste0(n_sets, " lotes de adquisición (Sets) | ",
                      nrow(metadatos), " canales TMT totales"),
    x        = "Tipo de muestra",
    y        = "Número de canales"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

print(p_dist)
ggsave("results/plots/00_exploracion_dataset.png",
       p_dist, width = 8, height = 5, dpi = 300)

# --- 5. INSPECCIÓN DE UNA CAPA CRUDA ---
message("\nPrimeras filas de la primera placa (PSM crudos):")
primera_placa <- assay(leduc_data[[1]])
print(head(primera_placa[, 1:5]))

# --- 6. GUARDAR ---
saveRDS(leduc_data, file = "data/leduc_original.rds")
message("\nScript 01 completado. Objeto guardado en data/leduc_original.rds")
message("Figura guardada en results/plots/00_exploracion_dataset.png")