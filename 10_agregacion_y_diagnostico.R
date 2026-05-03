### ============================================================
### SCRIPT 10: AGREGACIÓN A PROTEÍNAS Y DIAGNÓSTICO DE NAs
### ============================================================
### Agrega péptidos a proteínas, aplica normalización final y
### genera un diagnóstico del porcentaje de valores faltantes
### a nivel proteico.
### INPUT:  data/leduc_peptides_ready.rds
### OUTPUT: data/leduc_proteins_final.rds
###         results/plots/03_diagnostico_nas_proteinas.png
### ============================================================

library(scp)
library(matrixStats)
library(ggplot2)
library(dplyr)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_peptides_ready.rds")
message("Péptidos disponibles: ", nrow(leduc[["peptides_log"]]))

# --- 2. AGREGACIÓN A NIVEL DE PROTEÍNA ---
leduc <- aggregateFeatures(leduc,
                           i     = "peptides_log",
                           name  = "proteins",
                           fcol  = "Leading.razor.protein.symbol",
                           fun   = matrixStats::colMedians,
                           na.rm = TRUE)
message("Proteínas cuantificadas: ", nrow(leduc[["proteins"]]))

# --- 3. NORMALIZACIÓN FINAL (Escala logarítmica) ---

# Centrado por columna (célula): restamos la mediana de cada célula
col_meds_prot <- colMedians(assay(leduc[["proteins"]]), na.rm = TRUE)
leduc <- sweep(leduc,
               i      = "proteins",
               name   = "proteins_norm1",
               MARGIN = 2,
               FUN    = "-",
               STATS  = col_meds_prot)

# Centrado por fila (proteína): restamos la mediana de cada proteína
leduc <- sweep(leduc,
               i      = "proteins_norm1",
               name   = "proteins_norm2",
               MARGIN = 1,
               FUN    = "-",
               STATS  = rowMedians(assay(leduc[["proteins_norm1"]]),
                                   na.rm = TRUE))

message("Agregación y normalización finalizadas.")

# --- 4. OPTIMIZACIÓN DE MEMORIA ---
leduc <- removeAssay(leduc, "proteins")
leduc <- removeAssay(leduc, "proteins_norm1")
gc()
message("Assays intermedios eliminados de memoria.")

# --- 5. DIAGNÓSTICO DE VALORES FALTANTES ---
nna_data         <- nNA(leduc, "proteins_norm2")$nNAcols
nna_data$pNA_pct <- nna_data$pNA * 100

p_nas <- nna_data %>%
  ggplot(aes(x = pNA_pct)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(nna_data$pNA_pct),
             linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text",
           x = mean(nna_data$pNA_pct) + 2, y = Inf, vjust = 2,
           label = paste0("Media: ", round(mean(nna_data$pNA_pct), 1), "%"),
           color = "red", size = 3.5) +
  labs(
    title    = "Distribución de Valores Faltantes por Célula (nivel proteína)",
    subtitle = "Justificación del problema de sparsity en SCP",
    x        = "Porcentaje de NAs (%)",
    y        = "Número de células"
  ) +
  theme_minimal()

print(p_nas)
ggsave("results/plots/03_diagnostico_nas_proteinas.png",
       p_nas, width = 8, height = 6, dpi = 300)

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_proteins_final.rds")
message("Script 10 finalizado. Objeto y gráfico de diagnóstico guardados.")