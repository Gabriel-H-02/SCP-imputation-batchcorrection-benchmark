### ============================================================
### SCRIPT 05: NORMALIZACIÓN POR REFERENCIA Y AGREGACIÓN PSM-PÉPTIDO
### ============================================================
### Aplica dos transformaciones secuenciales:
### (1) Normalización por referencia interna: divide cada canal TMT
###     por su canal de referencia para corregir diferencias de
###     sensibilidad entre las 134 placas de adquisición.
### (2) Agregación PSM → péptido: agrupa los espectros por secuencia
###     peptídica modificada (modseq), reduciendo la redundancia y
###     acercando el dato a la unidad biológica.
### INPUT:  data/leduc_scr_filtered.rds
### OUTPUT: data/leduc_peptides_aggregated.rds
### ============================================================

library(scp)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_scr_filtered.rds")
message("Objeto cargado: ", length(names(leduc)), " assays de PSM.")

# --- 2. NORMALIZACIÓN POR REFERENCIA INTERNA ---
# Divide cada canal TMT por el canal Reference de su mismo set.
# Esto elimina variaciones multiplicativas entre placas debidas a
# diferencias en la cantidad de muestra inyectada o en la
# sensibilidad del instrumento en cada carrera de LC-MS/MS.
# samplePattern = "." aplica a TODAS las columnas (células + controles)
leduc <- divideByReference(leduc,
                           i             = names(leduc),
                           colvar        = "SampleType",
                           samplePattern = ".",
                           refPattern    = "Reference")
message("Paso 1: Normalización por referencia finalizada.")

# --- 3. FUNCIÓN DE DEDUPLICACIÓN ---
# Para péptidos detectados en múltiples PSMs dentro de un mismo set,
# seleccionamos el primer valor no-NA disponible.
# Esta estrategia conservadora es la empleada por Leduc et al.
remove.duplicates <- function(x) {
  apply(x, 2, function(xx) xx[which(!is.na(xx))[1]])
}

# --- 4. AGREGACIÓN PSM → PÉPTIDO ---
# Agrupamos por secuencia modificada (modseq), que identifica
# de forma única cada péptido incluyendo sus modificaciones
# post-traduccionales.
peptideAssays <- paste0("peptides_", names(leduc))

leduc <- aggregateFeatures(leduc,
                           i    = names(leduc),
                           fcol = "modseq",
                           name = peptideAssays,
                           fun  = remove.duplicates)

n_assays_final <- length(names(leduc))
message("Paso 2: Agregación completada.")
message("  Assays totales: ", n_assays_final,
        " (", length(peptideAssays), " PSM originales + ",
        length(peptideAssays), " péptidos agregados)")

# --- 5. GUARDAR ---
if (!dir.exists("data")) dir.create("data")
saveRDS(leduc, file = "data/leduc_peptides_aggregated.rds")
message("\nScript 05 finalizado.")
message("  data/leduc_peptides_aggregated.rds → objeto con ",
        n_assays_final, " assays")
message("  NOTA: Este archivo es el más pesado del pipeline.")
message("        Los 134 assays individuales de péptidos se eliminarán")
message("        en el Script 06 tras el joinAssays.")