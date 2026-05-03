### ============================================================
### SCRIPT 09: REFINAMIENTO DE PÉPTIDOS Y TRANSFORMACIÓN LOG2
### ============================================================
### Filtra péptidos y células con exceso de NAs y aplica
### transformación logarítmica para estabilizar la varianza.
### INPUT:  data/leduc_peptides_normalized.rds
### OUTPUT: data/leduc_peptides_ready.rds
### ============================================================

library(scp)
library(matrixStats)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_peptides_normalized.rds")

# --- 2. FILTRADO DE NAs EXTREMOS EN FILAS (Péptidos) ---
# Eliminamos péptidos con más del 99% de valores faltantes:
# son prácticamente no informativos y solo añaden ruido.
n_antes <- nrow(leduc[["peptides_norm2"]])
leduc   <- filterNA(leduc, i = "peptides_norm2", pNA = 0.99)
n_despues <- nrow(leduc[["peptides_norm2"]])
message("Péptidos eliminados por exceso de NAs (>99%): ", n_antes - n_despues)
message("Péptidos restantes: ", n_despues)

# --- 3. FILTRADO DE NAs EXTREMOS EN COLUMNAS (Células) ---
# CORRECCIÓN: nNA() devuelve pNA en escala 0-1 (fracción), no 0-100.
# El umbral original era pNA < 99, lo que nunca filtraba nada.
# El umbral correcto es pNA < 0.99.
nnaRes   <- nNA(leduc, "peptides_norm2")
sel_cols <- nnaRes$nNAcols$pNA < 0.99   # CORRECTO: escala 0-1

n_cel_antes  <- ncol(leduc[["peptides_norm2"]])
leduc[["peptides_norm2"]] <- leduc[["peptides_norm2"]][, sel_cols]
n_cel_despues <- ncol(leduc[["peptides_norm2"]])
message("Células eliminadas por exceso de NAs (>99%): ", 
        n_cel_antes - n_cel_despues)
message("Células restantes: ", n_cel_despues)

# --- 4. TRANSFORMACIÓN LOGARÍTMICA ---
# log2 estabiliza la varianza y permite trabajar con log-fold-changes
# lineales en el análisis diferencial posterior.
leduc <- logTransform(leduc,
                      base = 2,
                      i    = "peptides_norm2",
                      name = "peptides_log")
message("Transformación Log2 completada.")

# --- 5. OPTIMIZACIÓN DE MEMORIA ---
# peptides_norm2 ya no es necesario una vez que existe peptides_log
leduc <- removeAssay(leduc, "peptides_norm2")
gc()
message("Assay intermedio 'peptides_norm2' eliminado de memoria.")

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_peptides_ready.rds")
message("Script 09 finalizado. Matriz de péptidos lista para agregación a 
        proteínas.")
