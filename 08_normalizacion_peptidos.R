### ============================================================
### SCRIPT 08: NORMALIZACIÓN DE PÉPTIDOS
### ============================================================
### Normalización en dos etapas: por columna (célula) y por fila
### (péptido). Elimina el assay intermedio peptides_norm1 una vez
### que existe peptides_norm2.
### INPUT:  data/leduc_cells_filtered.rds
### OUTPUT: data/leduc_peptides_normalized.rds
### ============================================================

library(scp)
library(matrixStats)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_cells_filtered.rds")
message("Células disponibles: ", ncol(leduc[["peptides"]]))

# --- 2. NORMALIZACIÓN POR COLUMNAS (Células) ---
# Divide cada célula por su mediana global → iguala el contenido
# total de proteínas entre células y corrige diferencias de tamaño
# celular o eficiencia de carga.
mat_pep <- assay(leduc[["peptides"]])
col_meds <- colMedians(mat_pep, na.rm = TRUE)

leduc <- sweep(leduc,
               i      = "peptides",
               name   = "peptides_norm1",
               MARGIN = 2,
               FUN    = "/",
               STATS  = col_meds)
message("Paso 1: Normalización por columnas finalizada.")

# --- 3. NORMALIZACIÓN POR FILAS (Péptidos) ---
# Divide cada péptido por su mediana a través de todas las células
# → elimina sesgos sistemáticos por ionización diferencial.
leduc <- sweep(leduc,
               i      = "peptides_norm1",
               name   = "peptides_norm2",
               MARGIN = 1,
               FUN    = "/",
               STATS  = rowMedians(assay(leduc[["peptides_norm1"]]), 
                                   na.rm = TRUE))
message("Paso 2: Normalización por filas finalizada.")

# --- 4. OPTIMIZACIÓN DE MEMORIA ---
# peptides_norm1 ya no es necesario una vez que existe peptides_norm2
leduc <- removeAssay(leduc, "peptides_norm1")
gc()
message("Paso 3: Assay intermedio 'peptides_norm1' eliminado de memoria.")

# --- 5. GUARDAR ---
saveRDS(leduc, file = "data/leduc_peptides_normalized.rds")
message("Script 08 finalizado. Datos normalizados guardados.")
