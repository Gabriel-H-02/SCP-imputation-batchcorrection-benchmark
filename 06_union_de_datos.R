### ============================================================
### SCRIPT 06: UNIÓN GLOBAL DE DATOS
### ============================================================
### Consenso de proteínas por mayoría, limpieza de ceros/infinitos
### y fusión de los 134 assays individuales de péptidos en una
### única matriz global.
### OPTIMIZACIÓN DE MEMORIA: elimina los 134 assays individuales
### tras el joinAssays — son completamente redundantes a partir
### de este punto y representan el mayor consumo innecesario
### de RAM de todo el pipeline.
### INPUT:  data/leduc_peptides_aggregated.rds
### OUTPUT: data/leduc_peptides_joined.rds
### ============================================================

library(scp)
library(dplyr)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_peptides_aggregated.rds")
peptideAssays <- grep("^peptides_", names(leduc), value = TRUE)
message("Assays de péptidos individuales: ", length(peptideAssays))

# --- 2. CONSENSUS MAPPING (Voto por mayoría) ---
# Para cada péptido (modseq), asignamos la proteína más frecuente
# entre todos los lotes. Esto resuelve ambigüedades en la asignación
# péptido-proteína de forma reproducible.
rbindRowData(leduc, i = peptideAssays) %>%
  data.frame() %>%
  group_by(modseq) %>%
  mutate(Leading.razor.protein.symbol =
           names(sort(table(Leading.razor.protein), decreasing = TRUE))[1]) %>%
  select(modseq, Leading.razor.protein.symbol) %>%
  filter(!duplicated(modseq)) -> ppMap

consensus <- lapply(peptideAssays, function(i) {
  ind <- match(rowData(leduc[[i]])$modseq, ppMap$modseq)
  DataFrame(Leading.razor.protein.symbol = 
              ppMap$Leading.razor.protein.symbol[ind])
})
names(consensus) <- peptideAssays
rowData(leduc) <- consensus
message("Paso 1: Mapeo de consenso finalizado.")

# --- 3. LIMPIEZA DE CEROS E INFINITOS ---
leduc <- infIsNA(leduc, i = peptideAssays)
leduc <- zeroIsNA(leduc, i = peptideAssays)
message("Paso 2: Ceros e infinitos convertidos a NA.")

# --- 4. UNIÓN GLOBAL ---
leduc <- joinAssays(leduc, i = peptideAssays, name = "peptides")
message("Paso 3: Unión completada. Dimensiones de 'peptides': ",
        nrow(leduc[["peptides"]]), " péptidos x ",
        ncol(leduc[["peptides"]]), " células.")

# --- 5. OPTIMIZACIÓN DE MEMORIA ---

n_antes <- length(names(leduc))
assays_redundantes <- grep("^peptides_WAL", names(leduc), value = TRUE)
leduc <- leduc[, , !names(leduc) %in% assays_redundantes]
n_despues <- length(names(leduc))

message("Paso 4: Eliminados ", n_antes - n_despues,
        " assays redundantes de péptidos individuales.")
message("Assays restantes en el objeto: ", n_despues)
gc()

# --- 6. GUARDAR ---
saveRDS(leduc, file = "data/leduc_peptides_joined.rds")
message("Script 06 finalizado. Objeto guardado en 
        'data/leduc_peptides_joined.rds'")
