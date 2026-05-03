### ============================================================
### SCRIPT 12_3: ANÁLISIS DIFERENCIAL SOBRE MATRIZ KNN
### ============================================================
### Aplica el análisis de expresión diferencial (limma + eBayes)
### sobre la matriz imputada por KNN y corregida por lote
### (removeBatchEffect).
### INPUT:  data/leduc_final_processed.rds
### OUTPUT: results/biomarcadores_knn.csv
###         results/tables/TOP_20_BIOMARCADORES_KNN.csv
### ============================================================

library(scp)
library(limma)

if (!dir.exists("results"))        dir.create("results", recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables",
                                              recursive = TRUE)

# --- 1. CARGA ---
message("Cargando matriz KNN procesada...")
leduc <- readRDS("data/leduc_final_processed.rds")

assay_knn <- names(leduc)[length(names(leduc))]
message("Assays disponibles: ", paste(names(leduc), collapse = ", "))
message("Usando assay: ", assay_knn)

sce  <- getWithColData(leduc, assay_knn)
mat  <- assay(sce)
meta <- as.data.frame(colData(sce))

message("Dimensiones: ", nrow(mat), " proteínas x ", ncol(mat), " células")
message("Tipos celulares: ", paste(unique(meta$SampleType), collapse = ", "))

# --- 2. ANÁLISIS DIFERENCIAL ---
design    <- model.matrix(~ SampleType, data = meta)
fit       <- lmFit(mat, design)
fit       <- eBayes(fit)
tabla_knn <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
tabla_knn$is_significant <- tabla_knn$adj.P.Val < 0.05

n_sig <- sum(tabla_knn$is_significant, na.rm = TRUE)
message("\nProteínas significativas KNN (FDR < 0.05): ", n_sig)

# --- 3. GUARDAR ---
write.csv(tabla_knn, "results/biomarcadores_knn.csv")
top20 <- head(tabla_knn[order(tabla_knn$adj.P.Val), ], 20)
write.csv(top20, "results/tables/TOP_20_BIOMARCADORES_KNN.csv")

# --- 4. RESUMEN ---
message("\n--- RESUMEN ANÁLISIS DIFERENCIAL KNN ---")
message("Total proteínas analizadas:  ", nrow(tabla_knn))
message("Significativas (FDR < 0.05): ", n_sig)
message("Significativas (FDR < 0.01): ",
        sum(tabla_knn$adj.P.Val < 0.01, na.rm = TRUE))
message("Top biomarcador:             ", rownames(tabla_knn)[1])

message("\nScript 12_3 finalizado.")
message("  results/biomarcadores_knn.csv → alimenta el Script 15_2")