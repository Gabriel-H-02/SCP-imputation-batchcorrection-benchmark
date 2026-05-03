### ============================================================
### SCRIPT 16: VOLCANO PLOT — BIOMARCADORES LIMPA
### ============================================================
### Genera el volcano plot de expresión diferencial
### Melanoma vs Monocitos usando el pipeline LIMPA base.
### MEJORA: limpieza de nombres UniProt para etiquetas legibles.
### INPUT:  results/biomarcadores_limpa_vs_leduc.csv
### OUTPUT: results/FIGURA_16_VOLCANO_LIMPA.png
###         results/tables/TOP_20_BIOMARCADORES_LIMPA.csv
### ============================================================

library(EnhancedVolcano)
library(ggplot2)

# --- 1. CARGAR RESULTADOS ---
res_limpa <- read.csv("results/biomarcadores_limpa_vs_leduc.csv", row.names = 1)

# --- 2. LIMPIAR NOMBRES DE PROTEÍNAS ---
# Los nombres vienen en formato sp|ACCESSION|NAME_HUMAN
# Extraemos solo el nombre corto (ej: LMNA, PRDX1) para las etiquetas
limpiar_nombre <- function(nombre) {
  if (grepl("\\|", nombre)) {
    # Extraer la parte después del segundo |  y antes de _HUMAN
    parte <- gsub(".*\\|(.+)_HUMAN.*", "\\1", nombre)
    return(parte)
  }
  return(nombre)
}

etiquetas <- sapply(rownames(res_limpa), limpiar_nombre)

# --- 3. RESUMEN DIAGNÓSTICO ---
message("Total proteínas en el análisis: ", nrow(res_limpa))
message("Significativas FDR < 0.05: ", sum(res_limpa$adj.P.Val < 0.05, na.rm = TRUE))
message("Rango logFC: [", round(min(res_limpa$logFC, na.rm = TRUE), 2),
        ", ", round(max(res_limpa$logFC, na.rm = TRUE), 2), "]")

# --- 4. VOLCANO PLOT ---
p_adj_cut <- 1e-10
logFC_cut <- 0.5

p_volcano <- EnhancedVolcano(res_limpa,
                             lab          = etiquetas,
                             x            = "logFC",
                             y            = "adj.P.Val",
                             title        = "Proteoma Diferencial: Melanoma vs Monocitos",
                             subtitle     = paste0("Workflow LIMPA | Total proteínas: ", nrow(res_limpa)),
                             pCutoff      = p_adj_cut,
                             FCcutoff     = logFC_cut,
                             pointSize    = 2.0,
                             labSize      = 3.5,
                             col          = c("grey30", "forestgreen", "royalblue", "red2"),
                             colAlpha     = 0.5,
                             legendLabels = c("NS", "Log2 FC", "p-adj", "p-adj & Log2 FC"),
                             legendPosition = "right",
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             max.overlaps = 20)

print(p_volcano)

# --- 5. GUARDAR FIGURA ---
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
ggsave("results/FIGURA_16_VOLCANO_LIMPA.png",
       p_volcano, width = 10, height = 8, dpi = 300)

# --- 6. TOP 20 BIOMARCADORES ---
top_hits <- res_limpa[order(res_limpa$adj.P.Val), ]
top_hits$nombre_corto <- sapply(rownames(top_hits), limpiar_nombre)
top_hits <- top_hits[, c("nombre_corto", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "is_significant")]

if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
write.csv(head(top_hits, 20),
          "results/tables/TOP_20_BIOMARCADORES_LIMPA.csv")

message("Volcano Plot generado y guardado.")
message("Top 20 biomarcadores guardados en results/tables/TOP_20_BIOMARCADORES_LIMPA.csv")
