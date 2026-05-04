### ============================================================
### SCRIPT 16_1: VOLCANO PLOT — BIOMARCADORES LIMPA SET
### ============================================================
### Genera el volcano plot individual de LIMPA Set (máxima
### resolución, 133 placas) y un panel comparativo LIMPA base
### vs LIMPA Set para visualizar el impacto de la corrección
### de lote en la potencia estadística del análisis diferencial.
### INPUT:  results/biomarcadores_limpa_RESOLUCION_SET.csv
###         results/biomarcadores_limpa_vs_leduc.csv
### OUTPUT: results/plots/FIGURA_16_2_VOLCANO_LIMPA_SET.png
###         results/plots/FIGURA_16_3_PANEL_VOLCANOS_COMPARATIVO.png
###         results/tables/TOP_20_BIOMARCADORES_LIMPA_SET.csv
### ============================================================

library(EnhancedVolcano)
library(ggplot2)
library(patchwork)

if (!dir.exists("results/plots"))  dir.create("results/plots",  recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

# --- 1. CARGAR RESULTADOS ---
res_limpa_set  <- read.csv("results/biomarcadores_limpa_RESOLUCION_SET.csv",
                           row.names = 1)
res_limpa_base <- read.csv("results/biomarcadores_limpa_vs_leduc.csv",
                           row.names = 1)

# --- 2. LIMPIAR NOMBRES ---
limpiar_nombre <- function(nombre) {
  if (grepl("\\|", nombre)) gsub(".*\\|(.+)_HUMAN.*", "\\1", nombre)
  else nombre
}
etiquetas_set  <- sapply(rownames(res_limpa_set),  limpiar_nombre)
etiquetas_base <- sapply(rownames(res_limpa_base), limpiar_nombre)

# --- 3. RESUMEN ---
n_sig_set  <- sum(res_limpa_set$adj.P.Val  < 0.05, na.rm = TRUE)
n_sig_base <- sum(res_limpa_base$adj.P.Val < 0.05, na.rm = TRUE)
message("LIMPA base — significativas FDR < 0.05: ", n_sig_base)
message("LIMPA Set  — significativas FDR < 0.05: ", n_sig_set)

p_adj_cut <- 1e-10
logFC_cut <- 0.5

# --- 4. VOLCANO INDIVIDUAL SET ---
p_volcano_set <- EnhancedVolcano(
  res_limpa_set,
  lab             = etiquetas_set,
  x               = "logFC",
  y               = "adj.P.Val",
  title           = "Proteoma Diferencial: Melanoma vs Monocitos",
  subtitle        = paste0("LIMPA (corrección Set) | ",
                           n_sig_set, " proteínas significativas (FDR < 5%)"),
  pCutoff         = p_adj_cut,
  FCcutoff        = logFC_cut,
  pointSize       = 2.0,
  labSize         = 3.5,
  col             = c("grey30", "forestgreen", "royalblue", "red2"),
  colAlpha        = 0.5,
  legendLabels    = c("NS", "Log2 FC", "p-adj", "p-adj & Log2 FC"),
  legendPosition  = "right",
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  max.overlaps    = 20,
  xlim            = c(-3, 3)
)

print(p_volcano_set)
ggsave("results/plots/FIGURA_16_2_VOLCANO_LIMPA_SET.png",
       p_volcano_set, width = 10, height = 8, dpi = 300)

# --- 5. PANEL COMPARATIVO BASE vs SET ---
hacer_volcano_panel <- function(res, etiquetas, subtitulo) {
  EnhancedVolcano(
    res,
    lab             = etiquetas,
    x               = "logFC",
    y               = "adj.P.Val",
    title           = "",
    subtitle        = subtitulo,
    pCutoff         = p_adj_cut,
    FCcutoff        = logFC_cut,
    pointSize       = 1.8,
    labSize         = 3.0,
    col             = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha        = 0.5,
    legendLabels    = c("NS", "Log2 FC", "p-adj", "p-adj & Log2 FC"),
    legendPosition  = "bottom",
    drawConnectors  = TRUE,
    widthConnectors = 0.4,
    max.overlaps    = 15,
    xlim            = c(-3, 3)
  ) +
    theme(plot.subtitle   = element_text(size = 11, hjust = 0.5, face = "bold"),
          legend.position = "bottom",
          legend.text     = element_text(size = 9))
}

p_base_panel <- hacer_volcano_panel(
  res_limpa_base, etiquetas_base,
  paste0("LIMPA base (sin corrección de lote)\n",
         n_sig_base, " proteínas significativas (FDR < 5%)")
)

p_set_panel <- hacer_volcano_panel(
  res_limpa_set, etiquetas_set,
  paste0("LIMPA Set (corrección por 133 placas)\n",
         n_sig_set, " proteínas significativas (FDR < 5%)")
)

panel_volcanos <- (p_base_panel | p_set_panel) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Impacto de la corrección de lote en la potencia estadística (LIMPA)",
    theme = theme(plot.title = element_text(size = 14, face = "bold",
                                            hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

print(panel_volcanos)
ggsave("results/plots/FIGURA_16_3_PANEL_VOLCANOS_COMPARATIVO.png",
       panel_volcanos, width = 16, height = 9, dpi = 300)

# --- 6. TOP 20 ---
top_hits_set              <- res_limpa_set[order(res_limpa_set$adj.P.Val), ]
top_hits_set$nombre_corto <- sapply(rownames(top_hits_set), limpiar_nombre)
top_hits_set              <- top_hits_set[, c("nombre_corto", "logFC",
                                              "AveExpr", "t", "P.Value",
                                              "adj.P.Val", "is_significant")]
write.csv(head(top_hits_set, 20),
          "results/tables/TOP_20_BIOMARCADORES_LIMPA_SET.csv")

message("\nScript 16_1 finalizado.")
message("  results/plots/FIGURA_16_2_VOLCANO_LIMPA_SET.png")
message("  results/plots/FIGURA_16_3_PANEL_VOLCANOS_COMPARATIVO.png")
message("  results/tables/TOP_20_BIOMARCADORES_LIMPA_SET.csv")