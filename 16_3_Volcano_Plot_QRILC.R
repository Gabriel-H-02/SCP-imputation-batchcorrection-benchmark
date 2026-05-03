### ============================================================
### SCRIPT 16_3: VOLCANO PLOT — BIOMARCADORES QRILC
### ============================================================
### Genera el volcano plot individual de QRILC y un panel
### comparativo QRILC vs LIMPA Set para visualizar la diferencia
### entre imputación fija (MNAR clásico) y modelado probabilístico
### de la DPC. Ambos métodos asumen MNAR pero difieren en cómo
### propagan la incertidumbre al análisis diferencial.
### INPUT:  results/biomarcadores_qrilc.csv
###         results/biomarcadores_limpa_RESOLUCION_SET.csv
### OUTPUT: results/plots/FIGURA_16_4_VOLCANO_QRILC.png
###         results/plots/FIGURA_16_5_PANEL_QRILC_vs_LIMPA.png
###         results/tables/TOP_20_BIOMARCADORES_QRILC.csv
### ============================================================

library(EnhancedVolcano)
library(ggplot2)
library(patchwork)

if (!dir.exists("results/plots"))  dir.create("results/plots",  recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

# --- 1. CARGAR RESULTADOS ---
res_qrilc     <- read.csv("results/biomarcadores_qrilc.csv",              row.names = 1)
res_limpa_set <- read.csv("results/biomarcadores_limpa_RESOLUCION_SET.csv", row.names = 1)

# --- 2. LIMPIAR NOMBRES ---
limpiar_nombre <- function(nombre) {
  if (grepl("\\|", nombre)) gsub(".*\\|(.+)_HUMAN.*", "\\1", nombre)
  else nombre
}
etiquetas_qrilc <- sapply(rownames(res_qrilc),     limpiar_nombre)
etiquetas_limpa <- sapply(rownames(res_limpa_set), limpiar_nombre)

# --- 3. RESUMEN ---
n_sig_qrilc <- sum(res_qrilc$adj.P.Val     < 0.05, na.rm = TRUE)
n_sig_limpa <- sum(res_limpa_set$adj.P.Val < 0.05, na.rm = TRUE)
message("QRILC     — significativas FDR < 0.05: ", n_sig_qrilc)
message("LIMPA Set — significativas FDR < 0.05: ", n_sig_limpa)

p_adj_cut <- 1e-10
logFC_cut <- 0.5

# --- 4. VOLCANO INDIVIDUAL QRILC ---
p_volcano_qrilc <- EnhancedVolcano(
  res_qrilc,
  lab             = etiquetas_qrilc,
  x               = "logFC",
  y               = "adj.P.Val",
  title           = "Proteoma Diferencial: Melanoma vs Monocitos",
  subtitle        = paste0("Imputación QRILC | ",
                           n_sig_qrilc, " proteínas significativas (FDR < 5%)"),
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

print(p_volcano_qrilc)
ggsave("results/plots/FIGURA_16_4_VOLCANO_QRILC.png",
       p_volcano_qrilc, width = 10, height = 8, dpi = 300)

# --- 5. PANEL COMPARATIVO QRILC vs LIMPA Set ---
hacer_panel <- function(res, etiquetas, subtitulo) {
  EnhancedVolcano(
    res,
    lab             = etiquetas,
    x               = "logFC",
    y               = "adj.P.Val",
    title           = "",
    subtitle        = subtitulo,
    pCutoff         = p_adj_cut,
    FCcutoff        = logFC_cut,
    pointSize       = 1.5,
    labSize         = 2.8,
    col             = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha        = 0.5,
    legendLabels    = c("NS", "Log2 FC", "p-adj", "p-adj & Log2 FC"),
    legendPosition  = "bottom",
    drawConnectors  = TRUE,
    widthConnectors = 0.4,
    max.overlaps    = 12,
    xlim            = c(-3, 3)
  ) +
    theme(plot.subtitle   = element_text(size = 10, hjust = 0.5, face = "bold"),
          legend.position = "bottom",
          legend.text     = element_text(size = 8))
}

p1 <- hacer_panel(res_qrilc,
                  etiquetas_qrilc,
                  paste0("QRILC (imputación fija MNAR)\n",
                         n_sig_qrilc, " proteínas (FDR < 5%)"))

p2 <- hacer_panel(res_limpa_set,
                  etiquetas_limpa,
                  paste0("LIMPA Set (modelado probabilístico DPC)\n",
                         n_sig_limpa, " proteínas (FDR < 5%)"))

panel_qrilc_limpa <- (p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Comparativa de Biomarcadores: QRILC vs LIMPA Set",
    subtitle = "Ambos métodos asumen MNAR — diferencia: imputación fija vs modelado probabilístico",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  ) &
  theme(legend.position = "bottom")

print(panel_qrilc_limpa)
ggsave("results/plots/FIGURA_16_5_PANEL_QRILC_vs_LIMPA.png",
       panel_qrilc_limpa, width = 14, height = 9, dpi = 300)

# --- 6. TOP 20 QRILC ---
top_qrilc              <- res_qrilc[order(res_qrilc$adj.P.Val), ]
top_qrilc$nombre_corto <- sapply(rownames(top_qrilc), limpiar_nombre)
write.csv(head(top_qrilc[, c("nombre_corto", "logFC", "AveExpr",
                             "t", "P.Value", "adj.P.Val")], 20),
          "results/tables/TOP_20_BIOMARCADORES_QRILC.csv")

message("\nScript 16_3 finalizado.")
message("  results/plots/FIGURA_16_4_VOLCANO_QRILC.png")
message("  results/plots/FIGURA_16_5_PANEL_QRILC_vs_LIMPA.png")
message("  results/tables/TOP_20_BIOMARCADORES_QRILC.csv")