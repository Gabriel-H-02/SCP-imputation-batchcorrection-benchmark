### ============================================================
### SCRIPT 16_3: VOLCANO PLOT — BIOMARCADORES LIMPA LCBATCH
###              + PANEL COMPARATIVO 3 VARIANTES LIMPA
### ============================================================
### Genera el volcano plot individual de LIMPA lcbatch y un panel
### comparativo de las tres variantes de LIMPA (base, lcbatch, Set)
### para visualizar el impacto de la resolución de la corrección
### de lote en la potencia estadística del análisis diferencial.
### INPUT:  results/biomarcadores_limpa_vs_leduc.csv
###         results/biomarcadores_limpa_vs_leduc_lcbatch.csv
###         results/biomarcadores_limpa_RESOLUCION_SET.csv
### OUTPUT: results/plots/FIGURA_16_6_VOLCANO_LIMPA_LCBATCH.png
###         results/plots/FIGURA_16_7_PANEL_3_VARIANTES_LIMPA.png
###         results/tables/TOP_20_BIOMARCADORES_LIMPA_LCBATCH.csv
### ============================================================

library(EnhancedVolcano)
library(ggplot2)
library(patchwork)

if (!dir.exists("results/plots"))  dir.create("results/plots",  recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

# --- 1. CARGAR RESULTADOS ---
res_base    <- read.csv("results/biomarcadores_limpa_vs_leduc.csv",         row.names = 1)
res_lcbatch <- read.csv("results/biomarcadores_limpa_vs_leduc_lcbatch.csv", row.names = 1)
res_set     <- read.csv("results/biomarcadores_limpa_RESOLUCION_SET.csv",   row.names = 1)

# --- 2. LIMPIAR NOMBRES ---
limpiar_nombre <- function(nombre) {
  if (grepl("\\|", nombre)) gsub(".*\\|(.+)_HUMAN.*", "\\1", nombre)
  else nombre
}
etiquetas_base    <- sapply(rownames(res_base),    limpiar_nombre)
etiquetas_lcbatch <- sapply(rownames(res_lcbatch), limpiar_nombre)
etiquetas_set     <- sapply(rownames(res_set),     limpiar_nombre)

# --- 3. RESUMEN ---
n_base    <- sum(res_base$adj.P.Val    < 0.05, na.rm = TRUE)
n_lcbatch <- sum(res_lcbatch$adj.P.Val < 0.05, na.rm = TRUE)
n_set     <- sum(res_set$adj.P.Val     < 0.05, na.rm = TRUE)
message("LIMPA base    — significativas FDR < 0.05: ", n_base)
message("LIMPA lcbatch — significativas FDR < 0.05: ", n_lcbatch)
message("LIMPA Set     — significativas FDR < 0.05: ", n_set)

p_adj_cut <- 1e-10
logFC_cut <- 0.5

# --- 4. VOLCANO INDIVIDUAL LCBATCH ---
p_volcano_lcbatch <- EnhancedVolcano(
  res_lcbatch,
  lab             = etiquetas_lcbatch,
  x               = "logFC",
  y               = "adj.P.Val",
  title           = "Proteoma Diferencial: Melanoma vs Monocitos",
  subtitle        = paste0("LIMPA (corrección lcbatch) | ",
                           n_lcbatch, " proteínas significativas (FDR < 5%)"),
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

ggsave("results/plots/FIGURA_16_6_VOLCANO_LIMPA_LCBATCH.png",
       p_volcano_lcbatch, width = 10, height = 8, dpi = 300)

# --- 5. PANEL COMPARATIVO 3 VARIANTES ---
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
    labSize         = 2.5,
    col             = c("grey30", "forestgreen", "royalblue", "red2"),
    colAlpha        = 0.5,
    legendLabels    = c("NS", "Log2 FC", "p-adj", "p-adj & Log2 FC"),
    legendPosition  = "bottom",
    drawConnectors  = TRUE,
    widthConnectors = 0.4,
    max.overlaps    = 10,
    xlim            = c(-3, 3)
  ) +
    theme(plot.subtitle   = element_text(size = 10, hjust = 0.5, face = "bold"),
          legend.position = "bottom",
          legend.text     = element_text(size = 8))
}

p1 <- hacer_panel(res_base,    etiquetas_base,
                  paste0("LIMPA base (sin corrección)\n", n_base, " proteínas (FDR < 5%)"))

p2 <- hacer_panel(res_lcbatch, etiquetas_lcbatch,
                  paste0("LIMPA lcbatch (2 niveles)\n", n_lcbatch, " proteínas (FDR < 5%)"))

p3 <- hacer_panel(res_set,     etiquetas_set,
                  paste0("LIMPA Set (133 niveles)\n", n_set, " proteínas (FDR < 5%)"))

panel_3 <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Impacto de la resolución de corrección de lote en LIMPA",
    subtitle = "Base (sin corrección) → lcbatch (2 niveles) → Set (133 niveles)",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  ) &
  theme(legend.position = "bottom")

ggsave("results/plots/FIGURA_16_7_PANEL_3_VARIANTES_LIMPA.png",
       panel_3, width = 18, height = 8, dpi = 300)

# --- 6. TOP 20 LCBATCH ---
top_lcbatch              <- res_lcbatch[order(res_lcbatch$adj.P.Val), ]
top_lcbatch$nombre_corto <- sapply(rownames(top_lcbatch), limpiar_nombre)
top_lcbatch              <- top_lcbatch[, c("nombre_corto", "logFC",
                                            "AveExpr", "t", "P.Value",
                                            "adj.P.Val", "is_significant")]
write.csv(head(top_lcbatch, 20),
          "results/tables/TOP_20_BIOMARCADORES_LIMPA_LCBATCH.csv")

message("\nScript 16_3 finalizado.")
message("  results/plots/FIGURA_16_6_VOLCANO_LIMPA_LCBATCH.png")
message("  results/plots/FIGURA_16_7_PANEL_3_VARIANTES_LIMPA.png")
message("  results/tables/TOP_20_BIOMARCADORES_LIMPA_LCBATCH.csv")