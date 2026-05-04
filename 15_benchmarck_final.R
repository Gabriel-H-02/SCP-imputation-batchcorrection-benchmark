============================================================
  ### SCRIPT 15: BENCHMARKING VISUAL — 6 ESCENARIOS
  ### ============================================================
### Genera 3 figuras del benchmark comparativo:
###
### Figura 1 — Baseline individual:
###   1. Baseline (KNN + removeBE) — SCoPE2 pipeline original
###
### Figura 2 — Comparativa resolución gruesa:
###   2. ComBat lc    — KNN + ComBat (2 lotes)
###   3. LIMPA lc     — Probabilístico resolución gruesa
###   4. QRILC        — Imputación MNAR + removeBE
###
### Figura 3 — Comparativa resolución fina:
###   5. ComBat Set   — KNN + ComBat (133 placas)
###   6. LIMPA Set    — Probabilístico resolución fina
###
### INPUT:  data/leduc_final_baseline.rds
###         data/leduc_combat_lcbatch.rds
###         data/leduc_combat_SET.rds
###         data/leduc_qrilc_final.rds
###         data/matriz_final_limpa_lcbatch.rds
###         data/matriz_final_limpa_SET.rds
### OUTPUT: results/plots/FIGURA_BENCHMARK_BASELINE.png
###         results/plots/FIGURA_BENCHMARK_GRUESA.png
###         results/plots/FIGURA_BENCHMARK_FINA.png
###         results/tables/tabla_comparativa_varianza_final.csv
### ============================================================

library(ggplot2)
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(scp)

if (!dir.exists("results/plots"))  dir.create("results/plots",  recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

# --- 1. CARGA ---
leduc_baseline   <- readRDS("data/leduc_final_baseline.rds")
leduc_combat_lc  <- readRDS("data/leduc_combat_lcbatch.rds")
leduc_combat_set <- readRDS("data/leduc_combat_SET.rds")
leduc_qrilc      <- readRDS("data/leduc_qrilc_final.rds")
mat_limpa_lc     <- readRDS("data/matriz_final_limpa_lcbatch.rds")
mat_limpa_set    <- readRDS("data/matriz_final_limpa_SET.rds")

base_sce <- leduc_baseline[["proteins_with_pca"]]

# --- 2. FUNCIÓN AUXILIAR ---
generar_sce_pca <- function(matriz, base_sce, nombre_pca) {
  cols_comunes <- intersect(colnames(matriz), colnames(base_sce))
  matriz <- matriz[, cols_comunes]
  sce <- SingleCellExperiment(assays = list(reconstructed = matriz))
  colData(sce) <- colData(base_sce)[cols_comunes, ]
  sce <- runPCA(sce, name = nombre_pca,
                exprs_values = "reconstructed",
                ncomponents  = 10,
                scale        = TRUE)
  return(sce)
}

# --- 3. CALCULAR PCAs ---
message("Calculando PCAs de los 6 escenarios...")

sce_combat_lc <- getWithColData(leduc_combat_lc, "proteins_combat")
sce_combat_lc <- runPCA(sce_combat_lc, name = "PCA_COMBAT_LC",
                        exprs_values = "reconstituted", scale = TRUE)

nm_qrilc <- grep("processed", names(leduc_qrilc), value = TRUE, ignore.case = TRUE)
nm_qrilc <- if (length(nm_qrilc) == 0) names(leduc_qrilc)[length(names(leduc_qrilc))] else nm_qrilc[1]
message("Assay QRILC: ", nm_qrilc)
sce_qrilc <- getWithColData(leduc_qrilc, nm_qrilc)
sce_qrilc <- runPCA(sce_qrilc, name = "PCA_QRILC",
                    exprs_values = assayNames(sce_qrilc)[1], scale = TRUE)

sce_combat_set <- getWithColData(leduc_combat_set, "proteins_combat_SET")
sce_combat_set <- runPCA(sce_combat_set, name = "PCA_COMBAT_SET",
                         exprs_values = "reconstituted", scale = TRUE)

sce_limpa_lc  <- generar_sce_pca(mat_limpa_lc,  base_sce, "PCA_LIMPA_LC")
sce_limpa_set <- generar_sce_pca(mat_limpa_set, base_sce, "PCA_LIMPA_SET")

message("PCAs calculados correctamente.")

# --- 4. PLOTS ---
colores <- scale_color_manual(
  values = c("Monocyte" = "#1b9e77", "Melanoma" = "#d95f02")
)

tema <- theme_minimal() +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# Función para formatear los ejes PC1/PC2 con 2 decimales
fix_pc_labels <- function(p, sce, dimred) {
  pct <- attr(reducedDim(sce, dimred), "percentVar")
  if (is.null(pct)) return(p)
  p + labs(
    x = paste0(sub(" \\(.*", "", p$labels$x), " (", sprintf("%.2f", pct[1]), "%)"),
    y = paste0(sub(" \\(.*", "", p$labels$y), " (", sprintf("%.2f", pct[2]), "%)")
  )
}

# 1 — Baseline
p1 <- plotReducedDim(base_sce, dimred = "PCA_Leduc",
                     colour_by = "SampleType") +
  labs(title    = "1 — Baseline",
       subtitle = "KNN + removeBatchEffect | SCoPE2 pipeline original") +
  colores + tema
p1 <- fix_pc_labels(p1, base_sce, "PCA_Leduc")

# 2 — ComBat lc
p2 <- plotReducedDim(sce_combat_lc, dimred = "PCA_COMBAT_LC",
                     colour_by = "SampleType") +
  labs(title    = "2 — ComBat lc",
       subtitle = "KNN + ComBat (lcbatch, 2 lotes)") +
  colores + tema
p2 <- fix_pc_labels(p2, sce_combat_lc, "PCA_COMBAT_LC")

# 3 — LIMPA lc
p3 <- plotReducedDim(sce_limpa_lc, dimred = "PCA_LIMPA_LC",
                     colour_by = "SampleType") +
  labs(title    = "3 — LIMPA lc",
       subtitle = "DPC-Quant + lcbatch en modelo lineal") +
  colores + tema
p3 <- fix_pc_labels(p3, sce_limpa_lc, "PCA_LIMPA_LC")

# 4 — QRILC
p4 <- plotReducedDim(sce_qrilc, dimred = "PCA_QRILC",
                     colour_by = "SampleType") +
  labs(title    = "4 — QRILC",
       subtitle = "QRILC + removeBatchEffect") +
  colores + tema
p4 <- fix_pc_labels(p4, sce_qrilc, "PCA_QRILC")

# 5 — ComBat Set
p5 <- plotReducedDim(sce_combat_set, dimred = "PCA_COMBAT_SET",
                     colour_by = "SampleType") +
  labs(title    = "5 — ComBat Set",
       subtitle = "KNN + ComBat (Set, 133 placas)") +
  colores + tema
p5 <- fix_pc_labels(p5, sce_combat_set, "PCA_COMBAT_SET")

# 6 — LIMPA Set
p6 <- plotReducedDim(sce_limpa_set, dimred = "PCA_LIMPA_SET",
                     colour_by = "SampleType") +
  labs(title    = "6 — LIMPA Set",
       subtitle = "DPC-Quant + Set en modelo lineal") +
  colores + tema
p6 <- fix_pc_labels(p6, sce_limpa_set, "PCA_LIMPA_SET")

# --- 5. FIGURAS ---

# Figura 1: Baseline individual
ggsave("results/plots/FIGURA_BENCHMARK_BASELINE.png",
       p1, width = 7, height = 6, dpi = 300)

# Figura 2: Comparativa resolución gruesa (2, 3, 4)
panel_gruesa <- (p2 | p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Comparativa resolución gruesa (lcbatch + removeBE)",
    subtitle = "Escenarios 2 — ComBat lc | 3 — LIMPA lc | 4 — QRILC",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

ggsave("results/plots/FIGURA_BENCHMARK_GRUESA.png",
       panel_gruesa, width = 15, height = 5, dpi = 300)

# Figura 3: Comparativa resolución fina (5, 6)
panel_fina <- (p5 | p6) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Comparativa resolución fina (Set / modelo lineal)",
    subtitle = "Escenarios 5 — ComBat Set | 6 — LIMPA Set",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    )
  )

ggsave("results/plots/FIGURA_BENCHMARK_FINA.png",
       panel_fina, width = 10, height = 5, dpi = 300)

message("Figuras guardadas:")
message("  results/plots/FIGURA_BENCHMARK_BASELINE.png")
message("  results/plots/FIGURA_BENCHMARK_GRUESA.png")
message("  results/plots/FIGURA_BENCHMARK_FINA.png")

# --- 6. TABLA DE VARIANZA PC1 ---
extraer_var <- function(sce, dimred) {
  v <- attr(reducedDim(sce, dimred), "percentVar")
  if (is.null(v)) return(NA_real_)
  round(v[1], 4)
}

resumen_tfm <- data.frame(
  Escenario    = c("1 — Baseline", "2 — ComBat lc", "3 — LIMPA lc",
                   "4 — QRILC", "5 — ComBat Set", "6 — LIMPA Set"),
  Varianza_PC1 = c(
    extraer_var(base_sce,      "PCA_Leduc"),
    extraer_var(sce_combat_lc, "PCA_COMBAT_LC"),
    extraer_var(sce_limpa_lc,  "PCA_LIMPA_LC"),
    extraer_var(sce_qrilc,     "PCA_QRILC"),
    extraer_var(sce_combat_set,"PCA_COMBAT_SET"),
    extraer_var(sce_limpa_set, "PCA_LIMPA_SET")
  ),
  Metodo     = c("SCoPE2", "KNN+ComBat", "LIMPA",
                 "QRILC", "KNN+ComBat", "LIMPA"),
  Resolucion = c("Baseline", "Gruesa", "Gruesa",
                 "Gruesa", "Fina", "Fina")
)

message("\n--- VARIANZA PC1 POR ESCENARIO ---")
print(resumen_tfm)

write.csv(resumen_tfm, "results/tables/tabla_comparativa_varianza_final.csv",
          row.names = FALSE)

message("\nScript 15 finalizado.")