### ============================================================
### SCRIPT 15_1: MÉTRICAS CUANTITATIVAS DEL BENCHMARKING
### ============================================================
### Calcula ASW, ARI, N_sig y Overlap para 6 escenarios.
### INPUT:  data/leduc_final_baseline.rds
###         data/leduc_combat_lcbatch.rds
###         data/leduc_combat_SET.rds
###         data/matriz_final_limpa_lcbatch.rds
###         data/matriz_final_limpa_SET.rds
###         data/leduc_qrilc_final.rds
###         results/biomarcadores_*.csv
### OUTPUT: results/metrics/tabla_metricas_benchmark.csv
###         results/plots/FIGURA_METRICAS_BENCHMARK.png
### ============================================================

rm(list = ls()); gc()

library(ggplot2)
library(patchwork)
library(cluster)
library(mclust)
library(SingleCellExperiment)
library(dplyr)
library(scp)

if (!dir.exists("results/metrics")) dir.create("results/metrics", recursive = TRUE)
if (!dir.exists("results/plots"))   dir.create("results/plots",   recursive = TRUE)

# ── FUNCIONES ─────────────────────────────────────────────────────────────────

calcular_ASW <- function(mat, tipos) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx) > 0) mat[idx] <- row_means[idx[, 1]]
  sil <- silhouette(as.integer(as.factor(tipos)), dist(t(mat)))
  mean(sil[, "sil_width"])
}

calcular_ARI <- function(mat, tipos, n_pc = 10) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx) > 0) mat[idx] <- row_means[idx[, 1]]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE,
                rank. = min(n_pc, ncol(mat) - 1))
  set.seed(42)
  km  <- kmeans(pca$x, centers = 2, nstart = 25)
  adjustedRandIndex(km$cluster, as.integer(as.factor(tipos)))
}

extraer_accession <- function(nombres) {
  gsub(".*\\|(.+)\\|.*", "\\1", nombres)
}

calcular_N_sig <- function(csv_path) {
  if (is.na(csv_path) || !file.exists(csv_path)) return(NA_integer_)
  res <- read.csv(csv_path, row.names = 1)
  sum(res$adj.P.Val < 0.05, na.rm = TRUE)
}

calcular_overlap <- function(csv_path, top_autores, top_n = 50) {
  if (is.na(csv_path) || !file.exists(csv_path)) return(NA_real_)
  res     <- read.csv(csv_path, row.names = 1)
  res     <- res[order(res$adj.P.Val), ]
  mis_top <- rownames(res)[seq_len(min(top_n, nrow(res)))]
  mis_top <- ifelse(grepl("\\|", mis_top), extraer_accession(mis_top), mis_top)
  round(100 * length(intersect(mis_top, top_autores)) / top_n, 1)
}

calcular_escenario <- function(mat, tipos, csv_path, top_autores, nombre) {
  message("  → ", nombre, " (", nrow(mat), " prot x ", ncol(mat), " cel)")
  asw <- tryCatch(calcular_ASW(mat, tipos),
                  error = function(e) { message("    ASW error: ", e$message); NA_real_ })
  gc()
  ari <- tryCatch(calcular_ARI(mat, tipos),
                  error = function(e) { message("    ARI error: ", e$message); NA_real_ })
  gc()
  data.frame(Escenario     = nombre,
             ASW           = round(asw, 4),
             ARI           = round(ari, 4),
             N_sig_FDR5    = calcular_N_sig(csv_path),
             Overlap_top50 = calcular_overlap(csv_path, top_autores))
}

# ── REFERENCIAS COMUNES ───────────────────────────────────────────────────────
message("Cargando referencias comunes...")
obj <- readRDS("data/leduc_final_baseline.rds")
sce <- obj[["proteins_with_pca"]]
mat_manual     <- assay(sce, assayNames(sce)[1])
tipos_baseline <- colData(sce)$SampleType
ref_rows       <- rownames(mat_manual)
ref_cols       <- colnames(mat_manual)
rm(obj, sce); gc()

load("data/leduc_author_solutions.RData")
mat_aut     <- assay(proteins_processed_leduc)
var_aut     <- apply(mat_aut, 1, var, na.rm = TRUE)
top_autores <- names(sort(var_aut, decreasing = TRUE))[1:50]
rm(mat_aut, var_aut, proteins_processed_leduc, proteins_norm_leduc); gc()

message("Referencias listas. Iniciando cálculo por escenarios...")

# ── ESCENARIO 1: KNN + removeBE ───────────────────────────────────────────────
res1 <- calcular_escenario(mat_manual, tipos_baseline,
                           "results/biomarcadores_knn.csv",
                           top_autores, "KNN + removeBE")
gc()

# ── ESCENARIO 2: ComBat lcbatch ───────────────────────────────────────────────
message("Cargando ComBat lcbatch...")
obj      <- readRDS("data/leduc_combat_lcbatch.rds")
nm       <- grep("combat", names(obj), value = TRUE, ignore.case = TRUE)[1]
sce_tmp  <- getWithColData(obj, nm)
mat_tmp  <- assay(sce_tmp)
tipos_tmp <- colData(sce_tmp)$SampleType
rm(obj, sce_tmp); gc()

res2 <- calcular_escenario(mat_tmp, tipos_tmp,
                           "results/biomarcadores_combat_lcbatch.csv",
                           top_autores, "ComBat (lcbatch)")
rm(mat_tmp, tipos_tmp); gc()

# ── ESCENARIO 3: ComBat Set ───────────────────────────────────────────────────
message("Cargando ComBat Set...")
obj      <- readRDS("data/leduc_combat_SET.rds")
nm       <- grep("combat", names(obj), value = TRUE, ignore.case = TRUE)[1]
sce_tmp  <- getWithColData(obj, nm)
mat_tmp  <- assay(sce_tmp)
tipos_tmp <- colData(sce_tmp)$SampleType
rm(obj, sce_tmp); gc()

res3 <- calcular_escenario(mat_tmp, tipos_tmp,
                           "results/biomarcadores_combat_SET.csv",
                           top_autores, "ComBat (Set)")
rm(mat_tmp, tipos_tmp); gc()

# ── ESCENARIO 4: QRILC + removeBE ────────────────────────────────────────────
message("Cargando QRILC...")
obj      <- readRDS("data/leduc_qrilc_final.rds")
nm       <- grep("processed", names(obj), value = TRUE, ignore.case = TRUE)
nm       <- if (length(nm) == 0) names(obj)[length(names(obj))] else nm[1]
message("  Assay QRILC: ", nm)
sce_tmp  <- getWithColData(obj, nm)
mat_tmp  <- assay(sce_tmp)
tipos_tmp <- colData(sce_tmp)$SampleType
rm(obj, sce_tmp); gc()

res4 <- calcular_escenario(mat_tmp, tipos_tmp,
                           "results/biomarcadores_qrilc.csv",
                           top_autores, "QRILC + removeBE")
rm(mat_tmp, tipos_tmp); gc()

# ── ESCENARIO 5: LIMPA lcbatch ────────────────────────────────────────────────
message("Cargando LIMPA lcbatch...")
mat_tmp  <- readRDS("data/matriz_final_limpa_lcbatch.rds")
mat_tmp  <- mat_tmp[intersect(rownames(mat_tmp), ref_rows),
                    intersect(colnames(mat_tmp), ref_cols)]
tipos_tmp <- tipos_baseline[seq_len(ncol(mat_tmp))]

res5 <- calcular_escenario(mat_tmp, tipos_tmp,
                           "results/biomarcadores_limpa_vs_leduc_lcbatch.csv",
                           top_autores, "LIMPA (lcbatch)")
rm(mat_tmp, tipos_tmp); gc()

# ── ESCENARIO 6: LIMPA Set ────────────────────────────────────────────────────
message("Cargando LIMPA Set...")
mat_tmp  <- readRDS("data/matriz_final_limpa_SET.rds")
mat_tmp  <- mat_tmp[intersect(rownames(mat_tmp), ref_rows),
                    intersect(colnames(mat_tmp), ref_cols)]
tipos_tmp <- tipos_baseline[seq_len(ncol(mat_tmp))]

res6 <- calcular_escenario(mat_tmp, tipos_tmp,
                           "results/biomarcadores_limpa_RESOLUCION_SET.csv",
                           top_autores, "LIMPA (Set)")
rm(mat_tmp, tipos_tmp); gc()

# ── TABLA FINAL ───────────────────────────────────────────────────────────────
tabla <- do.call(rbind, list(res1, res2, res3, res4, res5, res6))
message("\n--- TABLA DE MÉTRICAS FINALES ---")
print(tabla)
write.csv(tabla, "results/metrics/tabla_metricas_benchmark.csv", row.names = FALSE)

# ── FIGURAS ───────────────────────────────────────────────────────────────────
tabla$Escenario <- factor(tabla$Escenario, levels = tabla$Escenario)

colores <- c("KNN + removeBE"   = "#4393C3",
             "ComBat (lcbatch)" = "#D6604D",
             "ComBat (Set)"     = "#F4A582",
             "QRILC + removeBE" = "#9970AB",
             "LIMPA (lcbatch)"  = "#74C476",
             "LIMPA (Set)"      = "#238B45")

tema <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        plot.title  = element_text(face = "bold", size = 11))

p_asw <- ggplot(tabla, aes(x = Escenario, y = ASW, fill = Escenario)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = colores) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "ASW — Separación biológica",
       subtitle = "Mayor = mejor separación Melanoma/Monocito",
       y = "Average Silhouette Width", x = NULL) + tema

p_ari <- ggplot(tabla, aes(x = Escenario, y = ARI, fill = Escenario)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = colores) +
  labs(title = "ARI — Concordancia de clustering",
       subtitle = "Mayor = clustering más fiel al tipo celular real",
       y = "Adjusted Rand Index", x = NULL) + tema

p_nsig <- ggplot(tabla %>% filter(!is.na(N_sig_FDR5)),
                 aes(x = Escenario, y = N_sig_FDR5, fill = Escenario)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = colores) +
  labs(title = "N proteínas significativas (FDR < 5%)",
       subtitle = "Mayor potencia estadística del análisis diferencial",
       y = "N proteínas", x = NULL) + tema

p_overlap <- ggplot(tabla %>% filter(!is.na(Overlap_top50)),
                    aes(x = Escenario, y = Overlap_top50, fill = Escenario)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = colores) +
  labs(title = "Overlap top 50 hits vs referencia (%)",
       subtitle = "Mayor = mejor concordancia con Leduc et al.",
       y = "Overlap (%)", x = NULL) + tema

panel <- (p_asw | p_ari) / (p_nsig | p_overlap) +
  plot_annotation(
    title    = "Benchmarking Cuantitativo — 6 Escenarios",
    subtitle = "ASW y ARI: separación biológica | N_sig y Overlap: potencia estadística",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"))
  )

print(panel)
ggsave("results/plots/FIGURA_METRICAS_BENCHMARK.png",
       panel, width = 14, height = 10, dpi = 300)

message("\nScript 15_1 finalizado.")
message("  results/metrics/tabla_metricas_benchmark.csv")
message("  results/plots/FIGURA_METRICAS_BENCHMARK.png")