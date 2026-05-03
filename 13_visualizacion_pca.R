### ============================================================
### SCRIPT 13: VISUALIZACIÓN PCA BASELINE
### ============================================================
### Genera el PCA del pipeline de referencia (Leduc et al.) en
### dos variantes:
###   (A) PCA pesado estilo SCoPE2: usa la matriz de similitud
###       entre células (t(mat) %*% mat) para reproducir el PCA
###       del paper original de Leduc et al.
###   (B) PCA estándar con scater: PCA clásico escalado para
###       comparación con los métodos de benchmarking.
### INPUT:  data/leduc_final_processed.rds
### OUTPUT: data/leduc_final_baseline.rds
###         results/plots/03_pca_weighted_leduc.png
###         results/plots/04_pca_standard.png
### ============================================================

library(scp)
library(ggplot2)
library(dplyr)
library(scater)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_final_processed.rds")
sce   <- getWithColData(leduc, "proteins_processed")

message("Dimensiones del assay de proteínas procesadas:")
print(dim(assay(sce)))
message("Nombre del assay en el SCE: ", assayNames(sce)[1])

# --- 2. PCA PESADO (Estilo SCoPE2) ---
pcaSCoPE2 <- function(sce) {
  mat       <- assay(sce)
  row_means <- rowMeans(mat, na.rm = TRUE)
  idx       <- which(is.na(mat), arr.ind = TRUE)
  mat[idx]  <- row_means[idx[, 1]]
  cell_sim  <- t(mat) %*% mat
  svd_res   <- svd(cell_sim)
  return(list(vectors = svd_res$u, values = svd_res$d))
}

pcaRes <- pcaSCoPE2(sce)

reducedDim(sce, "PCA_Leduc") <- pcaRes$vectors[, 1:5]
attr(reducedDim(sce, "PCA_Leduc"), "percentVar") <-
  (pcaRes$values / sum(pcaRes$values)) * 100

pcaPercentVar <- round(
  attr(reducedDim(sce, "PCA_Leduc"), "percentVar")[1:2], 1)

p_pca_weighted <- data.frame(
  PC1        = pcaRes$vectors[, 1],
  PC2        = pcaRes$vectors[, 2],
  SampleType = colData(sce)$SampleType
) %>%
  ggplot(aes(x = PC1, y = PC2, colour = SampleType)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Monocyte" = "#E41A1C", "Melanoma" = "#377EB8")) +
  labs(
    title = "PCA Pesado (SCoPE2): Reproducción de Leduc et al.",
    x     = paste0("PC1 (", pcaPercentVar[1], "%)"),
    y     = paste0("PC2 (", pcaPercentVar[2], "%)")
  ) +
  theme_minimal()

print(p_pca_weighted)

# --- 3. PCA ESTÁNDAR (scater) ---
# theme_minimal() corrige el fondo negro por defecto de plotPCA()
assay_name <- assayNames(sce)[1]

sce_standard <- runPCA(sce,
                       ncomponents  = 2,
                       ntop         = Inf,
                       scale        = TRUE,
                       exprs_values = assay_name)

p_pca_standard <- plotPCA(sce_standard, colour_by = "SampleType") +
  ggtitle("PCA Estándar (scater): Baseline SCoPE2") +
  scale_color_manual(
    values = c("Monocyte" = "#E41A1C", "Melanoma" = "#377EB8")) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

print(p_pca_standard)

# --- 4. GUARDAR FIGURAS ---
if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)

ggsave("results/plots/03_pca_weighted_leduc.png",
       p_pca_weighted, width = 8, height = 6, dpi = 300)
ggsave("results/plots/04_pca_standard.png",
       p_pca_standard, width = 8, height = 6, dpi = 300)

message("Figuras guardadas:")
message("  results/plots/03_pca_weighted_leduc.png")
message("  results/plots/04_pca_standard.png")

# --- 5. ACTUALIZAR OBJETO Y GUARDAR ---
leduc <- addAssay(leduc, y = sce, name = "proteins_with_pca")
saveRDS(leduc, "data/leduc_final_baseline.rds")

message("\nScript 13 finalizado. PCAs guardados y objeto baseline listo.")