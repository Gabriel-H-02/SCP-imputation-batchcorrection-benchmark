### ============================================================
### SCRIPT 13_1: PCA Y ANÁLISIS DIFERENCIAL QRILC
### ============================================================
### Genera el PCA de la matriz QRILC procesada y realiza el
### análisis diferencial Melanoma vs Monocitos con limma.
### INPUT:  data/leduc_qrilc_processed.rds
### OUTPUT: data/leduc_qrilc_final.rds
###         results/biomarcadores_qrilc.csv
###         results/plots/05_pca_qrilc.png
### ============================================================

library(scp)
library(limma)
library(ggplot2)
library(scater)
library(matrixStats)

# --- 1. CARGA ---
leduc <- readRDS("data/leduc_qrilc_processed.rds")
sce   <- getWithColData(leduc, "proteins_qrilc_processed")
message("Dimensiones QRILC: ", nrow(sce), " proteínas x ", ncol(sce), " células")

assay_name <- assayNames(sce)[1]

# --- 2. PCA PESADO (estilo SCoPE2) ---
pcaSCoPE2 <- function(sce) {
  mat <- assay(sce)
  # Imputar NAs residuales con media de fila
  row_means <- rowMeans(mat, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx) > 0) mat[idx] <- row_means[idx[, 1]]
  cell_sim <- t(mat) %*% mat
  svd_res  <- svd(cell_sim)
  return(list(vectors = svd_res$u, values = svd_res$d))
}

pcaRes <- pcaSCoPE2(sce)
reducedDim(sce, "PCA_QRILC") <- pcaRes$vectors[, 1:5]
attr(reducedDim(sce, "PCA_QRILC"), "percentVar") <-
  (pcaRes$values / sum(pcaRes$values)) * 100

pcaPercentVar <- round(attr(reducedDim(sce, "PCA_QRILC"), "percentVar")[1:2], 1)

p_pca_qrilc <- data.frame(
  PC1        = pcaRes$vectors[, 1],
  PC2        = pcaRes$vectors[, 2],
  SampleType = colData(sce)$SampleType
) %>%
  ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = SampleType)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.5) +
  ggplot2::scale_color_manual(
    values = c("Monocyte" = "#E41A1C", "Melanoma" = "#377EB8")) +
  ggplot2::labs(
    title = "PCA — Imputación QRILC + removeBatchEffect",
    x     = paste0("PC1 (", pcaPercentVar[1], "%)"),
    y     = paste0("PC2 (", pcaPercentVar[2], "%)")
  ) +
  ggplot2::theme_minimal()

print(p_pca_qrilc)
ggplot2::ggsave("results/plots/05_pca_qrilc.png",
                p_pca_qrilc, width = 8, height = 6, dpi = 300)

message("PCA QRILC generado. Varianza PC1: ", pcaPercentVar[1], "%")

# --- 3. ANÁLISIS DIFERENCIAL ---
mat  <- assay(sce, assay_name)
meta <- as.data.frame(colData(sce))

design <- model.matrix(~ SampleType, data = meta)
fit    <- lmFit(mat, design)
fit    <- eBayes(fit)

tabla_qrilc <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
tabla_qrilc$is_significant <- tabla_qrilc$adj.P.Val < 0.05

n_sig <- sum(tabla_qrilc$is_significant, na.rm = TRUE)
message("Proteínas significativas QRILC (FDR < 0.05): ", n_sig)

write.csv(tabla_qrilc, "results/biomarcadores_qrilc.csv")

# --- 4. GUARDAR OBJETO ---
leduc <- addAssay(leduc, y = sce, name = "proteins_qrilc_with_pca")
saveRDS(leduc, "data/leduc_qrilc_final.rds")

message("Script 13_1 finalizado.")
