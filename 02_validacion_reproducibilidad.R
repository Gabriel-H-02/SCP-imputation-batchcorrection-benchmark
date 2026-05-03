### ============================================================
### SCRIPT 02_1: VALIDACIÓN DE REPRODUCIBILIDAD VS AUTORES
### ============================================================
### Compara la matriz proteins_processed generada por el pipeline
### propio con la de Leduc et al. para demostrar que el
### preprocesamiento es fiel al original.
### INPUT:  data/leduc_final_processed.rds
###         data/leduc_author_solutions.RData
### OUTPUT: results/plots/00_validacion_reproducibilidad.png
###         results/tables/00_correlaciones_por_proteina.csv
### ============================================================

library(scp)
library(ggplot2)
library(dplyr)
library(patchwork)

# --- 1. CARGA ---
message("Cargando datos propios y soluciones de los autores...")
leduc <- readRDS("data/leduc_final_processed.rds")
load("data/leduc_author_solutions.RData")

mat_propio  <- assay(leduc[["proteins_processed"]])
mat_autores <- assay(proteins_processed_leduc)

message("Formato rownames propio:   ", rownames(mat_propio)[1])
message("Formato rownames autores:  ", rownames(mat_autores)[1])

# --- 2. NORMALIZACIÓN DE NOMBRES ---
# mat_propio:  sp|ACCESSION|NAME_HUMAN
# mat_autores: ACCESSION (solo código UniProt)
# Extraemos el accession del formato largo para poder alinear
extraer_acc <- function(nombres) {
  ifelse(grepl("\\|", nombres),
         gsub(".*\\|(.+)\\|.*", "\\1", nombres),
         nombres)
}

acc_propio  <- extraer_acc(rownames(mat_propio))
acc_autores <- rownames(mat_autores)

# --- 3. ALINEACIÓN ---
acc_comunes   <- intersect(acc_propio, acc_autores)
cells_comunes <- intersect(colnames(mat_propio), colnames(mat_autores))

message("\nProteínas comunes (tras normalizar nombres): ", length(acc_comunes))
message("Células comunes:                             ", length(cells_comunes))

if (length(acc_comunes) == 0) {
  stop("Sin proteínas comunes. Revisa el formato de los nombres.")
}

# Subconjunto alineado proteína a proteína
idx_p <- match(acc_comunes, acc_propio)
idx_a <- match(acc_comunes, acc_autores)

mat_A <- mat_propio[idx_p,  cells_comunes]
mat_B <- mat_autores[idx_a, cells_comunes]
rownames(mat_A) <- acc_comunes
rownames(mat_B) <- acc_comunes

# --- 4. CORRELACIÓN POR PROTEÍNA ---
cor_por_proteina <- sapply(seq_len(nrow(mat_A)), function(i) {
  x   <- mat_A[i, ]
  y   <- mat_B[i, ]
  idx <- !is.na(x) & !is.na(y)
  if (sum(idx) < 5) return(NA_real_)
  cor(x[idx], y[idx], method = "pearson")
})

cor_df <- data.frame(
  proteina    = acc_comunes,
  correlacion = cor_por_proteina
)

message("\n--- RESUMEN DE REPRODUCIBILIDAD ---")
message("Correlación media:     ",
        round(mean(cor_df$correlacion,   na.rm = TRUE), 4))
message("Correlación mediana:   ",
        round(median(cor_df$correlacion, na.rm = TRUE), 4))
message("Proteínas con r > 0.9: ",
        sum(cor_df$correlacion > 0.9, na.rm = TRUE), " (",
        round(100 * mean(cor_df$correlacion > 0.9, na.rm = TRUE), 1), "%)")
message("Proteínas con r < 0.5: ",
        sum(cor_df$correlacion < 0.5, na.rm = TRUE), " (",
        round(100 * mean(cor_df$correlacion < 0.5, na.rm = TRUE), 1), "%)")

# --- 5. FIGURAS ---
set.seed(42)
idx_muestra <- sample(length(acc_comunes), min(200, length(acc_comunes)))
scatter_df  <- data.frame(
  propio  = rowMeans(mat_A[idx_muestra, ], na.rm = TRUE),
  autores = rowMeans(mat_B[idx_muestra, ], na.rm = TRUE)
)

p_scatter <- ggplot(scatter_df, aes(x = autores, y = propio)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE,
              color = "darkblue", linewidth = 0.8) +
  labs(
    title    = "Reproducibilidad: Pipeline propio vs Leduc et al.",
    subtitle = paste0("Muestra de ", nrow(scatter_df),
                      " proteínas | Media de células"),
    x        = "Log2-intensidad media (Leduc et al.)",
    y        = "Log2-intensidad media (Pipeline propio)"
  ) +
  theme_minimal()

cor_media <- round(mean(cor_df$correlacion, na.rm = TRUE), 3)

p_hist <- ggplot(cor_df, aes(x = correlacion)) +
  geom_histogram(bins = 40, fill = "steelblue",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = cor_media,
             linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text", x = cor_media - 0.05, y = Inf, vjust = 2,
           label = paste0("Media = ", cor_media),
           color = "red", size = 3.5) +
  labs(
    title    = "Distribución de correlaciones por proteína",
    subtitle = "Pearson entre pipeline propio y Leduc et al. (célula a célula)",
    x        = "Correlación de Pearson (r)",
    y        = "Número de proteínas"
  ) +
  theme_minimal()

panel_repro <- p_scatter + p_hist +
  plot_annotation(
    title = "Validación de Reproducibilidad del Preprocesamiento",
    theme = theme(plot.title = element_text(size = 13, face = "bold",
                                            hjust = 0.5))
  )

print(panel_repro)

if (!dir.exists("results/plots"))  dir.create("results/plots",  recursive = TRUE)
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

ggsave("results/plots/00_validacion_reproducibilidad.png",
       panel_repro, width = 12, height = 5, dpi = 300)

# --- 6. GUARDAR TABLA ---
write.csv(cor_df[order(cor_df$correlacion), ],
          "results/tables/00_correlaciones_por_proteina.csv",
          row.names = FALSE)

message("\nScript 02_1 finalizado.")
message("  results/plots/00_validacion_reproducibilidad.png")
message("  results/tables/00_correlaciones_por_proteina.csv")