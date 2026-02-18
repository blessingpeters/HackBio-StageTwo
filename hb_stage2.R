install.packages("ggplot2")
install.packages("pheatmap")
install.packages("dplyr")
library(readxl)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(igraph)
library(gridExtra)
library(grid)
library(png)

url_hbr_uhr_normalized_counts <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"

url_hbr_uhr_deg_chr22 <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"

url_breast_cancer_wisconsin <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"

hbr_uhr_normalized_counts <- read.csv(url_hbr_uhr_normalized_counts, header = TRUE)

deg_results_chr22 <- read.csv(url_hbr_uhr_deg_chr22, header = TRUE)

breast_cancer_data <- read.csv(url_breast_cancer_wisconsin, header = TRUE)

# Part 1 – Gene Expression Analysis
# a. Heatmap

colnames(hbr_uhr_normalized_counts[, 2:7])
gene_names <- hbr_uhr_normalized_counts$X

blues <- colorRampPalette(c("#FFFFFF", "#08519C"))(100)
png("panel_heatmap.png", width = 800, height = 600, res = 100)
pheatmap(
  mat = hbr_uhr_normalized_counts[, 2:7],
  border_color = "black",
  legend = T,
  labels_row = gene_names,
  fontsize_row = 6,
  color = blues,
  cluster_rows = T,
  cluster_cols = T
)
dev.off()

# b. Volcano Plot
deg_results_chr22$significance <- as.factor(deg_results_chr22$significance)

plot_volcano <- ggplot(
  deg_results_chr22,
  aes(
    x = log2FoldChange,
    y = X.log10PAdj,
    color = significance
  )
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(
    values =
      c("up" = "green", "down" = "orange", "ns" = "grey")
  ) +
  labs(
    title = "Volcano Plot: HBR vs UHR (Chr22)",
    x = "log2FoldChange",
    y = "-log10(Adjusted P-value)",
    color = "Significance"
  ) +
  theme_minimal()

ggsave("panel_volcano.png", plot_volcano, width = 8, height = 6, dpi = 100)

# ============================================================================
# ASSEMBLE: Gene Expression Analysis Grid
# ============================================================================

grob_heatmap <- rasterGrob(readPNG("panel_heatmap.png"), interpolate = TRUE)
grob_volcano <- rasterGrob(readPNG("panel_volcano.png"), interpolate = TRUE)

png("gene_expression_analysis_grid.png", width = 1600, height = 800, res = 150)
grid.arrange(
  grob_heatmap, grob_volcano,
  ncol = 2,
  top = textGrob("Part 1: Gene Expression Analysis",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()

# Part 2 – Breast Cancer Data Exploration
# c. Scatter Plot (radius vs texture)
breast_cancer_data$diagnosis <- as.factor(breast_cancer_data$diagnosis)

png("panel_scatter_radius_texture.png", width = 800, height = 600, res = 100)
plot(
  x = breast_cancer_data$radius_mean,
  y = breast_cancer_data$texture_mean,
  xlab = "radius_mean",
  ylab = "texture_mean",
  main = "Scatter Plot: Texture vs Radius",
  col = as.factor(breast_cancer_data$diagnosis),
  pch = 16,
  cex = 0.8
)
legend("topright",
  legend = c("B (Benign)", "M (Malignant)"),
  col = 1:2,
  pch = 16,
  cex = 0.8
)
dev.off()

# d. Correlation Heatmaps
features <- breast_cancer_data[, c(
  "radius_mean", "texture_mean", "perimeter_mean",
  "area_mean", "smoothness_mean", "compactness_mean"
)]

corr_mat <- cor(features, use = "complete.obs")

png("panel_corr_heatmap.png", width = 800, height = 600, res = 100)
pheatmap(
  mat = corr_mat,
  color = blues9,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Correlation Heatmap (6 Features)",
  border_color = "black"
)
dev.off()

# e. Scatter Plot: smoothness vs compactness
png("panel_scatter_smooth_compact.png", width = 800, height = 600, res = 100)
plot(
  x = breast_cancer_data$smoothness_mean,
  y = breast_cancer_data$compactness_mean,
  xlab = "smoothness_mean",
  ylab = "compactness_mean",
  main = "Scatter Plot: Compactness vs Smoothness",
  col = as.factor(breast_cancer_data$diagnosis),
  pch = 16,
  cex = 0.8
)
grid()
legend("topright",
  legend = c("B (Benign)", "M (Malignant)"),
  col = 1:2,
  pch = 16,
  cex = 0.8
)
dev.off()

# f. Density Plot (area distribution)
area_b <- breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "B"]
area_m <- breast_cancer_data$area_mean[breast_cancer_data$diagnosis == "M"]

dens_b <- density(area_b)
dens_m <- density(area_m)

png("panel_density.png", width = 800, height = 600, res = 100)
plot(dens_b,
  main = "Density Plot: area_mean (B vs M)",
  xlab = "area_mean",
  ylab = "Density",
  col = "blue",
  lwd = 2
)
lines(dens_m, col = "red", lwd = 2)
legend("topright",
  legend = c("B (Benign)", "M (Malignant)"),
  col = c("blue", "red"),
  lwd = 2,
  cex = 0.8
)
dev.off()

# ============================================================================
# ASSEMBLE: Breast Cancer Data Exploration Grid
# ============================================================================

grob_scatter1 <- rasterGrob(readPNG("panel_scatter_radius_texture.png"), interpolate = TRUE)
grob_corr     <- rasterGrob(readPNG("panel_corr_heatmap.png"), interpolate = TRUE)
grob_scatter2 <- rasterGrob(readPNG("panel_scatter_smooth_compact.png"), interpolate = TRUE)
grob_density  <- rasterGrob(readPNG("panel_density.png"), interpolate = TRUE)

png("breast_cancer_exploration_grid.png", width = 1600, height = 1400, res = 150)
grid.arrange(
  grob_scatter1, grob_corr,
  grob_scatter2, grob_density,
  ncol = 2, nrow = 2,
  top = textGrob("Part 2: Breast Cancer Data Exploration",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()
