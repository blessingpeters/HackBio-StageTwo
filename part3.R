# ============================================================================
# BIOINFORMATICS PROJECT - FIGURE 2 REPRODUCTION
# HackBio Data Visualization Internship 2026
# Part 3: Stage Two
# ============================================================================

# Task 0: Orientation and data hygiene (mandatory, graded)

install.packages("readxl")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("igraph")
install.packages("tidyr")

# Load libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(tidyr)

# Helper function for transparent colors
transparent_color <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name
  )
  invisible(t.col)
}

# Custom color palette
hb_pal <- c(
  "#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8",
  "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
  "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
  "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942"
)

# Download data
url <- "https://github.com/HackBio-Internship/2025_project_collection/raw/refs/heads/main/hb_stage_2.xlsx"
download.file(url, "hb_stage_2.xlsx", mode = "wb")


# TASK 1: PANEL 2A - Cell-type ratio distributions
# Goal: Compare distributions across immune cell types
# Requirements: Boxplot with rotated labels, relative scaling, outlier visibility

data_2a <- read_excel("hb_stage_2.xlsx", sheet = "a")

data_2a$new_ratio <- as.numeric(data_2a$new_ratio)
data_2a$cell_type <- as.factor(data_2a$cell_type)

pdf("panel_2a.pdf", width = 12, height = 12)
boxplot(new_ratio ~ cell_type,
  data = data_2a,
  ylim = c(0, 0.5),
  las = 2,
  outline = TRUE,
  pch = 1,
  col = hb_pal[seq_along(levels(data_2a$cell_type))],
  ylab = "Ratio",
  xlab = "",
  main = "Panel 2a: Cell-type Ratio Distributions"
)
dev.off()

# ============================================================================
# TASK 2: PANEL 2B - Half-life vs alpha-life scatter
# ============================================================================
# Goal: Identify kinetic regimes
# Requirements: log2 scale, cutoff lines, color-coded quadrants, labeled genes

data_2b <- read_excel("hb_stage_2.xlsx", sheet = "b")

data_2b$alpha <- as.numeric(data_2b$alpha)
data_2b$half_life <- as.numeric(data_2b$half_life)
data_2b$log2_half_life <- log2(data_2b$half_life)
data_2b$log2_alpha <- log2(data_2b$alpha)

# Define cutoffs
x_cut <- 2.5
y_cut <- -3.5

# Classify genes into groups
data_2b$group <- "grey"
data_2b$group[data_2b$log2_half_life < x_cut & data_2b$log2_alpha > y_cut] <- "#8cd17d"
data_2b$group[data_2b$log2_half_life >= x_cut] <- "#4e79a7"
data_2b$group[data_2b$log2_half_life >= x_cut & data_2b$log2_alpha > y_cut] <- "#e15759"

pdf("panel_2b.pdf", width = 12, height = 12)
plot(data_2b$log2_half_life, data_2b$log2_alpha,
  pch = 16, cex = 0.7, col = data_2b$group,
  xlab = "log2(Half Life)", ylab = "log2(Alpha Life)",
  main = "Panel 2b: Half-life vs Alpha-life Scatter"
)
abline(v = x_cut, lty = 2, lwd = 2, col = "black")
abline(h = y_cut, lty = 2, lwd = 2, col = "black")
legend("topright",
  legend = c("Ccr2", "Camp"),
  col = c("#8cd17d", "#4e79a7"), pch = 16, pt.cex = 1, bty = "n"
)
dev.off()

# CONCEPTUAL CHECKS:
# Why log2?
# - Compresses wide ranges (gene expression spans orders of magnitude)
# - Normalizes skewed distributions
# - Makes fold-changes interpretable: log2(2x) = 1, log2(4x) = 2

# What do the four quadrants mean?
# Top-left (green): Low half-life, high alpha - rapid turnover genes (Ccr2)
# Top-right (red): High half-life, high alpha - stable, highly expressed genes
# Bottom-left (grey): Low half-life, low alpha - unstable, low expression
# Bottom-right (blue): High half-life, low alpha - stable but low expression (Camp)

# ============================================================================
# TASK 3: PANEL 2C - Heatmap across cell types and time
# ============================================================================
# Goal: Visualize temporal structure across immune compartments
# Requirements: Cluster rows only, annotate CellType and Time

data_2c <- read_excel("hb_stage_2.xlsx", sheet = "c")

rownames(data_2c) <- data_2c$genes
mat_c <- as.matrix(data_2c[, -1])
mat_c <- apply(mat_c, 2, as.numeric)
rownames(mat_c) <- data_2c$genes

col_names <- colnames(mat_c)
cell_type <- sub("_[0-9]+h.*", "", col_names)
time <- sub(".*_([0-9]+h).*", "\\1", col_names)

annotation_col <- data.frame(CellType = cell_type, Time = time, row.names = col_names)

heat_cols <- colorRampPalette(c("#e15759", "white", "#4e79a7"))(100)
max_abs <- max(abs(mat_c), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(mat_c,
  color = heat_cols, breaks = breaks, border_color = "grey",
  annotation_col = annotation_col,
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_colnames = FALSE, show_rownames = FALSE,
  legend = TRUE, fontsize = 8,
  main = "Panel 2c: Heatmap Across Cell Types and Time",
  filename = "panel_2c.pdf", width = 12, height = 10
)

# KEY THINKING CHECK:
# Why cluster genes but not time?
# - Time has natural order (0h → 2h → 6h → ... → 72h)
# - Clustering would destroy this meaningful temporal sequence
# - Clustering genes reveals co-regulated modules with similar patterns
# - Preserves biological interpretation of progression over time

# ============================================================================
# TASK 4: PANEL 2D - Pathway enrichment heatmap
# ============================================================================
# Goal: Compare pathway-level responses across timepoints
# Requirements: No clustering, diverging colors centered at zero

data_2d <- read_excel("hb_stage_2.xlsx", sheet = "d_1")

data_2d <- as.data.frame(data_2d)
rownames(data_2d) <- data_2d$pathway
time_cols <- colnames(data_2d)[colnames(data_2d) != "pathway"]
mat_d <- as.matrix(data_2d[, time_cols])

max_abs <- max(abs(mat_d), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
pal <- colorRampPalette(c("red", "white", "royalblue"))(100)

pheatmap(mat_d,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = TRUE,
  border_color = "grey60", color = pal, breaks = breaks,
  legend = TRUE, fontsize_row = 9,
  main = "Panel 2d: Pathway Enrichment Heatmap",
  filename = "panel_2d.pdf", width = 8, height = 10
)

# CONCEPTUAL CHECKS:
# Why no clustering?
# - Pathways may be ordered by biological relevance or curated groups
# - Clustering could destroy meaningful organization
# - Preserves narrative structure from original analysis

# Why diverging palette?
# - Enrichment scores have direction (activation vs suppression)
# - Centered at zero distinguishes up-regulated (blue) from down-regulated (red)
# - White = no change, makes neutral values immediately visible

# ============================================================================
# TASK 5: PANEL 2E - Bubble plot of kinetic regimes
# ============================================================================
# Goal: Show count-weighted kinetic behavior

data_2e <- read_excel("hb_stage_2.xlsx", sheet = "e")

plot_2e <- ggplot(data_2e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = hb_pal) +
  scale_size_continuous(range = c(2, 15)) +
  labs(title = "Panel 2e: Kinetic Regimes", x = "Half Life", y = "Alpha") +
  theme_bw()

ggsave("panel_2e.pdf", plot_2e, width = 8, height = 6, dpi = 100)

# ============================================================================
# TASK 6: PANEL 2F - Stacked proportions
# ============================================================================
# Goal: Compare B vs Plasma cell proportions over time

data_2f <- read_excel("hb_stage_2.xlsx", sheet = "f")
data_2f <- data_2f[data_2f$stage %in% c("s00h", "s72h"), ]
data_2f$proportion <- as.numeric(data_2f$proportion)

plot_2f <- ggplot(data_2f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("B" = hb_pal[1], "Plasma" = hb_pal[2])) +
  labs(
    title = "Panel 2f: B vs Plasma Proportions",
    x = "Time", y = "Proportion", fill = "Cell Type"
  ) +
  ylim(0, 0.3) +
  theme_bw()

ggsave("panel_2f.pdf", plot_2f, width = 7, height = 6, dpi = 100)

# CONCEPTUAL CHECK:
# Why stacked instead of side-by-side?
# - Shows part-to-whole relationship (composition)
# - Easy to compare total proportion at each timepoint
# - Emphasizes how the ratio changes over time
# - Stacked format better for compositional data

# ============================================================================
# TASK 7: PANEL 2G - Directed cell-cell interaction network
# ============================================================================
# Goal: Visualize directional communication changes

data_2g <- read_excel("hb_stage_2.xlsx", sheet = "g")

# Convert to proper adjacency matrix
mat_g <- as.data.frame(data_2g)

# First column contains row names (cell types)
rownames(mat_g) <- mat_g[[1]]
mat_g <- mat_g[, -1]  # Remove first column

# Convert to numeric matrix
mat_g <- as.matrix(mat_g)
mode(mat_g) <- "numeric"

# Remove self-loops (diagonal)
diag(mat_g) <- 0

# Create graph from adjacency matrix
g <- graph_from_adjacency_matrix(mat_g,
                                 mode = "directed",
                                 weighted = TRUE,
                                 diag = FALSE)

# Remove zero-weight edges
g <- delete_edges(g, E(g)[E(g)$weight == 0])

# Set visual properties
E(g)$color <- "grey70"
w <- E(g)$weight
E(g)$width <- 0.8 + 2.0 * (w / max(w))
E(g)$arrow.size <- 0.35

V(g)$color <- hb_pal[12]
V(g)$frame.color <- "grey40"
V(g)$size <- 18
V(g)$label.color <- "blue4"
V(g)$label.cex <- 1.0
V(g)$label.font <- 2

# Create layout
set.seed(1)
layout <- layout_with_fr(g)

png("panel_2g.png", width = 800, height = 600, res = 100)
plot(g,
     layout = layout,
     edge.curved = 0,
     main = "Panel 2g: Directed Cell-Cell Interaction Network")
dev.off()


# CONCEPTUAL CHECKS:
# Why directed?
# - Cell signaling has directionality (sender → receiver)
# What does edge weight encode biologically?
# - Strength of interaction, communication frequency

# CONCEPTUAL CHECKS:
# Why directed?
# - Cell signaling has directionality (sender → receiver)
# - Not all interactions are reciprocal
# - Shows information flow through the immune system

# What does edge weight encode biologically?
# - Strength of interaction (ligand-receptor binding affinity)
# - Frequency of communication
# - Signaling intensity between cell types
# - Often derived from ligand-receptor expression levels

# ============================================================================
# TASK 8: FINAL ASSEMBLY (MANDATORY)
# ============================================================================
# Combine all panels into a single publication-quality figure

# --- RECREATE PANELS AS OBJECTS ---
install.packages("png")
library(gridExtra)
library(grid)
library(png)

grob_2a <- rasterGrob(readPNG("panel_2a.png"), interpolate = TRUE)
grob_2b <- rasterGrob(readPNG("panel_2b.png"), interpolate = TRUE)
grob_2c <- rasterGrob(readPNG("panel_2c.png"), interpolate = TRUE)
grob_2d <- rasterGrob(readPNG("panel_2d.png"), interpolate = TRUE)
grob_2e <- rasterGrob(readPNG("panel_2e.png"), interpolate = TRUE)
grob_2f <- rasterGrob(readPNG("panel_2f.png"), interpolate = TRUE)
grob_2g <- rasterGrob(readPNG("panel_2g.png"), interpolate = TRUE)

# Assemble into final figure
pdf("Figure_2_final.pdf", width = 12, height = 16)
grid.arrange(
  grob_2a, grob_2b,
  grob_2c, grob_2d,
  grob_2e, grob_2f,
  grob_2g,
  ncol = 2,
  top = textGrob("Figure 2: Global Immune Response Kinetics",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()
