# HackBio Stage Two – Bioinformatics Data Visualization in R

## Overview

This project performs exploratory data analysis and visualization on gene expression, breast cancer, and immune response datasets using R. It is part of the [HackBio Internship](https://thehackbio.com/) Stage Two tasks.

## Datasets

| Dataset | Source | Description |
|---------|--------|-------------|
| HBR/UHR Normalized Counts | [GitHub](https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv) | Top differentially expressed genes (normalized counts) |
| HBR/UHR DEG Chr22 | [GitHub](https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv) | Differential expression results for chromosome 22 |
| Breast Cancer Wisconsin | [GitHub](https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv) | Breast cancer diagnostic features |
| HackBio Stage 2 Excel | [GitHub](https://github.com/HackBio-Internship/2025_project_collection/raw/refs/heads/main/hb_stage_2.xlsx) | Multi-sheet immune response kinetics data (panels 2a–2g) |

## Tasks & Outputs

### Part 1 – Gene Expression Analysis (`hb_stage2.R`)

| Panel | Description | Method |
|-------|-------------|--------|
| **a. Heatmap** | Visualizes normalized expression of top DEGs across HBR and UHR samples | `pheatmap` with hierarchical clustering |
| **b. Volcano Plot** | Highlights significantly up/down-regulated genes on Chr22 | `ggplot2` scatter with fold-change and p-value cutoffs |

### Part 2 – Breast Cancer Data Exploration (`hb_stage2.R`)

| Panel | Description | Method |
|-------|-------------|--------|
| **c. Scatter Plot** | Radius vs texture, colored by diagnosis (B/M) | Base R `plot()` |
| **d. Correlation Heatmap** | Pairwise correlations of 6 mean features | `pheatmap` with numeric annotations |
| **e. Scatter Plot** | Smoothness vs compactness, colored by diagnosis | Base R `plot()` |
| **f. Density Plot** | Distribution of `area_mean` for benign vs malignant | Base R `density()` + `lines()` |

### Part 3 – Figure 2 Reproduction: Global Immune Response Kinetics (`part3.R`)

Reproduces a multi-panel publication figure from immune response time-course data.

| Panel | Description | Method | Output |
|-------|-------------|--------|--------|
| **2a. Boxplot** | Cell-type ratio distributions across immune compartments | Base R `boxplot()` with custom palette | `panel_2a.pdf` |
| **2b. Scatter Plot** | Half-life vs alpha-life on log2 scale with quadrant classification | Base R `plot()` with cutoff lines | `panel_2b.pdf` |
| **2c. Heatmap** | Gene expression across cell types and time; rows clustered, columns ordered by time | `pheatmap` with column annotations (CellType, Time) | `panel_2c.pdf` |
| **2d. Heatmap** | Pathway enrichment scores across timepoints; no clustering, diverging palette centered at zero | `pheatmap` with diverging red–white–blue palette | `panel_2d.pdf` |
| **2e. Bubble Plot** | Kinetic regimes showing count-weighted half-life vs alpha by stage | `ggplot2` with `geom_point()` + size aesthetic | `panel_2e.pdf` |
| **2f. Stacked Bar** | B vs Plasma cell proportions at 0h and 72h | `ggplot2` with `geom_bar(position = "stack")` | `panel_2f.pdf` |
| **2g. Network** | Directed cell-cell interaction network from adjacency matrix | `igraph` with Fruchterman-Reingold layout | `panel_2g.png` |
| **Final Assembly** | All 7 panels combined into a single publication-quality figure | `gridExtra::grid.arrange()` with PNG grobs | `Figure_2_final.pdf` |

#### Key Conceptual Points (Part 3)

- **Log2 scaling (Panel 2b):** Compresses wide-range kinetic data; makes fold-changes symmetric
- **Cluster rows only (Panel 2c):** Time has natural order and should not be reordered; clustering genes reveals co-regulated modules
- **No clustering (Panel 2d):** Preserves curated pathway ordering and biological narrative
- **Diverging palette (Panel 2d):** Enrichment scores are directional—red (suppressed) ↔ white (neutral) ↔ blue (activated)
- **Stacked bars (Panel 2f):** Emphasizes part-to-whole composition change over time
- **Directed network (Panel 2g):** Cell signaling is directional (sender → receiver); edge weight encodes interaction strength

## Requirements

### R Packages

```r
install.packages(c("readxl", "ggplot2", "pheatmap", "dplyr", "igraph", "tidyr", "gridExtra", "png"))
```

| Package | Purpose |
|---------|---------|
| `readxl` | Reading Excel files |
| `ggplot2` | Volcano plot, bubble plot, stacked bar chart |
| `pheatmap` | Heatmaps (gene expression, pathway enrichment) |
| `dplyr` | Data manipulation |
| `igraph` | Directed network graph |
| `tidyr` | Data reshaping |
| `gridExtra` | Multi-panel figure assembly |
| `png` | Reading PNG files for final assembly |

## Project Structure

```
HackBio-StageTwo/
├── hb_stage2.R                    # Part 1 & 2: Gene expression + breast cancer analysis
├── part3.R                        # Part 3: Figure 2 reproduction (panels 2a–2g + assembly)
├── hb_stage_2.xlsx                # Downloaded data for Part 3
├── README.md                      # This file
├── gene_expression_analysis_heatmap.png
├── panel_2a.pdf
├── panel_2b.pdf
├── panel_2c.pdf
├── panel_2d.pdf
├── panel_2e.pdf
├── panel_2f.pdf
├── panel_2g.png
└── Figure_2_final.pdf             # Final multi-panel figure
```

## Author

**Blessing** — HackBio Internship, Stage Two (2025)