# DyCaP (Dynamic Correlation Patterns)

**DyCaP** is a lightweight R tool designed to identify recurrent dynamic correlation patterns in single-cell RNA-seq data. 
Unlike conventional methods that treat gene-gene correlations as static values, DyCaP models them as functions along pseudotime and detects similarities in their temporal trajectories.

## Features
- **Fast:** Runs on a standard laptop within minutes.
- **Simple:** Pure R implementation with minimal dependencies. No supercomputer required.
- **Unique:** Focuses on the "shape" of correlation changes (dynamic coherence) rather than static magnitude.

## Installation
You can use DyCaP by simply sourcing the R script provided in this repository.

1. Download `DyCaP.R` from this repository.
2. Load it in your R session.

```r
library(tidyverse)
library(splines)

# Load the core function
source("DyCaP.R")
```

## Input Data Format (Important)
DyCaP requires a data frame (or tibble) formatted as follows:

- **Rows:** Cells
- **Columns:**
    - **`t` (Required):** A numeric column representing pseudotime.
        - The range does not need to be 0-1 (it is automatically normalized internally).
    - **Other columns:** Numeric gene expression values.

**⚠️ Note on Metadata:**
By default, DyCaP treats **all columns other than `t`** as genes to be analyzed. 
If your data frame contains non-numeric metadata (e.g., `ClusterID`, `Batch`, `SampleName`), please **remove them** before running `dycap_run()`, or specify the target genes explicitly using the `genes` argument.

## Data Preparation (Seurat Users)

DyCaP is designed to work seamlessly with Seurat workflows.
Since DyCaP requires a specific input format (a data frame with a `t` column), we provide a helper script to convert your Seurat object.

<details>
<summary><b>👇 Click here to show the Seurat conversion script</b></summary>

```r
library(Seurat)
library(dplyr)
library(tibble)
library(slingshot)
library(SingleCellExperiment)

# -----------------------------
# 0) Read10X -> Seurat
# -----------------------------
data_dir <- "path/to/your/data" # Update this path
counts <- Read10X(data.dir = data_dir, gene.column = 1)
seu <- CreateSeuratObject(counts = counts)

# -----------------------------
# 1) Standard preprocessing (memory-friendly)
# -----------------------------
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000)

vf <- VariableFeatures(seu)
stopifnot(length(vf) > 1)

seu <- ScaleData(seu, features = vf) %>%
  RunPCA(features = vf) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# -----------------------------
# 2) Slingshot pseudotime -> store in Seurat (robust)
# -----------------------------
sce <- as.SingleCellExperiment(seu)
reducedDims(sce)$PCA <- Embeddings(seu, "pca")[, 1:20, drop = FALSE]
colData(sce)$cluster <- factor(seu$seurat_clusters)

sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "PCA")
t_raw <- slingPseudotime(sce)[, 1]

# Ensure name matching between slingshot output and Seurat cells
names(t_raw) <- colnames(sce)
seu$pseudotime_raw <- t_raw[colnames(seu)]

# 0-1 scaling (keep NAs as is)
rng <- range(seu$pseudotime_raw, na.rm = TRUE)
seu$pseudotime01 <- (seu$pseudotime_raw - rng[1]) / (rng[2] - rng[1])

FeaturePlot(seu, features = "pseudotime01", reduction = "umap")

# -----------------------------
# 3) Subsample cells AFTER pseudotime (recommended)
# -----------------------------
set.seed(1)

n_target <- 3000
nbin <- 40

pt <- seu$pseudotime01
idx_all <- which(!is.na(pt))

bin <- cut(pt[idx_all], breaks = nbin, labels = FALSE)
idx_by_bin <- split(idx_all, bin)

k_per_bin <- max(1, floor(n_target / nbin))

keep_idx <- unlist(lapply(idx_by_bin, function(i) {
  if (length(i) == 0) return(integer(0))
  sample(i, size = min(length(i), k_per_bin), replace = FALSE)
}))

seu_sub <- seu[, keep_idx]

# -----------------------------
# 4) Build DyCaP input dat (cells x genes + t)
# -----------------------------
pt_sub <- seu_sub$pseudotime01[colnames(seu_sub)]
cells_keep <- colnames(seu_sub)[!is.na(pt_sub)]

expr_mat <- GetAssayData(seu_sub, slot = "data")[, cells_keep, drop = FALSE] %>%
  t() %>%
  as.matrix()

genes_use <- VariableFeatures(seu_sub) %>%
  intersect(colnames(expr_mat))

# Fallback: if variable features are empty for some reason, use all genes
if (length(genes_use) < 2) {
  genes_use <- colnames(expr_mat)
}

genes_use <- head(genes_use, 100)

dat <- tibble::as_tibble(expr_mat[, genes_use, drop = FALSE]) %>%
  dplyr::mutate(t = as.numeric(pt_sub[cells_keep]))

# Sanity checks
stopifnot(nrow(dat) == length(cells_keep))
stopifnot("t" %in% colnames(dat))
stopifnot(length(genes_use) > 1)

# dat is now ready for dycap_run(dat, ...)

```
</details>

## Usage

### 1. Quick Start (Synthetic Data)
You can try DyCaP immediately using the provided example script.

```r
# 1. Load required libraries and DyCaP
library(tidyverse)
library(splines)
source("DyCaP.R")

# 2. Generate synthetic data & Run DyCaP
# (This script generates dummy data and runs the analysis automatically)
source("example_run.R")

# 3. Check the results
print(head(results$pairpair_tbl))
```

### 2. Run on Your Data
Here is a minimal example of how to run DyCaP on your own dataset.

```r
# Prepare your data
# df should have a column named "t" and numeric gene columns
# e.g., df <- data.frame(t = pseudotime, GeneA = ..., GeneB = ...)

# Run DyCaP
res <- dycap_run(
  dat = df,
  tau_grid = seq(0, 1, length.out = 50), # Resolution of time points
  h_band   = 0.05,                       # Bandwidth for kernel smoothing
  traj_cor_threshold = 0.99              # Threshold for trajectory similarity
)

# View detected recurrent patterns
print(res$pairpair_tbl)
```
### 🔬 Advanced Recipe: Filtering "Dynamic" Patterns

Sometimes, you want to find gene pairs that not only have high correlation but also show **dynamic changes** (peaks and valleys) along the pseudotime, rather than stable high correlation.

Use the following script to filter pairs based on the **number of extrema (peaks/valleys)** and visualize them.

<details>
<summary><b>👇 Click here to show the filtering & visualization script</b></summary>

```r
# ==============================================================================
# Step 0: Setup & Libraries
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# Assuming 'res' is the output from dycap_run(dat)
# res <- dycap_run(dat) 

# --- Tuning Knobs ---
EXTREMA_RHO_MIN <- 0.30  # Minimum |rho| to be considered a "strong" peak/valley
N_EXT_MIN       <- 1     # Minimum number of strong extrema required (1 = at least one peak)

# ==============================================================================
# Step 1: Helper Functions
# ==============================================================================
# Function to count "strong" extrema (peaks/valleys) to identify dynamic shapes
count_strong_extrema <- function(rho_vec, thr = 0.30) {
  rho_vec <- as.numeric(rho_vec)
  # Return NA if too short or contains NA
  if (any(!is.finite(rho_vec)) || length(rho_vec) < 5) return(NA_integer_)
  
  d1 <- diff(rho_vec)
  s1 <- sign(d1)
  
  # Handle zeros (flat regions) by filling with previous non-zero slope
  for (i in seq_along(s1)) {
    if (s1[i] == 0) {
      prev_nz <- if (i > 1) s1[rev(seq_len(i - 1))][which(s1[rev(seq_len(i - 1))] != 0)[1]] else NA
      next_nz <- if (i < length(s1)) s1[(i + 1):length(s1)][which(s1[(i + 1):length(s1)] != 0)[1]] else NA
      s1[i] <- dplyr::coalesce(prev_nz, next_nz, 0)
    }
  }
  
  # Detect sign changes (potential extrema)
  idx <- which(diff(s1) != 0) + 1
  idx <- idx[idx > 2 & idx < (length(rho_vec) - 1)] # Exclude endpoints
  
  # Count only extrema that exceed the threshold
  sum(abs(rho_vec[idx]) >= thr, na.rm = TRUE)
}

# Visualization helper
plot_dycap_pair <- function(traj_tbl, pair1, pair2, title_text) {
  pairs_show <- c(pair1, pair2)
  plot_df <- traj_tbl %>%
    dplyr::mutate(pair = paste(gene_i, gene_j, sep = "__")) %>%
    dplyr::filter(pair %in% pairs_show)
  
  ggplot(plot_df, aes(x = tau, y = rho, color = pair)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw(base_size = 14) +
    labs(title = title_text, x = "Pseudotime", y = expression(rho(tau))) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# ==============================================================================
# Step 2: Metrics Calculation & Filtering
# ==============================================================================
message("--- Calculating metrics for all pairs... ---")

# 1. Calculate metrics (max correlation, number of extrema) for each pair
pair_metrics <- res$traj_tbl %>%
  dplyr::mutate(pair = paste(gene_i, gene_j, sep = "__")) %>%
  dplyr::group_by(pair) %>%
  dplyr::summarise(
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    n_strong_extrema = count_strong_extrema(rho, thr = EXTREMA_RHO_MIN),
    .groups = "drop"
  )

# 2. Filter pairpair_tbl to keep only "interesting" dynamic patterns
pairpair_filt <- res$pairpair_tbl %>%
  dplyr::filter(share_gene) %>% # Focus on shared-gene pairs (hub analysis)
  dplyr::left_join(pair_metrics, by = c("pair1" = "pair")) %>%
  dplyr::rename(max_abs_rho1 = max_abs_rho, n_ext1 = n_strong_extrema) %>%
  dplyr::left_join(pair_metrics, by = c("pair2" = "pair")) %>%
  dplyr::rename(max_abs_rho2 = max_abs_rho, n_ext2 = n_strong_extrema) %>%
  dplyr::filter(
    max_abs_rho1 >= EXTREMA_RHO_MIN, max_abs_rho2 >= EXTREMA_RHO_MIN,
    n_ext1 >= N_EXT_MIN, n_ext2 >= N_EXT_MIN
  ) %>%
  dplyr::arrange(dplyr::desc(cor_traj))

message(paste("Filtered pairs:", nrow(pairpair_filt)))

# ==============================================================================
# Step 3: Export Significant Patterns to PDF
# ==============================================================================
if (nrow(pairpair_filt) > 0) {
  pdf_filename <- "DyCaP_Significant_Patterns.pdf"
  message(paste("Exporting plots to", pdf_filename, "..."))
  
  pdf(file = pdf_filename, width = 6, height = 4)
  
  for (i in seq_len(nrow(pairpair_filt))) {
    pp <- pairpair_filt[i, ]
    
    title_str <- paste0(
      "Rank ", i, " / ", nrow(pairpair_filt), "\n",
      pp$pair1, " vs ", pp$pair2, "\n",
      "(cor=", round(pp$cor_traj, 3), ", n_ext=", min(pp$n_ext1, pp$n_ext2), ")"
    )
    
    p <- plot_dycap_pair(res$traj_tbl, pp$pair1, pp$pair2, title_str)
    print(p)
  }
  
  dev.off()
  message("Done! Check the PDF.")
} else {
  warning("No pairs passed the filter. Try relaxing EXTREMA_RHO_MIN or N_EXT_MIN.")
}
```
</details>

## Citation
If you use this tool, please cite our preprint:

> **Aoyama T. (2025). DyCaP: Identifying recurrent dynamic correlation patterns in single-cell RNA-seq data.** Zenodo. DOI: [10.5281/zenodo.18028361](https://doi.org/10.5281/zenodo.18028361)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
