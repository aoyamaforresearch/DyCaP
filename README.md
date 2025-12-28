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

# -----------------------------
# 0) Read10X -> Seurat
# -----------------------------
data_dir <- "path/to/your/data" # Update this path
counts <- Read10X(data.dir = data_dir, gene.column = 1)

seu <- CreateSeuratObject(counts = counts)

# -----------------------------
# 1) Standard preprocessing (lighter ScaleData)
# -----------------------------
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000)

# Run ScaleData on variable features only (to save memory)
seu <- ScaleData(seu, features = VariableFeatures(seu)) %>%
  RunPCA(features = VariableFeatures(seu)) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# -----------------------------
# 2) Slingshot pseudotime -> store in Seurat
# -----------------------------
rd <- Embeddings(seu, "pca")[, 1:20, drop = FALSE]
cl <- seu$seurat_clusters

sds <- slingshot(rd, clusterLabels = cl)
t_raw <- slingPseudotime(sds)[, 1]

# Store pseudotime matching Seurat cell order
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

# Target number of cells per bin (approximate)
k_per_bin <- max(1, floor(n_target / nbin))

keep_idx <- unlist(lapply(idx_by_bin, function(i) {
  if (length(i) == 0) return(integer(0))
  sample(i, size = min(length(i), k_per_bin), replace = FALSE)
}))


seu_sub <- seu[, keep_idx]

# -----------------------------
# 4) Build DyCaP input dat (cells x genes + t)
# -----------------------------
keep <- !is.na(seu_sub$pseudotime01)

expr_mat <- GetAssayData(seu_sub, slot = "data")[, keep] %>%
  t() %>%
  as.matrix()

genes_use <- VariableFeatures(seu_sub) %>%
  intersect(colnames(expr_mat)) %>%
  head(100)

dat <- tibble::as_tibble(expr_mat[, genes_use, drop = FALSE]) %>%
  dplyr::mutate(t = seu_sub$pseudotime01[keep])

# Sanity checks
stopifnot(nrow(dat) == sum(keep))
stopifnot("t" %in% colnames(dat))

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

## Citation
If you use this tool, please cite our preprint:

> **Aoyama T. (2025). DyCaP: Identifying recurrent dynamic correlation patterns in single-cell RNA-seq data.** Zenodo. DOI: [10.5281/zenodo.18028361](https://doi.org/10.5281/zenodo.18028361)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
