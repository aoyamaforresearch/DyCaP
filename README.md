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
