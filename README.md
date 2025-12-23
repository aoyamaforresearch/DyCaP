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

## Usage (Quick Start)
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
# The 'results' object is created by the example script
print(head(results$pairpair_tbl))
```

## Citation
If you use this tool, please cite our preprint:

> **Aoyama T. (2025). DyCaP: Identifying recurrent dynamic correlation patterns in single-cell RNA-seq data.** Zenodo. DOI: [10.5281/zenodo.18028361](https://doi.org/10.5281/zenodo.18028361)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
