# ------------------------------------------------------------
# Example script to run DyCaP with synthetic data
# ------------------------------------------------------------

# 1. Load DyCaP functions
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("splines")) install.packages("splines")
library(tidyverse)
library(splines)

source("DyCaP.R") 

# 2. Generate Dummy Data (Synthetic) containing a "Recurrent Pattern"
set.seed(123)
n_cells <- 300
t_vec <- runif(n_cells, 0, 1)

# --- Pattern 1: Peak correlation in the middle (t=0.5) ---
# Pair 1: GeneA & GeneB
gene_A <- sin(2 * pi * t_vec) + rnorm(n_cells, 0, 0.3)
gene_B <- sin(2 * pi * t_vec) * (1 - 2*abs(t_vec - 0.5)) + rnorm(n_cells, 0, 0.3)

# Pair 2: GeneC & GeneD (Mimics Pair 1 behavior)
gene_C <- cos(2 * pi * t_vec) + rnorm(n_cells, 0, 0.3) 
gene_D <- cos(2 * pi * t_vec) * (1 - 2*abs(t_vec - 0.5)) + rnorm(n_cells, 0, 0.3)

# Noise Gene
gene_E <- rnorm(n_cells)

dat <- tibble(
  t = t_vec,
  GeneA = gene_A,
  GeneB = gene_B,
  GeneC = gene_C,
  GeneD = gene_D,
  GeneE = gene_E
)

print("Data generated. Running DyCaP...")

# 3. Run DyCaP
results <- dycap_run(
  dat = dat,
  tau_grid = seq(0, 1, length.out = 20),
  h_band = 0.1,
  K_basis = 5,
  traj_cor_threshold = 0.8, 
  genes = colnames(dat)[-1] 
)

# 4. Show results
print("--- Pair-Pair Similarity Table ---")

print(results$pairpair_tbl)

print("Done!")
