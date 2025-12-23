library(tidyverse)
library(splines)

# ------------------------------------------------------------
# DyCaP v1.0: core functions
# ------------------------------------------------------------

gaussian_kernel <- function(x) exp(-0.5 * x^2)

local_cor_generic <- function(t_vec, x, y, tau, h) {
  w_raw <- gaussian_kernel((t_vec - tau) / h)
  sw <- sum(w_raw)
  if (is.na(sw) || sw == 0) return(NA_real_)
  w <- w_raw / sw
  
  mu_x <- sum(w * x)
  mu_y <- sum(w * y)
  
  num <- sum(w * (x - mu_x) * (y - mu_y))
  den_x <- sqrt(sum(w * (x - mu_x)^2))
  den_y <- sqrt(sum(w * (y - mu_y)^2))
  
  if (den_x == 0 || den_y == 0) return(NA_real_)
  num / (den_x * den_y)
}

dycap_fit_traj <- function(dat,
                           tau_grid = seq(0, 1, length.out = 50),
                           h_band   = 0.05,
                           K_basis  = 12,
                           genes    = NULL) {
  
  stopifnot("t" %in% colnames(dat))
  
  t_vec <- dat$t
  
  if (is.null(genes)) {
    genes <- setdiff(colnames(dat), "t")
  } else {
    genes <- intersect(genes, colnames(dat))
    genes <- setdiff(genes, "t")
  }
  
  # 安全策：NA除外
  keep <- !is.na(t_vec)
  t_vec <- t_vec[keep]
  dat2 <- dat[keep, , drop = FALSE]
  
  pair_index <- combn(genes, 2, simplify = FALSE)
  
  # 1) kernel local correlation
  cor_long <- purrr::map_dfr(
    pair_index,
    function(pair) {
      gi <- pair[1]; gj <- pair[2]
      x <- dat2[[gi]]
      y <- dat2[[gj]]
      
      rho_vec <- purrr::map_dbl(
        tau_grid,
        ~ local_cor_generic(t_vec, x, y, tau = .x, h = h_band)
      )
      
      tibble(
        gene_i = gi,
        gene_j = gj,
        tau    = tau_grid,
        rho    = rho_vec
      )
    }
  )
  
  # 2) Fisher z
  z_long <- cor_long %>%
    mutate(
      rho_clipped = pmin(pmax(rho, -0.999), 0.999),
      z = 0.5 * log((1 + rho_clipped) / (1 - rho_clipped))
    )
  
  # 3) spline basis + beta
  B_mat <- splines::bs(
    tau_grid,
    df        = K_basis,
    degree    = 3,
    intercept = TRUE
  )
  BtB_inv_Bt <- solve(t(B_mat) %*% B_mat) %*% t(B_mat)
  
  beta_tbl <- z_long %>%
    group_by(gene_i, gene_j) %>%
    summarise(z_vec = list(z), .groups = "drop") %>%
    mutate(beta = purrr::map(z_vec, ~ as.numeric(BtB_inv_Bt %*% .x)))
  
  # 4) reconstruct smooth rho(tau)
  traj_tbl <- beta_tbl %>%
    mutate(
      rho_smooth = purrr::map(
        beta,
        ~ {
          z_smooth <- as.numeric(B_mat %*% .x)
          (exp(2 * z_smooth) - 1) / (exp(2 * z_smooth) + 1)
        }
      )
    ) %>%
    dplyr::select(gene_i, gene_j, rho_smooth) %>%
    tidyr::unnest_longer(rho_smooth, values_to = "rho") %>%
    group_by(gene_i, gene_j) %>%
    mutate(tau = tau_grid) %>%
    ungroup()
  
  list(
    traj_tbl = traj_tbl,
    beta_tbl = beta_tbl,
    tau_grid = tau_grid,
    params = list(h_band = h_band, K_basis = K_basis)
  )
}

dycap_pairpair <- function(traj_tbl,
                           traj_cor_threshold = 0.99,
                           require_share_gene = FALSE) {
  
  # pair x tau wide
  rho_wide <- traj_tbl %>%
    mutate(pair = paste(gene_i, gene_j, sep = "__")) %>%
    dplyr::select(pair, tau, rho) %>%
    tidyr::pivot_wider(names_from = tau, values_from = rho)
  
  rho_mat <- rho_wide %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix()
  
  # row-scale (shape only)
  row_sd <- apply(rho_mat, 1, sd, na.rm = TRUE)
  keep_pairs <- row_sd > 0 & !is.na(row_sd)
  
  rho_mat_filt   <- rho_mat[keep_pairs, , drop = FALSE]
  rho_mat_scaled <- t(scale(t(rho_mat_filt)))
  
  # similarity among trajectories
  cor_mat <- cor(t(rho_mat_scaled), use = "pairwise.complete.obs")
  cor_mat[cor_mat >  1] <-  1
  cor_mat[cor_mat < -1] <- -1
  
  pair_ids <- rownames(cor_mat)
  idx <- which(cor_mat >= traj_cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)
  
  pairpair_tbl <- tibble(
    pair1    = pair_ids[idx[, 1]],
    pair2    = pair_ids[idx[, 2]],
    cor_traj = cor_mat[idx]
  ) %>%
    tidyr::separate(pair1, into = c("gene_i1", "gene_j1"), sep = "__", remove = FALSE) %>%
    tidyr::separate(pair2, into = c("gene_i2", "gene_j2"), sep = "__", remove = FALSE) %>%
    mutate(
      share_gene = (gene_i1 == gene_i2) |
        (gene_i1 == gene_j2) |
        (gene_j1 == gene_i2) |
        (gene_j1 == gene_j2)
    ) %>%
    arrange(desc(cor_traj))
  
  if (isTRUE(require_share_gene)) {
    pairpair_tbl <- pairpair_tbl %>% dplyr::filter(share_gene)
  }
  
  list(
    pairpair_tbl = pairpair_tbl,
    rho_mat_scaled = rho_mat_scaled,
    cor_mat = cor_mat,
    threshold = traj_cor_threshold
  )
}

# ------------------------------------------------------------
# Wrapper: DyCaP v1.0 run
# ------------------------------------------------------------

dycap_run <- function(dat,
                      tau_grid = seq(0, 1, length.out = 50),
                      h_band   = 0.05,
                      K_basis  = 12,
                      traj_cor_threshold = 0.99,
                      require_share_gene = FALSE,
                      genes = NULL) {
  
  fit <- dycap_fit_traj(
    dat      = dat,
    tau_grid = tau_grid,
    h_band   = h_band,
    K_basis  = K_basis,
    genes    = genes
  )
  
  pp <- dycap_pairpair(
    traj_tbl = fit$traj_tbl,
    traj_cor_threshold = traj_cor_threshold,
    require_share_gene = require_share_gene
  )
  
  list(
    traj_tbl     = fit$traj_tbl,
    pairpair_tbl = pp$pairpair_tbl,
    params = list(
      tau_grid = tau_grid,
      h_band = h_band,
      K_basis = K_basis,
      traj_cor_threshold = traj_cor_threshold,
      require_share_gene = require_share_gene
    )
  )
}
