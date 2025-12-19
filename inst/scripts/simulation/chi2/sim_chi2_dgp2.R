## NOTE: Replication / post-processing script.
## Expects the global results file chi2_DETAIL_GLOBAL_dgp2.csv produced by
## inst/scripts/simulations/sim_chi2_dgp2.R. Not run during R CMD check.

# ================================================================
# BOOTEI χ² (Sobol QMC) – DGP2 (PGD2)
# Simulations over n for 2x2, 3x3 and 4x4 contingency tables
#
# - Test: Pearson χ² test of independence (one-sided: greater)
#
# - Methods compared:
#     * Permutation p-value (add-one correction; B = 1, midp = FALSE)
#     * mid-p permutation p-value (B = 1, midp = TRUE)
#     * BOOTEI p-value (B > 1) with lexicographic tie-breaking
#
# - Engine: bootei() (C++/Rcpp) with test = "chisq"
#     * boot_type = "sobol"
#       (Sobol QMC resampling indices; resampling keys are then fixed and reused across permutations via coupling)
#
# - DGP (DGP2 / PGD2): RC1 model with NON-uniform monotone margins (rows and columns)
#     * 2x2: p_row = (0.25, 0.75),  q_col = (0.25, 0.75)
#     * 3x3: p_row ∝ (1, 2, 3),  q_col ∝ (1, 2, 3)
#     * 4x4: p_row ∝ (1, 2, 3, 4),  q_col ∝ (1, 2, 3, 4)
#
# - Design:
#     * Table sizes: k_row x k_col ∈ {2x2, 3x3, 4x4}
#     * H0: θ = 0  (independence; effe = FALSE)
#     * H1: θ = θ_ref(k) * n^{-gamma}, with gamma < 0.5
#           (local alternative: effect decreases with n, but effective distance increases)
#
#     * Sample size grids:
#         - 2x2: n ∈ {10, 20, 30, 50, 100, 200, 500, 1000}
#         - 3x3: n ∈ {10, 20, 30, 50}
#         - 4x4: n ∈ {10, 20, 30, 50}
#
# - Monte Carlo / computational budget:
#     * sim_MC_per_job = 250  Monte Carlo datasets per job
#     * n_rep          = 20   jobs per (k, n, scenario) cell
#     * R_perm         = 5000 permutations per dataset
#     * B_boot         = 200  BOOTEI bootstrap replicates (tie-breaker; Sobol-shifted QMC)
#
# - Parallelisation: batchtools (one job per rep_id; multicore backend)
#
# - Output (single global file; per-dataset records, H0 + H1):
#       - chi2_DETAIL_GLOBAL_pgd2.csv
#
#   NOTE: batchtools registries are created in temporary folders
#         and removed right after use; only the GLOBAL file stays in ROOT.
# ================================================================


suppressPackageStartupMessages({
  library(batchtools)
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
})

## ----------------------------------------------------------------
## 0) PATH, C++ engine and global parameters
## ----------------------------------------------------------------

ROOT <- "/media/data/corradini/prove_bootei/chi2/"
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

# C++ QMC source 
library(bootei)

set.seed(123)


alpha_test       <- 0.05
theta_H0         <- 0

sim_MC_per_job   <- 100L      # Monte Carlo simulations per job
n_rep            <- 50L     # number of jobs per cell
B_boot           <- 200L    # QMC bootstrap replicates
R_perm           <- 5000L   # permutations

# Local alternative: θ decreases with n but effective distance increases
gamma_theta      <- 0.35     # < 0.5: effective deviation grows with n

# Reference θ for each table size (can be tuned)
theta_H1_ref_2x2 <- 2.5
theta_H1_ref_3x3 <- 3
theta_H1_ref_4x4 <- 3.5

theta_H1_fun <- function(n, theta_ref, n_ref = n_ref_theta, gamma = gamma_theta) {
  theta_ref * n ** (-1*gamma)
}

# Table sizes: 2x2, 3x3, 4x4
k_pairs <- tibble::tribble(
  ~k_row, ~k_col,
  2L,     2L,
  3L,     3L,
  4L,     4L
)

# Sample-size grid for each dimension (can be adjusted)
get_n_grid <- function(k_row, k_col) {
  if (k_row == 2L && k_col == 2L) {
    as.integer(c(10, 20, 30, 50, 100, 200, 500, 1000))
  } else if (k_row == 3L && k_col == 3L) {
    as.integer(c(10, 20, 30, 50))
  } else if (k_row == 4L && k_col == 4L) {
    as.integer(c(10, 20, 30, 50))
  } else {
    stop("Sample size grid not defined for k_row = ", k_row, ", k_col = ", k_col)
  }
}

# θ_ref for each dimension
get_theta_ref <- function(k_row, k_col) {
  if (k_row == 2L && k_col == 2L) {
    theta_H1_ref_2x2
  } else if (k_row == 3L && k_col == 3L) {
    theta_H1_ref_3x3
  } else if (k_row == 4L && k_col == 4L) {
    theta_H1_ref_4x4
  } else {
    stop("theta_ref not defined for k_row = ", k_row, ", k_col = ", k_col)
  }
}

## ----------------------------------------------------------------
## 1) χ² DGP (PGD2): RC1 with NON-uniform monotone margins
## ----------------------------------------------------------------

std_scores <- function(k) as.numeric(scale(0:(k - 1), center = TRUE, scale = TRUE))

# Non-uniform monotone margins for PGD2
get_margins_pgd2 <- function(k_row, k_col) {
  
  if (k_row == 2L && k_col == 2L) {
    # 2x2 strongly unbalanced
    p <- c(0.25, 0.75)
    q <- c(0.25, 0.75)
    
  } else if (k_row == 3L && k_col == 3L) {
    # 3x3 with monotone margins
    p <- c(1, 2, 3)
    p <- p / sum(p)
    
    q <- c(1, 2, 3)
    q <- q / sum(q)
    
  } else if (k_row == 4L && k_col == 4L) {
    # 4x4 with monotone margins
    p <- c(1, 2, 3, 4)
    p <- p / sum(p)
    
    q <- c(1, 2, 3, 4)
    q <- q / sum(q)
    
  } else {
    stop("PGD2 margins not defined for k_row = ", k_row, ", k_col = ", k_col)
  }
  
  list(p = p, q = q)
}

simulate_rc1_rect_base <- function(n = 100,
                                   k_row = 3,
                                   k_col = 3,
                                   theta = 0,
                                   effe  = TRUE) {
  stopifnot(k_row >= 2L, k_col >= 2L, n >= 2L)
  cat_x <- paste0("R", seq_len(k_row))
  cat_y <- paste0("C", seq_len(k_col))
  
  # Non-uniform monotone margins (PGD2)
  marg <- get_margins_pgd2(k_row, k_col)
  p <- marg$p
  q <- marg$q
  
  u <- std_scores(k_row)
  v <- std_scores(k_col)
  
  M <- outer(p, q, "*")
  if (effe && theta != 0) {
    M <- M * exp(theta * outer(u, v, "*"))
  }
  P <- M / sum(M)
  
  P_vec <- as.vector(P)  # column-major
  idx <- sample.int(k_row * k_col, size = n, replace = TRUE, prob = P_vec)
  
  ii <- ((idx - 1L) %% k_row) + 1L
  jj <- ((idx - 1L) %/% k_row) + 1L
  
  list(
    x = factor(cat_x[ii], levels = cat_x),
    y = factor(cat_y[jj], levels = cat_y)
  )
}

# Safe wrapper: re-simulate until both margins have at least 2 levels
simulate_rc1_rect_safe <- function(..., max_retry = 100L) {
  for (i in seq_len(max_retry)) {
    dat <- simulate_rc1_rect_base(...)
    if (length(unique(dat$x)) > 1L && length(unique(dat$y)) > 1L)
      return(dat)
  }
  stop("simulate_rc1_rect_safe: exceeded max_retry without variability")
}

## ----------------------------------------------------------------
## 2) Safe wrapper around bootei() for χ² QMC
## ----------------------------------------------------------------

bootei_safe_p_chi2 <- function(x, y,
                               B,
                               R,
                               midp,
                               seed) {
  tryCatch({
    bootei(
      x           = as.character(x),
      y           = as.character(y),
      test        = "chisq",
      B           = as.integer(B),
      R           = as.integer(R),
      alternative = "greater",      # χ² ≥ 0, permutation tail is "greater"
      perm_seed   = as.numeric(seed),
      midp        = midp,
      boot_type   = 'sobol'
    )$p.value
  }, error = function(e) NA_real_)
}

## ----------------------------------------------------------------
## 3) One base cell (given n, k_row x k_col, theta)
## ----------------------------------------------------------------

run_cell_chi2_rect <- function(n,
                               k_row,
                               k_col,
                               alpha,
                               simulations,
                               theta,
                               effe,
                               B_boot,
                               R_perm,
                               seed_base = 910L,
                               max_retry_per_sim = 50L) {
  B_perm     <- 1L
  B_midp     <- 1L
  midp_perm  <- FALSE
  midp_midp  <- TRUE
  
  p_perm <- rep(NA_real_, simulations)
  p_midp <- rep(NA_real_, simulations)
  p_boot <- rep(NA_real_, simulations)
  
  i <- 1L
  while (i <= simulations) {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      if (attempts > max_retry_per_sim) {
        p_perm[i] <- NA_real_
        p_midp[i] <- NA_real_
        p_boot[i] <- NA_real_
        i <- i + 1L
        break
      }
      
      dat <- simulate_rc1_rect_safe(
        n     = n,
        k_row = k_row,
        k_col = k_col,
        theta = theta,
        effe  = effe
      )
      x <- dat$x
      y <- dat$y
      
      perm_seed_i <- as.numeric(seed_base + i * 999 + attempts)
      
      p_perm_i <- bootei_safe_p_chi2(
        x, y, B = B_perm, R = R_perm, midp = midp_perm, seed = perm_seed_i
      )
      
      p_midp_i <- bootei_safe_p_chi2(
        x, y, B = B_midp, R = R_perm, midp = midp_midp, seed = perm_seed_i
      )
      
      p_boot_i <- bootei_safe_p_chi2(
        x, y, B = B_boot, R = R_perm, midp = FALSE, seed = perm_seed_i
      )
      
      if (anyNA(c(p_perm_i, p_midp_i, p_boot_i))) next
      
      p_perm[i] <- p_perm_i
      p_midp[i] <- p_midp_i
      p_boot[i] <- p_boot_i
      i <- i + 1L
      break
    }
  }
  
  keep <- which(!is.na(p_perm) & !is.na(p_midp) & !is.na(p_boot))
  
  if (length(keep) == 0L) {
    return(tibble(
      p_perm   = numeric(0),
      p_midp   = numeric(0),
      p_boot   = numeric(0),
      sig_perm = logical(0),
      sig_midp = logical(0),
      sig_boot = logical(0)
    ))
  }
  
  tibble(
    p_perm   = p_perm[keep],
    p_midp   = p_midp[keep],
    p_boot   = p_boot[keep],
    sig_perm = p_perm[keep] < alpha,
    sig_midp = p_midp[keep] < alpha,
    sig_boot = p_boot[keep] < alpha
  )
}

## ----------------------------------------------------------------
## 4) Job function for batchtools
## ----------------------------------------------------------------

chi2_rep_task <- function(rep_id,
                          n,
                          k_row,
                          k_col,
                          alpha,
                          theta,
                          effe,
                          sim_MC_per_job,
                          B_boot,
                          R_perm) {
  Sys.setenv(
    OMP_NUM_THREADS        = "1",
    OPENBLAS_NUM_THREADS   = "1",
    MKL_NUM_THREADS        = "1",
    BLIS_NUM_THREADS       = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS    = "1",
    GOTO_NUM_THREADS       = "1"
  )
  set.seed(100 + as.integer(rep_id))
  
  detail <- run_cell_chi2_rect(
    n           = n,
    k_row       = k_row,
    k_col       = k_col,
    alpha       = alpha,
    simulations = sim_MC_per_job,
    theta       = theta,
    effe        = effe,
    B_boot      = B_boot,
    R_perm      = R_perm,
    seed_base   = as.numeric(123 + as.integer(rep_id) * 234 +
                               n + 345 * k_row + 456 * k_col)
  )
  
  scenario_lab <- if (abs(theta) < 1e-12 || !effe) "H0 (Type-I)" else "H1 (Power)"
  
  detail %>%
    mutate(
      n        = n,
      k_row    = k_row,
      k_col    = k_col,
      alpha    = alpha,
      theta    = theta,
      rep      = rep_id,
      scenario = scenario_lab
    )
}

## ----------------------------------------------------------------
## 5) Batch execution for ONE cell (k_row, k_col, n, theta)
##    (registry in tempdir(), nothing per-cell saved in ROOT)
## ----------------------------------------------------------------

run_one_n_batch <- function(k_row,
                            k_col,
                            n,
                            alpha,
                            theta,
                            effe,
                            sim_MC_per_job,
                            n_rep,
                            B_boot,
                            R_perm,
                            root_dir = ROOT) {
  
  scen_lab <- if (abs(theta) < 1e-12 || !effe) "H0" else "H1"
  
  # Registry in a TEMPORARY folder, not in ROOT
  regdir <- tempfile(
    pattern = sprintf("chi2_k%dx%d_n%03d_%s_pgd2_", k_row, k_col, n, scen_lab),
    tmpdir  = tempdir()
  )
  
  N_CORES <- min(parallel::detectCores(), n_rep)
  
  reg <- makeRegistry(
    file.dir  = regdir,
    packages  = c("Rcpp", "dplyr", "tibble", "tidyr", "readr", "purrr")
  )
  reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = N_CORES)
  
  batchMap(
    fun      = chi2_rep_task,
    rep_id   = 1:n_rep,
    more.args = list(
      n              = n,
      k_row          = k_row,
      k_col          = k_col,
      alpha          = alpha,
      theta          = theta,
      effe           = effe,
      sim_MC_per_job = sim_MC_per_job,
      B_boot         = B_boot,
      R_perm         = R_perm
    ),
    reg = reg
  )
  
  ids <- batchtools::findJobs(reg = reg)
  batchtools::submitJobs(ids = ids, resources = list(ncpus = 1L), reg = reg)
  cat(" Running (PGD2):", regdir, "on", N_CORES, "cores...\n")
  
  ok <- batchtools::waitForJobs(ids = ids, reg = reg)
  cat("waitForJobs result:", ok, "for", regdir, "\n")
  
  if (!ok) {
    cat(" Failed jobs in", regdir, ":\n")
    print(batchtools::getErrorMessages(reg = reg))
    # Try to remove registry anyway
    try(batchtools::removeRegistry(reg = reg), silent = TRUE)
    stop("Stopping: fix job errors before aggregating.")
  }
  
  done_ids <- batchtools::findDone(reg = reg)
  res_list <- batchtools::reduceResultsList(ids = done_ids, reg = reg)
  
  # Remove registry folder (and everything inside)
  try(batchtools::removeRegistry(reg = reg), silent = TRUE)
  
  if (length(res_list) == 0L) {
    warning("No results in ", regdir)
    return(tibble())
  }
  
  detail_all <- dplyr::bind_rows(res_list)
  
  cat(" Completed cell (PGD2): k=",
      k_row, "x", k_col,
      ", n=", n,
      " — detail=", nrow(detail_all), "\n",
      sep = "")
  
  detail_all
}

## ----------------------------------------------------------------
## 6) MAIN LOOP over (k_row, k_col) and n (H0 + H1)
## ----------------------------------------------------------------

detail_global_list  <- list()

for (idx in seq_len(nrow(k_pairs))) {
  k_row_i <- k_pairs$k_row[idx]
  k_col_i <- k_pairs$k_col[idx]
  
  n_grid      <- get_n_grid(k_row_i, k_col_i)
  theta_ref_i <- get_theta_ref(k_row_i, k_col_i)
  
  message("\n=======================")
  message(sprintf(">>> χ² (PGD2): k = %dx%d", k_row_i, k_col_i))
  message("=======================\n")
  
  for (n_cur in n_grid) {
    message(sprintf(">>> [H1] k=%dx%d, n=%d", k_row_i, k_col_i, n_cur))
    
    # H1: theta > 0, "smoothed" effect that decreases with n
    theta_H1_cur <- theta_H1_fun(n_cur, theta_ref = theta_ref_i)
    
    res_H1 <- run_one_n_batch(
      k_row          = k_row_i,
      k_col          = k_col_i,
      n              = n_cur,
      alpha          = alpha_test,
      theta          = theta_H1_cur,
      effe           = TRUE,
      sim_MC_per_job = sim_MC_per_job,
      n_rep          = n_rep,
      B_boot         = B_boot,
      R_perm         = R_perm,
      root_dir       = ROOT
    )
    
    message(sprintf(">>> [H0] k=%dx%d, n=%d", k_row_i, k_col_i, n_cur))
    
    # H0: theta = 0, no RC1 effect
    res_H0 <- run_one_n_batch(
      k_row          = k_row_i,
      k_col          = k_col_i,
      n              = n_cur,
      alpha          = alpha_test,
      theta          = theta_H0,
      effe           = FALSE,
      sim_MC_per_job = sim_MC_per_job,
      n_rep          = n_rep,
      B_boot         = B_boot,
      R_perm         = R_perm,
      root_dir       = ROOT
    )
    
    detail_global_list <- c(detail_global_list, list(res_H1, res_H0))
  }
}

cat("\n All combinations (k_row, k_col, n) completed (PGD2).\n")

## -----------------------------------------------------------------
## 7) Aggregate everything into a single GLOBAL DETAIL (PGD2)
## -----------------------------------------------------------------

detail_global  <- dplyr::bind_rows(detail_global_list)

if (nrow(detail_global) > 0L) {
  write_csv(detail_global, file.path(ROOT, "chi2_DETAIL_GLOBAL_dgp2.csv"))
}

cat("    - Detail global rows (PGD2):", nrow(detail_global),  "\n")
