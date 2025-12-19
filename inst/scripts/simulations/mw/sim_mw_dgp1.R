## NOTE: Replication / post-processing script.
## Expects the global results file mw_DETAIL_GLOBAL_dgp1.csv produced by
## inst/scripts/simulations/sim_mw_dgp1.R Not run during R CMD check.

# ================================================================
# BOOTEI Mann–Whitney (Sobol QMC) – DGP1
# Simulations over (n1, n2) for ordinal outcomes with k = 3, 5, 7 levels
#
# - Test: Mann–Whitney (Wilcoxon rank-sum) test
# - Methods compared:
#     * Permutation p-value (add-one correction; B = 1, midp = FALSE)
#     * mid-p permutation p-value (B = 1, midp = TRUE)
#     * BOOTEI p-value (B > 1) with lexicographic tie-breaking
#
# - Engine: bootei() (C++/Rcpp) with test = "mannwhitney"
#     * boot_type = "sobol"
#       (Sobol QMC resampling indices, resampling keys are then fixed and coupled across permutations)
#
# - DGP (latent-Likert):
#     * Two independent groups of sizes n1 and n2
#     * Latent variables:
#         Z1 ~ N(0, 1),  Z2 ~ N(assoc, 1)
#       discretised via common Normal-quantile cutpoints into k ordered categories
#     * H0: assoc = 0  (no location shift; effe = FALSE)
#     * H1: assoc = assoc_ref(k) * (n1 + n2)^(-gamma), with gamma < 0.5
#           (local alternative: effect decreases with n1+n2, but effective distance increases)
#
# - Design:
#     * Likert levels: k ∈ {3, 5, 7}
#     * Sample-size grids are k-specific (balanced designs: n1 = n2):
#         - k = 3: (n1,n2) ∈ {(5,5), (10,10), (25,25), (50,50), (100,100), (250,250)}
#         - k = 5: (n1,n2) ∈ {(5,5), (10,10), (15,15), (25,25), (50,50)}
#         - k = 7: (n1,n2) ∈ {(5,5), (10,10), (15,15), (25,25), (50,50)}
#
# - Monte Carlo / computational budget:
#     * sim_MC_per_job = 250  Monte Carlo datasets per job
#     * n_rep          = 20   jobs per (k, n1, n2, scenario) cell
#     * R_perm         = 5000 permutations per dataset
#     * B_boot         = 200  BOOTEI bootstrap replicates (tie-breaker; Sobol-shifted QMC)
#
# - Parallelisation: batchtools (one job per rep_id; multicore backend)
#
# - Output (single global file; per-dataset records, H0 + H1):
#       - mw_DETAIL_GLOBAL_dgp1.csv
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
## 0) PATH, C++ and GLOBAL PARAMETERS
## ----------------------------------------------------------------

ROOT <- "."
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

# install.packages("remotes")
# remotes::install_github("corradinigiovanni/bootei")  
library(bootei)

set.seed(123)

alpha_test     <- 0.05
assoc_H0       <- 0           # H0: no effect

sim_MC_per_job   <- 100L      # Monte Carlo simulations per job
n_rep            <- 50L      # number of jobs per cell
B_boot           <- 200L     # QMC bootstrap replicates
R_perm           <- 5000L    # permutations

# Local alternative "smoothed":
gamma_assoc    <- 0.35

# Reference effects for k = 3, 5, 7
assoc_H1_ref_3 <- 7/3
assoc_H1_ref_5 <- 8/3
assoc_H1_ref_7 <- 9/3

# List of k-values (Likert levels)
k_list <- c(3L, 5L, 7L)

# ----------------------------------------------------------------
# NEW: sample-size grid specific per k
# ----------------------------------------------------------------

get_mw_pairs <- function(k) {
  if (k == 3L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      25L,  25L,
      50L,  50L,
      100L, 100L,
      250L, 250L
    )
  } else if (k == 5L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      15L,  15L,
      25L,  25L,
      50L,  50L
    )
  } else if (k == 7L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      15L,  15L,
      25L,  25L,
      50L,  50L
    )
  } else {
    stop("sample-size grid not defined for k = ", k)
  }
}

# assoc_ref for each k
get_assoc_ref <- function(k) {
  if (k == 3L) {
    assoc_H1_ref_3
  } else if (k == 5L) {
    assoc_H1_ref_5
  } else if (k == 7L) {
    assoc_H1_ref_7
  } else {
    stop("assoc_ref not defined for k = ", k)
  }
}

# Alternative H1: depends only on n1+n2 and assoc_ref_k
assoc_H1_fun <- function(n1, n2,
                         assoc_ref,
                         gamma = gamma_assoc) {
  assoc_ref * (n1 + n2) ** (-1*gamma)
}

## ----------------------------------------------------------------
## 1) Latent Likert DGP and BOOTEI wrapper for Mann–Whitney
## ----------------------------------------------------------------

# Utility: pin BLAS/OMP threads to 1
pin_threads <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    GOTO_NUM_THREADS = "1"
  )
}

# DGP: Z ~ N(mu, 1), same normal-quantile cutpoints for both groups
simulate_group_likert_latent <- function(n, k = 5, mu = 0) {
  z <- rnorm(n, mean = mu, sd = 1)
  cuts <- qnorm(seq(0, 1, length.out = k + 1))
  as.numeric(cut(z, breaks = cuts, labels = FALSE, include.lowest = TRUE))
}

# MW generator: latent shift on Y of magnitude "assoc"
simulate_mw <- function(n1, n2, k, assoc = 0, effe = TRUE) {
  if (!effe || assoc == 0) {
    mu1 <- 0
    mu2 <- 0
  } else {
    mu1 <- 0
    mu2 <- assoc   # shift on Y only
  }
  g1 <- simulate_group_likert_latent(n1, k, mu = mu1)
  g2 <- simulate_group_likert_latent(n2, k, mu = mu2)
  list(x = g1, y = g2)
}

# Safe wrapper around bootei() for MW
bootei_safe_p_mw <- function(x, y,
                             B,
                             R,
                             midp,
                             seed) {
  tryCatch({
    bootei(
      x           = as.numeric(x),
      y           = as.numeric(y),
      test        = "mannwhitney",
      B           = as.integer(B),
      R           = as.integer(R),
      alternative = "two.sided",
      perm_seed   = as.numeric(seed),
      midp        = midp,
      boot_type   = 'sobol'
    )$p.value
  }, error = function(e) NA_real_)
}

## ----------------------------------------------------------------
## 2) One MW cell (one (n1,n2), one k, one assoc)
## ----------------------------------------------------------------

run_cell_mw <- function(n1,
                        n2,
                        k,
                        alpha,
                        assoc,
                        effe,
                        B_boot,
                        R_perm,
                        simulations = sim_MC_per_job,
                        seed_base = 910L,
                        max_retry_per_sim = 50L) {
  B_perm    <- 1L
  B_midp    <- 1L
  midp_perm <- FALSE
  midp_midp <- TRUE
  
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
      
      dat <- simulate_mw(n1 = n1, n2 = n2, k = k, assoc = assoc, effe = effe)
      x <- as.numeric(dat$x)
      y <- as.numeric(dat$y)
      
      seed_i <- as.numeric(seed_base + i * 999 + attempts)
      
      p_perm_i <- bootei_safe_p_mw(
        x, y, B = B_perm, R = R_perm, midp = midp_perm, seed = seed_i
      )
      p_midp_i <- bootei_safe_p_mw(
        x, y, B = B_midp, R = R_perm, midp = midp_midp, seed = seed_i
      )
      p_boot_i <- bootei_safe_p_mw(
        x, y, B = B_boot, R = R_perm, midp = FALSE, seed = seed_i
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
## 3) Job function for batchtools
## ----------------------------------------------------------------

mw_rep_task <- function(rep_id,
                        n1,
                        n2,
                        k,
                        alpha,
                        assoc,
                        effe,
                        sim_MC_per_job,
                        B_boot,
                        R_perm) {
  pin_threads()
  set.seed(100 + as.integer(rep_id))
  
  detail <- run_cell_mw(
    n1          = n1,
    n2          = n2,
    k           = k,
    alpha       = alpha,
    assoc       = assoc,
    effe        = effe,
    B_boot      = B_boot,
    R_perm      = R_perm,
    simulations = sim_MC_per_job,
    seed_base   = as.numeric(123 + as.integer(rep_id) * 234 +
                               n1 + 345 * n2 + 456 * k)
  )
  
  scenario_lab <- if (assoc == 0 || !effe) "H0 (Type-I)" else "H1 (Power)"
  
  detail %>%
    mutate(
      n1       = n1,
      n2       = n2,
      k        = k,
      alpha    = alpha,
      assoc    = assoc,
      rep      = rep_id,
      scenario = scenario_lab
    )
}

## ----------------------------------------------------------------
## 4) Batch execution for ONE cell (k, n1, n2, assoc)
##    (registry in tempdir(), no per-cell files in ROOT)
## ----------------------------------------------------------------

run_one_pair_batch_mw <- function(k,
                                  n1,
                                  n2,
                                  alpha,
                                  assoc,
                                  effe,
                                  sim_MC_per_job,
                                  n_rep,
                                  B_boot,
                                  R_perm,
                                  root_dir = ROOT) {
  
  scen_lab <- if (assoc == 0 || !effe) "H0" else "H1"
  
  regdir <- tempfile(
    pattern = sprintf("mw_k%02d_n1%03d_n2%03d_%s_", k, n1, n2, scen_lab),
    tmpdir  = tempdir()
  )
  
  N_CORES <- min(parallel::detectCores(), n_rep)
  
  reg <- batchtools::makeRegistry(
    file.dir  = regdir,
    packages  = c("Rcpp", "dplyr", "tibble", "tidyr", "readr", "purrr")
  )
  reg$cluster.functions <- batchtools::makeClusterFunctionsMulticore(ncpus = N_CORES)
  
  batchtools::batchMap(
    fun = mw_rep_task,
    rep_id = 1:n_rep,
    more.args = list(
      n1              = n1,
      n2              = n2,
      k               = k,
      alpha           = alpha,
      assoc           = assoc,
      effe            = effe,
      sim_MC_per_job  = sim_MC_per_job,
      B_boot          = B_boot,
      R_perm          = R_perm
    ),
    reg = reg
  )
  
  ids <- batchtools::findJobs(reg = reg)
  batchtools::submitJobs(ids = ids, resources = list(ncpus = 1L), reg = reg)
  cat(" Running MW:", regdir, "on", N_CORES, "cores...\n")
  
  ok <- batchtools::waitForJobs(ids = ids, reg = reg)
  cat("waitForJobs result:", ok, "for", regdir, "\n")
  
  if (!ok) {
    cat(" Failed jobs in", regdir, ":\n")
    print(batchtools::getErrorMessages(reg = reg))
    # Try to remove the registry anyway
    try(batchtools::removeRegistry(reg = reg), silent = TRUE)
    stop("Abort: fix job errors before aggregating.")
  }
  
  done_ids <- batchtools::findDone(reg = reg)
  res_list <- batchtools::reduceResultsList(ids = done_ids, reg = reg)
  
  # Remove the registry folder (and everything inside)
  try(batchtools::removeRegistry(reg = reg), silent = TRUE)
  
  if (length(res_list) == 0L) {
    warning("No results in ", regdir)
    return(tibble())
  }
  
  detail_all <- dplyr::bind_rows(res_list)
  
  cat(" Completed MW cell: k=", k,
      ", (n1,n2)=(", n1, ",", n2, ")",
      " — detail=", nrow(detail_all), "\n", sep = "")
  
  detail_all
}

## ----------------------------------------------------------------
## 5) MAIN LOOP over k and (n1, n2) (H0 + H1)
## ----------------------------------------------------------------

detail_global_list <- list()

for (k_i in k_list) {
  assoc_ref_i <- get_assoc_ref(k_i)
  
  message("\n=======================")
  message(sprintf(">>> MW: k = %d", k_i))
  message("=======================\n")
  
  # k-specific (n1,n2) grid
  mw_pairs_k <- get_mw_pairs(k_i)
  
  for (row_idx in seq_len(nrow(mw_pairs_k))) {
    n1_cur <- mw_pairs_k$n1[row_idx]
    n2_cur <- mw_pairs_k$n2[row_idx]
    
    message(sprintf(">>> [H1] k=%d, n1=%d, n2=%d", k_i, n1_cur, n2_cur))
    
    # H1: assoc > 0, effect decreases with n1+n2 but power increases
    assoc_H1_cur <- assoc_H1_fun(n1_cur, n2_cur, assoc_ref = assoc_ref_i)
    
    res_H1 <- run_one_pair_batch_mw(
      k                = k_i,
      n1               = n1_cur,
      n2               = n2_cur,
      alpha            = alpha_test,
      assoc            = assoc_H1_cur,
      effe             = TRUE,
      sim_MC_per_job   = sim_MC_per_job,
      n_rep            = n_rep,
      B_boot           = B_boot,
      R_perm           = R_perm,
      root_dir         = ROOT
    )
    
    message(sprintf(">>> [H0] k=%d, n1=%d, n2=%d", k_i, n1_cur, n2_cur))
    
    # H0: assoc = 0, no effect
    res_H0 <- run_one_pair_batch_mw(
      k                = k_i,
      n1               = n1_cur,
      n2               = n2_cur,
      alpha            = alpha_test,
      assoc            = assoc_H0,
      effe             = FALSE,
      sim_MC_per_job   = sim_MC_per_job,
      n_rep            = n_rep,
      B_boot           = B_boot,
      R_perm           = R_perm,
      root_dir         = ROOT
    )
    
    detail_global_list <- c(detail_global_list, list(res_H1, res_H0))
  }
}

cat("\n All combinations (k, n1, n2) completed (MW).\n")

## -----------------------------------------------------------------
## 6) Aggregate everything into a single GLOBAL DETAIL file (MW)
## -----------------------------------------------------------------

detail_global <- if (length(detail_global_list)) dplyr::bind_rows(detail_global_list) else tibble()

if (nrow(detail_global) > 0L) {
  write_csv(detail_global, file.path(ROOT, "mw_DETAIL_GLOBAL_dgp1.csv"))
}

cat("    - Detail global rows (MW):", nrow(detail_global), "\n")
