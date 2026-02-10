## Replication script (not run during R CMD check).
## Runs BOOTEI / permutation Welch t-test simulations for DGP 1 and writes
## a single per-dataset “detail” file: welch_DETAIL_GLOBAL_dgp1.csv.
## Requires bootei.cpp and bootei.R at the paths set below.

# ================================================================
# BOOTEI Welch t-test (Sobol QMC) – DGP 1
# Simulations over balanced (n1 = n2) designs for ordinal (Likert-type) outcomes
# with k ∈ {3, 5, 7} ordered categories.
#
# - Test statistic: Welch two-sample t statistic (bootei test = "welch")
# - Methods recorded (per dataset):
#     * Permutation p-value via bootei with B = 1, midp = FALSE
#     * mid-p permutation p-value via bootei with B = 1, midp = TRUE
#     * BOOTEI p-value via bootei with B = B_boot > 1 (Sobol QMC bootstrap),
#       with internal tie-breaking as implemented in bootei.cpp
#
# - Engine: bootei() (C++/Rcpp) with:
#     * test = "welch"
#     * boot_type = "sobol"
#
# - DGP (latent-Likert):
#     * Two independent groups of sizes n1 and n2
#     * Latent variables:
#         Z1 ~ N(0, 1),  Z2 ~ N(mu2, 1)
#       Discretised using common Normal-quantile cutpoints into k ordered categories,
#       where cutpoints are qnorm(seq(0,1,length.out=k+1)) (i.e., equal-probability
#       bins under N(0,1)).
#     * H0 implementation: mu2 = 0 (effe = FALSE and assoc = 0)
#     * H1 implementation: mu2 = assoc (effe = TRUE), with
#         assoc = assoc_ref(k) * (n1 + n2)^(-gamma_assoc)
#       where gamma_assoc = 0.4 in this script.
#
# - Design (balanced grids; n1 = n2):
#     * k = 3: n ∈ {5, 10, 15, 25, 50, 100, 250, 500}
#     * k = 5: n ∈ {5, 10, 15, 25, 50, 100, 250}
#     * k = 7: n ∈ {5, 10, 15, 25, 50, 100, 250}
#
# - Monte Carlo / computational budget (defaults below):
#     * sim_MC_per_job = 200 datasets per job (replicate)
#     * n_rep          = 25 jobs per (k, n1, n2, scenario) cell
#       → nominally 200 * 25 = 5000 datasets per cell (before dropping failures)
#     * R_perm         = 5000 permutations per dataset
#     * B_boot         = 200 BOOTEI bootstrap replicates (Sobol QMC)
#
# - Parallelisation:
#     * batchtools registry created in a temporary folder per cell
#     * multicore backend with up to min(detectCores(), n_rep) workers
#     * each rep_id is a separate batchtools job (submitted with ncpus = 1)
#
# - Output:
#     * One global CSV with per-dataset records for both H0 and H1:
#         welch_DETAIL_GLOBAL_dgp1.csv
#     * Temporary batchtools registries are deleted after each cell finishes.
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
## 0) ROOT PATHS, COMPILATION, AND GLOBAL SETTINGS
## ----------------------------------------------------------------

ROOT <- "_"
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

library(bootei)
  
set.seed(123)

alpha_test     <- 0.05
assoc_H0       <- 0           # H0 uses no latent shift (also enforced via effe = FALSE)

sim_MC_per_job   <- 200L      # Monte Carlo datasets generated within each batch job
n_rep            <- 25L       # number of independent batch jobs per (k, n1, n2, scenario) cell
B_boot           <- 200L      # BOOTEI bootstrap replicates (Sobol QMC)
R_perm           <- 5000L     # number of permutations per dataset (Welch statistic)

# Local-alternative scaling exponent used in assoc_H1_fun()
gamma_assoc    <- 0.4

# Reference effect magnitudes (one per k) used to define assoc under H1
assoc_H1_ref_3 <- 7/3
assoc_H1_ref_5 <- 8/3
assoc_H1_ref_7 <- 9/3

# Likert levels (number of ordered categories)
k_list <- c(3L, 5L, 7L)

# ----------------------------------------------------------------
# Balanced sample-size grids (n1 = n2), specified separately per k
# ----------------------------------------------------------------

get_pairs <- function(k) {
  if (k == 3L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      15L,  15L,
      25L,  25L,
      50L,  50L,
      100L, 100L,
      250L, 250L,
      500L, 500L
    )
  } else if (k == 5L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      15L,  15L,
      25L,  25L,
      50L,  50L,
      100L, 100L,
      250L, 250L
    )
  } else if (k == 7L) {
    tibble::tribble(
      ~n1, ~n2,
      5L,   5L,
      10L,  10L,
      15L,  15L,
      25L,  25L,
      50L,  50L,
      100L, 100L,
      250L, 250L
    )
  } else {
    stop("sample-size grid not defined for k = ", k)
  }
}

# Return assoc_ref(k) used to scale the H1 latent shift for each k
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

# Local-alternative shift under H1 (depends on total sample size and assoc_ref)
assoc_H1_fun <- function(n1, n2,
                         assoc_ref,
                         gamma = gamma_assoc) {
  assoc_ref * (n1 + n2) ** (-1*gamma)
}

## ----------------------------------------------------------------
## 1) LATENT-LIKERT DGP AND A SAFE bootei() WRAPPER FOR WELCH
## ----------------------------------------------------------------

# Utility: force single-threaded BLAS/OMP to avoid oversubscription under multicore jobs
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

# Simulate one group from a latent N(mu,1) and discretise into k ordered categories
# using common Normal-quantile cutpoints (includes -Inf and +Inf endpoints).
simulate_group_likert_latent <- function(n, k = 5, mu = 0) {
  z <- rnorm(n, mean = mu, sd = 1)
  cuts <- qnorm(seq(0, 1, length.out = k + 1))
  as.numeric(cut(z, breaks = cuts, labels = FALSE, include.lowest = TRUE))
}

# Two-group generator: shift is applied to group 2 only when effe = TRUE and assoc != 0
simulate_two_group_likert <- function(n1, n2, k, assoc = 0, effe = TRUE) {
  if (!effe || assoc == 0) {
    mu1 <- 0
    mu2 <- 0
  } else {
    mu1 <- 0
    mu2 <- assoc   # latent location shift for group 2
  }
  g1 <- simulate_group_likert_latent(n1, k, mu = mu1)
  g2 <- simulate_group_likert_latent(n2, k, mu = mu2)
  list(x = g1, y = g2)
}

# Safe wrapper: call bootei(test="welch") and return p.value; return NA on any error
bootei_safe_p_welch <- function(x, y,
                                B,
                                R,
                                midp,
                                seed) {
  tryCatch({
    bootei(
      x           = as.numeric(x),
      y           = as.numeric(y),
      test        = "welch",
      B           = as.integer(B),
      R           = as.integer(R),
      alternative = "two.sided",
      perm_seed   = as.numeric(seed),
      midp        = midp,
      boot_type   = "sobol"
    )$p.value
  }, error = function(e) NA_real_)
}

## ----------------------------------------------------------------
## 2) RUN ONE (k, n1, n2, assoc) CELL LOCALLY (ONE JOB'S SHARE)
##    - Generates 'simulations' datasets
##    - Computes three p-values per dataset (perm, mid-p, BOOTEI)
##    - Retries simulation up to max_retry_per_sim times if bootei() errors/returns NA
##    - Drops datasets with any NA p-values before returning
## ----------------------------------------------------------------

run_cell_welch <- function(n1,
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
      
      dat <- simulate_two_group_likert(n1 = n1, n2 = n2, k = k, assoc = assoc, effe = effe)
      x <- as.numeric(dat$x)
      y <- as.numeric(dat$y)
      
      # Dataset-specific seed passed to bootei() (perm_seed); varies with i and retry counter
      seed_i <- as.numeric(seed_base + i * 999 + attempts)
      
      p_perm_i <- bootei_safe_p_welch(
        x, y, B = B_perm, R = R_perm, midp = midp_perm, seed = seed_i
      )
      p_midp_i <- bootei_safe_p_welch(
        x, y, B = B_midp, R = R_perm, midp = midp_midp, seed = seed_i
      )
      p_boot_i <- bootei_safe_p_welch(
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
## 3) batchtools JOB FUNCTION (ONE rep_id)
##    - Pins threads, sets an R seed for this job
##    - Runs run_cell_welch() for sim_MC_per_job datasets
##    - Attaches scenario metadata (H0 vs H1) and design parameters
## ----------------------------------------------------------------

welch_rep_task <- function(rep_id,
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
  
  detail <- run_cell_welch(
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
## 4) RUN ALL reps FOR ONE (k, n1, n2, scenario) CELL VIA batchtools
##    - Creates a temporary registry per cell (deleted after completion)
##    - Maps welch_rep_task over rep_id = 1:n_rep
##    - Collects all per-dataset records into a single tibble for this cell
## ----------------------------------------------------------------

run_one_pair_batch_welch <- function(k,
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
    pattern = sprintf("welch_k%02d_n1%03d_n2%03d_%s_", k, n1, n2, scen_lab),
    tmpdir  = tempdir()
  )
  
  N_CORES <- min(parallel::detectCores(), n_rep)
  
  reg <- batchtools::makeRegistry(
    file.dir  = regdir,
    packages  = c("Rcpp", "dplyr", "tibble", "tidyr", "readr", "purrr")
  )
  reg$cluster.functions <- batchtools::makeClusterFunctionsMulticore(ncpus = N_CORES)
  
  batchtools::batchMap(
    fun = welch_rep_task,
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
  cat(" Running Welch:", regdir, "on", N_CORES, "cores...\n")
  
  ok <- batchtools::waitForJobs(ids = ids, reg = reg)
  cat("waitForJobs result:", ok, "for", regdir, "\n")
  
  if (!ok) {
    cat(" Failed jobs in", regdir, ":\n")
    print(batchtools::getErrorMessages(reg = reg))
    try(batchtools::removeRegistry(reg = reg), silent = TRUE)
    stop("Abort: fix job errors before aggregating.")
  }
  
  done_ids <- batchtools::findDone(reg = reg)
  res_list <- batchtools::reduceResultsList(ids = done_ids, reg = reg)
  
  try(batchtools::removeRegistry(reg = reg), silent = TRUE)
  
  if (length(res_list) == 0L) {
    warning("No results in ", regdir)
    return(tibble())
  }
  
  detail_all <- dplyr::bind_rows(res_list)
  
  cat(" Completed Welch cell: k=", k,
      ", (n1,n2)=(", n1, ",", n2, ")",
      " — detail=", nrow(detail_all), "\n", sep = "")
  
  detail_all
}

## ----------------------------------------------------------------
## 5) MAIN LOOP OVER k AND (n1, n2)
##    For each balanced pair:
##      - Run H1 with assoc = assoc_ref(k) * (n1+n2)^(-gamma_assoc) and effe = TRUE
##      - Run H0 with assoc = 0 and effe = FALSE
##    Accumulate all per-dataset records in memory before writing the global file.
## ----------------------------------------------------------------

detail_global_list <- list()

for (k_i in k_list) {
  assoc_ref_i <- get_assoc_ref(k_i)
  
  message("\n=======================")
  message(sprintf(">>> Welch: k = %d", k_i))
  message("=======================\n")
  
  pairs_k <- get_pairs(k_i)
  
  for (row_idx in seq_len(nrow(pairs_k))) {
    n1_cur <- pairs_k$n1[row_idx]
    n2_cur <- pairs_k$n2[row_idx]
    
    message(sprintf(">>> [H1] k=%d, n1=%d, n2=%d", k_i, n1_cur, n2_cur))
    
    assoc_H1_cur <- assoc_H1_fun(n1_cur, n2_cur, assoc_ref = assoc_ref_i)
    
    res_H1 <- run_one_pair_batch_welch(
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
    
    res_H0 <- run_one_pair_batch_welch(
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

cat("\n All combinations (k, n1, n2) completed (Welch).\n")

## -----------------------------------------------------------------
## 6) AGGREGATE ALL CELLS AND WRITE A SINGLE GLOBAL DETAIL CSV
##    The output contains per-dataset p-values and rejection indicators for all
##    (k, n1, n2) under both H0 and H1, with rep_id labels.
## -----------------------------------------------------------------

detail_global <- if (length(detail_global_list)) dplyr::bind_rows(detail_global_list) else tibble()

if (nrow(detail_global) > 0L) {
  write_csv(detail_global, file.path(ROOT, "welch_DETAIL_GLOBAL_dgp1.csv"))
}

cat("    - Detail global rows (Welch):", nrow(detail_global), "\n")
