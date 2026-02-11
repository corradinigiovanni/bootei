## NOTE: Script di simulazione/replica (non eseguito durante R CMD check).
## Questo script LANCIA le simulazioni e PRODUCE il file globale
##   welch_DETAIL_GLOBAL_dgp2.csv
## (non è un post-processing che legge un file già esistente).

# ================================================================
# BOOTEI / Permutation per Welch t-test (Sobol QMC) – DGP2
# Simulazioni su esiti ordinali (Likert) con k = 3, 5, 7 livelli.
#
# Test:
# - Welch two-sample t-test implementato via bootei() con test = "welch"
# - alternative = "two.sided"
#
# Metodi confrontati (p-value per dataset):
# 1) Permutation p-value con B = 1 (add-one correction; midp = FALSE)
# 2) mid-p permutation p-value con B = 1 (midp = TRUE)
# 3) BOOTEI p-value con B = B_boot > 1 (bootstrap/tie-breaking interno),
#    con boot_type = "sobol" (QMC Sobol per le repliche bootstrap)
#
# DGP2 (alternativa di tipo Lehmann su scala latente, poi discretizzata):
# - Due gruppi indipendenti di dimensione n1 e n2
# - Generazione latente:
#     U ~ Unif(0,1)
#     U_theta = U^(1/theta)
#     Z = F^{-1}(U_theta) con F^{-1} = qnorm() (dist = "normal")
# - Discretizzazione in k categorie comuni (cutpoints su quantili Normali):
#     cuts = qnorm(seq(0,1,length.out = k+1))
#     Y = cut(Z; breaks = cuts, labels = FALSE, include.lowest = TRUE)
#
# Scenari:
# - H0 (Type-I): effe = FALSE e assoc = 0  => theta = 1
# - H1 (Power):  effe = TRUE  e assoc = assoc_ref(k) usato come "strength"
#   dove theta è determinata da:
#     theta = 1 + assoc_ref(k) * ((n1+n2)/n_min)^(-gamma_assoc)
#   con:
#     gamma_assoc = 0.495
#     n_min = 10
#     assoc_ref(k): k=3 -> 2, k=5 -> 3, k=7 -> 4
#   (alternativa locale: l'effetto decresce con n1+n2 secondo la potenza -gamma)
#
# Design (griglie campionarie, bilanciate n1=n2, specifiche per k):
# - k = 3: (5,5), (10,10), (15,15), (25,25), (50,50), (100,100), (250,250), (500,500)
# - k = 5: (5,5), (10,10), (15,15), (25,25), (50,50), (100,100), (250,250)
# - k = 7: (5,5), (10,10), (15,15), (25,25), (50,50), (100,100), (250,250)
#
# Budget computazionale (come in questo script):
# - sim_MC_per_job = 200  dataset Monte Carlo per job
# - n_rep          = 25   job per cella (k, n1, n2, scenario)
# - R_perm         = 5000 permutazioni per dataset
# - B_boot         = 200  repliche bootstrap BOOTEI (Sobol QMC)
#
# Parallelizzazione:
# - batchtools: una registry temporanea per ogni cella (k,n1,n2,scenario)
# - backend multicore con N_CORES = min(detectCores(), n_rep)
# - le registry sono create in cartelle temporanee e rimosse a fine cella
#
# Output:
# - Un unico file globale in ROOT, contenente record per-dataset per H0 e H1:
#     welch_DETAIL_GLOBAL_dgp2.csv
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
## 0) PERCORSI, COMPILAZIONE C++ E PARAMETRI GLOBALI
## ----------------------------------------------------------------


ROOT <- "_"
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

library(bootei)


set.seed(123)

alpha_test     <- 0.05
assoc_H0       <- 0           # Scenario H0: nessun effetto (usato insieme a effe = FALSE)

sim_MC_per_job   <- 200L      # Numero di dataset Monte Carlo generati per ciascun job (replica)
n_rep            <- 25L       # Numero di job (repliche) per ogni cella (k, n1, n2, scenario)
B_boot           <- 200L      # Numero di repliche bootstrap BOOTEI (QMC Sobol) quando B > 1
R_perm           <- 5000L     # Numero di permutazioni per il calcolo del p-value (per dataset)

# Parametro dell'alternativa locale: controllo del decadimento dell'effetto con n1+n2
gamma_assoc    <- 0.495

# Valori di riferimento assoc_ref(k): controllano l'intensità dell'allontanamento da H0 per ciascun k
assoc_H1_ref_3 <- 2
assoc_H1_ref_5 <- 3
assoc_H1_ref_7 <- 4

# Lista dei livelli Likert considerati
k_list <- c(3L, 5L, 7L)

# ----------------------------------------------------------------
# Griglia (n1,n2) bilanciata e specifica per k
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


# assoc_ref specifico per ciascun k (usato per costruire theta in H1)
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

theta_fun <- function(n1, n2,
                      assoc_ref,           # Intensità base dell'effetto (dipende da k)
                      gamma = gamma_assoc,
                      n_min = 10) {
  n_tot <- n1 + n2
  1 + assoc_ref * (n_tot / n_min)^(-gamma)
}

## ----------------------------------------------------------------
## 1) DGP latente + discretizzazione Likert e wrapper BOOTEI per Welch
## ----------------------------------------------------------------

# Utility: forza BLAS/OMP/MKL ecc. a 1 thread per evitare oversubscription in multicore
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

# Generatore Lehmann su scala latente: U^(1/theta) e poi quantile della distribuzione scelta
r_lehmann <- function(n, theta = 1, dist = c("normal", "logistic")) {
  dist <- match.arg(dist)
  U <- runif(n)
  U2 <- if (theta == 1) U else U^(1/theta)
  if (dist == "normal") {
    qnorm(U2)
  } else {
    qlogis(U2)
  }
}

# Discretizzazione in k categorie usando cutpoints comuni basati su quantili Normali
simulate_group_likert_lehmann <- function(n, k = 5, theta = 1) {
  z    <- r_lehmann(n, theta = theta)
  cuts <- qnorm(seq(0, 1, length.out = k + 1))
  as.numeric(cut(z, breaks = cuts, labels = FALSE, include.lowest = TRUE))
}

# DGP a due gruppi:
# - se effe = FALSE oppure assoc = 0 => theta = 1 (H0)
# - altrimenti theta viene costruito via theta_fun usando assoc come assoc_ref (H1)
simulate_two_group <- function(n1, n2, k, assoc = 0, effe = TRUE) {
  if (!effe || assoc == 0) {
    theta <- 1
  } else {
    theta <- theta_fun(n1, n2, assoc_ref = assoc)  # assoc qui è l'assoc_ref(k) per H1
  }
  
  x <- simulate_group_likert_lehmann(n1, k, theta = 1)      # gruppo controllo (theta=1)
  y <- simulate_group_likert_lehmann(n2, k, theta = theta)  # gruppo trattamento (theta calcolato)
  
  list(x = x, y = y)
}

# Wrapper robusto: chiama bootei() con test="welch" e restituisce p.value; NA in caso di errore
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
## 2) Esecuzione di UNA cella Welch: (n1,n2,k,scenario) per 'simulations' dataset
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
  # Qui B=1 corrisponde al caso "solo permutazioni" (no bootstrap BOOTEI)
  B_perm    <- 1L
  B_midp    <- 1L
  midp_perm <- FALSE
  midp_midp <- TRUE
  
  # Vettori p-value (uno per dataset Monte Carlo; NAs se errori persistenti)
  p_perm <- rep(NA_real_, simulations)
  p_midp <- rep(NA_real_, simulations)
  p_boot <- rep(NA_real_, simulations)
  
  # Loop sui dataset Monte Carlo con retry (max_retry_per_sim) per gestire eventuali NA/edge cases
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
      
      dat <- simulate_two_group(n1 = n1, n2 = n2, k = k, assoc = assoc, effe = effe)
      x <- as.numeric(dat$x)
      y <- as.numeric(dat$y)
      
      # Seed per dataset e tentativo (riproducibile, ma differenziato per i tentativi)
      seed_i <- as.numeric(seed_base + i * 999 + attempts)
      
      # p-value permutazionale (add-one correction; midp=FALSE)
      p_perm_i <- bootei_safe_p_welch(
        x, y, B = B_perm, R = R_perm, midp = midp_perm, seed = seed_i
      )
      # p-value permutazionale mid-p (midp=TRUE)
      p_midp_i <- bootei_safe_p_welch(
        x, y, B = B_midp, R = R_perm, midp = midp_midp, seed = seed_i
      )
      # BOOTEI p-value con B = B_boot (bootstrap/tie-breaking interno; midp=FALSE)
      p_boot_i <- bootei_safe_p_welch(
        x, y, B = B_boot, R = R_perm, midp = FALSE, seed = seed_i
      )
      
      # Se uno dei tre p-value è NA, rigenera il dataset (retry)
      if (anyNA(c(p_perm_i, p_midp_i, p_boot_i))) next
      
      p_perm[i] <- p_perm_i
      p_midp[i] <- p_midp_i
      p_boot[i] <- p_boot_i
      i <- i + 1L
      break
    }
  }
  
  # Tieni solo i dataset per cui tutti e tre i p-value sono disponibili
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
  
  # Output per-dataset: p-value + indicatori di significatività a livello alpha
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
## 3) Funzione "job" per batchtools: una replica (rep_id) di una cella
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
  
  # Esegue sim_MC_per_job dataset per questa replica e produce record dettagliati
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
  
  # Etichetta scenario coerente con gli argomenti passati (H0 se effe=FALSE o assoc=0)
  scenario_lab <- if (assoc == 0 || !effe) "H0 (Type-I)" else "H1 (Power)"
  
  # Aggiunge metadati della cella e della replica
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
## 4) Esecuzione batch di UNA cella (k, n1, n2, scenario) con batchtools
##    - registry in tempdir()
##    - nessun file per-cella in ROOT: si aggrega in memoria e si rimuove la registry
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
  
  # Etichetta sintetica scenario usata solo per nominare la directory temporanea
  scen_lab <- if (assoc == 0 || !effe) "H0" else "H1"
  
  regdir <- tempfile(
    pattern = sprintf("welch_k%02d_n1%03d_n2%03d_%s_", k, n1, n2, scen_lab),
    tmpdir  = tempdir()
  )
  
  # Numero di core usati dal backend multicore (limitato da n_rep)
  N_CORES <- min(parallel::detectCores(), n_rep)
  
  reg <- batchtools::makeRegistry(
    file.dir  = regdir,
    packages  = c("Rcpp", "dplyr", "tibble", "tidyr", "readr", "purrr")
  )
  reg$cluster.functions <- batchtools::makeClusterFunctionsMulticore(ncpus = N_CORES)
  
  # Mappa: un job per rep_id = 1..n_rep (ogni job genera sim_MC_per_job dataset)
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
  
  # Se ci sono errori, stampa i messaggi e interrompe (evitando aggregazioni inconsistenti)
  if (!ok) {
    cat(" Failed jobs in", regdir, ":\n")
    print(batchtools::getErrorMessages(reg = reg))
    try(batchtools::removeRegistry(reg = reg), silent = TRUE)
    stop("Abort: fix job errors before aggregating.")
  }
  
  # Aggrega i risultati dei job completati (lista di tibble) e rimuove la registry temporanea
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
## 5) LOOP PRINCIPALE su k e (n1,n2): esegue H1 e H0 per ogni coppia
## ----------------------------------------------------------------

detail_global_list <- list()

for (k_i in k_list) {
  # assoc_ref per questo k (usato solo in H1 come parametro "strength")
  assoc_ref_i <- get_assoc_ref(k_i)
  
  message("\n=======================")
  message(sprintf(">>> Welch: k = %d", k_i))
  message("=======================\n")
  
  pairs_k <- get_pairs(k_i)
  
  for (row_idx in seq_len(nrow(pairs_k))) {
    n1_cur <- pairs_k$n1[row_idx]
    n2_cur <- pairs_k$n2[row_idx]
    
    message(sprintf(">>> [H1] k=%d, n1=%d, n2=%d", k_i, n1_cur, n2_cur))
    
    # H1: assoc è impostato all'assoc_ref(k) (da cui theta viene calcolato in simulate_two_group)
    assoc_H1_cur <- assoc_ref_i
    
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
    
    # H0: nessun effetto (theta=1) forzato da effe=FALSE e assoc=0
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
    
    # Accumula i dettagli (H1 + H0) per poi fare bind_rows finale
    detail_global_list <- c(detail_global_list, list(res_H1, res_H0))
  }
}

cat("\n All combinations (k, n1, n2) completed (Welch).\n")

## -----------------------------------------------------------------
## 6) AGGREGAZIONE FINALE: salva il file globale con tutti i record per-dataset
## -----------------------------------------------------------------

detail_global <- if (length(detail_global_list)) dplyr::bind_rows(detail_global_list) else tibble()

if (nrow(detail_global) > 0L) {
  write_csv(detail_global, file.path(ROOT, "welch_DETAIL_GLOBAL_dgp2.csv"))
}

cat("    - Detail global rows (Welch):", nrow(detail_global), "\n")
