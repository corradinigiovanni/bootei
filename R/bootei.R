#' BOOTEI permutation tests with coupled bootstrap tie-breaking
#'
#' @param x,y Inputs (see README).
#' @param test One of "chisq", "mannwhitney", "spearman", "kruskalwallis".
#' @param B Number of bootstrap replicates (B <= 1 gives classical permutation p-value).
#' @param R Number of permutations.
#' @param alternative "two.sided", "greater", or "less".
#' @param perm_seed Optional seed used ONLY inside bootei (does not affect caller RNG).
#' @param midp Logical; mid-p correction (only used when B <= 1).
#' @param boot_type "sobol", "sobol_shift", or "efron".
#' @param keep_perm_stats Logical; keep permutation stats.
#' @export
bootei <- function(x, y,
                   test = "chisq",
                   B = 100L,
                   R = 1000L,
                   alternative = "two.sided",
                   perm_seed = NA_real_,
                   midp = FALSE,
                   boot_type = "sobol",
                   keep_perm_stats = FALSE) {

  if (!is.na(perm_seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit({
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
    }

    set.seed(as.integer(perm_seed))
  }

  # chiama la funzione Rcpp generata (non editare RcppExports.R)
  bootei_cpp(x, y,
             test = test,
             B = as.integer(B),
             R = as.integer(R),
             alternative = alternative,
             perm_seed = as.numeric(perm_seed),
             midp = midp,
             boot_type = boot_type,
             keep_perm_stats = keep_perm_stats)
}
