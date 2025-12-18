#' BOOTEI permutation tests
#'
#' Permutation-based inference with optional coupled bootstrap tie-breaking (BOOTEI).
#'
#' @param x,y Input vectors (numeric/integer/factor/character depending on `test`).
#' @param test One of 'chisq', 'spearman', 'mannwhitney', 'kruskalwallis'.
#' @param B Bootstrap budget for tie-breaking (B<=1 = classical permutation).
#' @param R Number of Monte Carlo permutations.
#' @param alternative 'two.sided', 'greater', or 'less'.
#' @param perm_seed Optional seed for reproducibility (does NOT affect caller RNG).
#' @param midp Logical; mid-p only when B<=1.
#' @param boot_type 'sobol', 'sobol_shift', or 'efron'.
#' @param keep_perm_stats Logical; keep permutation stats.
#' @return A list with statistic_raw, statistic, tie_breaker, p.value, method, alternative.
#' @export
bootei <- function(x, y, test = "chisq", B = 100L, R = 1000L,
                   alternative = "two.sided", perm_seed = NA_real_,
                   midp = FALSE, boot_type = "sobol", keep_perm_stats = FALSE) {

  if (!is.na(perm_seed)) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL

    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }
    }, add = TRUE)

    set.seed(as.integer(perm_seed))
  }

  bootei_cpp(
    x, y,
    test = test, B = as.integer(B), R = as.integer(R),
    alternative = alternative,
    perm_seed = NA_real_,          # IMPORTANT: C++ must NOT touch RNG state
    midp = midp, boot_type = boot_type,
    keep_perm_stats = keep_perm_stats
  )
}
