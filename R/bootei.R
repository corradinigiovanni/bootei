#' BOOTEI permutation tests
#'
#' Permutation-based inference with optional coupled bootstrap tie-breaking (BOOTEI).
#'
#' @param x,y Input vectors (numeric/integer/factor/character depending on `test`).
#' @param test One of 'chisq', 'spearman', 'mannwhitney', 'kruskalwallis'.
#' @param B Bootstrap budget for tie-breaking (B<=1 = classical permutation).
#' @param R Number of Monte Carlo permutations.
#' @param alternative 'two.sided', 'greater', or 'less'.
#' @param perm_seed Optional seed for reproducibility.
#' @param midp Logical; mid-p only when B<=1.
#' @param boot_type 'sobol', 'sobol_shift', or 'efron'.
#' @param keep_perm_stats Logical; keep permutation stats.
#' @return A list with statistic_raw, statistic, tie_breaker, p.value, method, alternative (and optionally perm stats).
#' @export
bootei <- function(x, y, test = "chisq", B = 100L, R = 1000L,
                   alternative = "two.sided", perm_seed = NA_real_,
                   midp = FALSE, boot_type = "sobol", keep_perm_stats = FALSE) {
  .Call(`_bootei_bootei`, x, y, test, as.integer(B), as.integer(R), alternative,
        as.numeric(perm_seed), midp, boot_type, keep_perm_stats)
}

