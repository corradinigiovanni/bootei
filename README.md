# bootei

BOOTEI (BOOTstrap Ensemble Inference) provides permutation tests with deterministic,
coupled bootstrap tie-breaking, designed for low-resolution regimes such as small
samples, sparse tables, few outcome levels, and pervasive ties.

The method preserves the original permutation test statistic and uses a fixed
bootstrap ensemble score solely to break permutation ties via a lexicographic
ordering. Under exchangeability, BOOTEI retains finite-sample validity and yields
permutation p-values that are never larger than those from the corresponding
classical permutation test.


## Install

```r
install.packages("remotes")
remotes::install_github("corradinigiovanni/bootei")
library(bootei)
```

## Examples (classical vs BOOTEI)

Classical permutation: `B = 1`  
BOOTEI tie-breaking: `B > 1` (e.g., `B = 200`)

These toy examples use *discrete or ordinal data* to induce many ties under permutations,
a setting in which BOOTEI deterministically refines the permutation ordering and
yields p-values that are never larger than the corresponding classical permutation p-values.

### 1) Chi-square (categorical)

```r
x <- sample(c(0,1), 10, replace = TRUE)
y <- sample(c(0,1), 10, replace = TRUE)

p_perm   <- bootei::bootei(x, y, test = "chisq", B = 1,   R = 5000, perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, y, test = "chisq", B = 200, R = 5000, perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)


```

### 2) Mann–Whitney (ordinal, ties)

```r
x <- sample(c(0,1,2), 5, replace = TRUE)
y <- sample(c(0,1,2), 5, replace = TRUE)

p_perm   <- bootei::bootei(x, y, test = "mannwhitney", B = 1,   R = 5000,
                           alternative = "two.sided", perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, y, test = "mannwhitney", B = 200, R = 5000,
                           alternative = "two.sided", perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)

```

### 3) Spearman (rank ties)

```r
x <- sample(c(0,1,2), 10, replace = TRUE)
y <- sample(c(0,1,2), 10, replace = TRUE)

p_perm   <- bootei::bootei(x, y, test = "spearman", B = 1,   R = 5000,
                           alternative = "two.sided", perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, y, test = "spearman", B = 200, R = 5000,
                           alternative = "two.sided", perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)

```

### 4) Kruskal–Wallis (multiple groups, ties)

```r
g <- factor(rep(c("G1","G2","G3"), each = 5))
x <- c(sample(c(0,1,2), 5, replace = TRUE),    # G1
       sample(c(0,1,2), 5, replace = TRUE),    # G2
       sample(c(0,1,2), 5, replace = TRUE))    # G3

p_perm   <- bootei::bootei(x, g, test = "kruskalwallis", B = 1,   R = 5000,
                           alternative = "greater", perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, g, test = "kruskalwallis", B = 200, R = 5000,
                           alternative = "greater", perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)

```

## Vignettes

An introductory vignette illustrating the use of BOOTEI on the Soybean dataset
is available via `browseVignettes("bootei")`.

## Notes

- `R` = number of Monte Carlo permutations (increase for final runs).
- `B > 1` activates BOOTEI tie-breaking; `B = 1` is classical permutation (optional `midp = TRUE`).
- `boot_type`: `"sobol"`, `"sobol_shift"`, `"efron"`.
- Use `perm_seed` for reproducibility.

## License

MIT License (see `LICENSE`).
