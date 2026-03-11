# bootei

BOOTEI (BOOTstrap Ensemble Inference) implements add-one Monte Carlo permutation
inference with deterministic bootstrap-based tie-breaking, designed for low-resolution
settings such as small samples, sparse or imbalanced contingency tables, few outcome
levels, and pervasive ties.

BOOTEI leaves the primary permutation test statistic unchanged and uses a fixed
bootstrap-plan ensemble score solely to break permutation ties via a lexicographic
ordering. Under exchangeability and the coupled fixed-plan construction, BOOTEI
retains finite-sample validity and, when evaluated on the same realised permutation
sample, yields permutation p-values that are never larger than those from the
corresponding classical add-one permutation test.

The paper focuses on Pearson's chi-square and Welch two-sample t permutation tests.
The package also includes implementations for the Mann–Whitney rank-sum test,
Spearman's rank correlation, and the Kruskal–Wallis test.

## Install

```r
install.packages("remotes")
remotes::install_github("corradinigiovanni/bootei")
library(bootei)
```

## Vignette
To install **bootei** *including vignettes*, use:
```r
remotes::install_github("corradinigiovanni/bootei", build_vignettes = TRUE)
```
Then you can open the vignette in R with:
```r
vignette("soybean", package = "bootei")
```

## Examples (classical vs BOOTEI)

Classical permutation: `B = 1`  
BOOTEI tie-breaking: `B > 1` (e.g., `B = 200`)

These toy examples use *discrete or ordinal data* to induce many ties under permutations,
a setting in which BOOTEI deterministically refines the permutation ordering and
yields p-values that are never larger (but still valid) than the corresponding classical permutation p-values.

### 1) Chi-square (categorical)

```r
x <- sample(c(0,1), 10, replace = TRUE)
y <- sample(c(0,1), 10, replace = TRUE)

p_perm   <- bootei::bootei(x, y, test = "chisq", B = 1,   R = 5000, perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, y, test = "chisq", B = 200, R = 5000, perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)


```

### 2) Welch t (two-sample, unequal variances)

```r
x <- sample(c(0,1,2), 5, replace = TRUE)
y <- sample(c(0,1,2), 5, replace = TRUE)

p_perm   <- bootei::bootei(x, y, test = "welch", B = 1,   R = 5000,
                           alternative = "two.sided", perm_seed = 910)$p.value
p_bootei <- bootei::bootei(x, y, test = "welch", B = 200, R = 5000,
                           alternative = "two.sided", perm_seed = 910,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)

```

### 3) Mann–Whitney (two-sample, ordinal)

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

### 4) Spearman (rank correlation)

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

### 5) Kruskal–Wallis (multiple groups, ordinal)

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


## Notes

- `R` = number of Monte Carlo permutations.
- `B > 1` activates BOOTEI tie-breaking; `B = 1` is classical permutation (optional `midp = TRUE`).
- `boot_type`: `"sobol"`, `"sobol_shift"`, `"efron"`.
- Use `perm_seed` for reproducibility.

## License

MIT License (see `LICENSE`).
