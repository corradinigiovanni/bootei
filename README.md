# bootei

BOOTEI permutation tests with coupled bootstrap tie-breaking (Rcpp/C++).

## Install

```r
install.packages("remotes")
remotes::install_github("condor-machine/bootei")  # change if repo name differs
```

## Examples (classical vs BOOTEI)

Classical permutation: `B = 1`  
BOOTEI tie-breaking: `B > 1` (e.g., `B = 200`)

These toy examples use *discrete/ordinal data* to induce many ties under permutations (often making BOOTEI < classical).

### 1) Chi-square (categorical)

```r
x <- factor(c(rep("A", 14), rep("B", 14), rep("C", 12)))
y <- factor(c(rep("yes", 10), rep("no", 30)))

p_perm   <- bootei::bootei(x, y, test = "chisq", B = 1,   R = 5000, perm_seed = 1)$p.value
p_bootei <- bootei::bootei(x, y, test = "chisq", B = 200, R = 5000, perm_seed = 1,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)
```

### 2) Mann–Whitney (ordinal, ties)

```r
x <- c(rep(0, 18), rep(1, 6), rep(2, 6))
y <- c(rep(0,  8), rep(1, 8), rep(2, 14))

p_perm   <- bootei::bootei(x, y, test = "mannwhitney", B = 1,   R = 5000,
                           alternative = "two.sided", perm_seed = 1)$p.value
p_bootei <- bootei::bootei(x, y, test = "mannwhitney", B = 200, R = 5000,
                           alternative = "two.sided", perm_seed = 1,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)
```

### 3) Spearman (rank ties)

```r
x <- c(rep(1, 12), rep(2, 12), rep(3, 12), rep(4, 12))
y <- c(rep(1, 18), rep(2,  6), rep(3, 14), rep(4, 10))

p_perm   <- bootei::bootei(x, y, test = "spearman", B = 1,   R = 5000,
                           alternative = "two.sided", perm_seed = 1)$p.value
p_bootei <- bootei::bootei(x, y, test = "spearman", B = 200, R = 5000,
                           alternative = "two.sided", perm_seed = 1,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)
```

### 4) Kruskal–Wallis (multiple groups, ties)

```r
g <- factor(rep(c("G1","G2","G3"), each = 20))
x <- c(rep(0, 14), rep(1, 6),    # G1
       rep(0, 10), rep(1, 10),   # G2
       rep(0,  6), rep(1, 14))   # G3

p_perm   <- bootei::bootei(x, g, test = "kruskalwallis", B = 1,   R = 5000,
                           alternative = "greater", perm_seed = 1)$p.value
p_bootei <- bootei::bootei(x, g, test = "kruskalwallis", B = 200, R = 5000,
                           alternative = "greater", perm_seed = 1,
                           boot_type = "sobol")$p.value
c(p_perm = p_perm, p_bootei = p_bootei)
```

## Notes

- `R` = number of Monte Carlo permutations (increase for final runs).
- `B > 1` activates BOOTEI tie-breaking; `B = 1` is classical permutation (optional `midp = TRUE`).
- `boot_type`: `"sobol"`, `"sobol_shift"`, `"efron"`.
- Use `perm_seed` for reproducibility.

## License

MIT License (see `LICENSE`).
