# ================================================================
# BOOTEI – Pearson chi-square: analysis of chi2_DETAIL_GLOBAL_dgp2.csv
#   - Power (H1): Perm vs BOOTEI, by n and table size
#   - McNemar tests under H1: BOOTEI > Perm; BOOTEI vs mid-p
#   - Size under H0: binom.test vs alpha (BH-adjusted and raw)
#   - Summary power figure with log10(n) x-axis
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(glue)
  library(conflicted)
  library(patchwork)
})

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")

# ---------------------------- PATH -----------------------------

ROOT <- "."
file_global <- file.path(ROOT, "chi2_DETAIL_GLOBAL_dgp2.csv")
stopifnot("chi2_DETAIL_GLOBAL_dgp2.csv not found!" = file.exists(file_global))

# --------------------- LOAD + CELL SUMMARIES -------------------

chi2_all <- read_csv(file_global, show_col_types = FALSE) %>%
  mutate(
    k_row = as.integer(k_row),
    k_col = as.integer(k_col),
    n     = as.integer(n),
    k_lab = glue("{k_row}x{k_col}")
  )

chi2_summary_all <- chi2_all %>%
  group_by(scenario, k_row, k_col, k_lab, n, alpha, theta) %>%
  summarise(
    total_sim = n(),
    rej_perm  = sum(sig_perm),
    rej_midp  = sum(sig_midp),
    rej_boot  = sum(sig_boot),
    rate_perm = rej_perm / total_sim,
    rate_midp = rej_midp / total_sim,
    rate_boot = rej_boot / total_sim,
    .groups   = "drop"
  )

alpha_nominal <- chi2_summary_all$alpha %>%
  unique() %>%
  sort() %>%
  tail(1)

cat(glue("Nominal alpha: {alpha_nominal}\n"))

# ------------------- LONG FORMAT (POWER, H1) -------------------

chi2_long_H1 <- chi2_summary_all %>%
  filter(grepl("H1", scenario)) %>%
  select(scenario, k_row, k_col, k_lab, n, theta, alpha,
         rate_perm, rate_boot) %>%
  pivot_longer(c(rate_perm, rate_boot),
               names_to = "method", values_to = "power") %>%
  mutate(
    method = recode(method, rate_perm = "Perm", rate_boot = "BOOTEI"),
    method = factor(method, levels = c("Perm", "BOOTEI"))
  )

# ---------------------- MCNEMAR TESTS (H1) ---------------------

mcnemar_binom <- function(b, c) {
  disc <- b + c
  if (disc <= 0) return(1.0)
  pbinom(b - 1, disc, 0.5, lower.tail = FALSE)
}

H1_mcnemar <- chi2_all %>%
  filter(grepl("H1", scenario)) %>%
  group_by(k_row, k_col, k_lab, n, theta) %>%
  summarise(
    total_sim  = n(),
    power_perm = mean(sig_perm),
    power_boot = mean(sig_boot),
    a          = sum(sig_boot &  sig_perm),
    b          = sum(sig_boot & !sig_perm),
    c          = sum(!sig_boot &  sig_perm),
    d          = sum(!sig_boot & !sig_perm),
    .groups    = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_mcn        = mcnemar_binom(b, c),   # BOOTEI > Perm
    diff_power   = power_boot - power_perm,
    sig_mcn_0.05 = is.finite(p_mcn) & p_mcn < 0.05
  ) %>%
  ungroup() %>%
  arrange(k_row, k_col, n)

H1_mcnemar_summary <- H1_mcnemar %>%
  group_by(k_row, k_col, k_lab) %>%
  summarise(
    n_cells        = n(),
    n_sig_raw      = sum(sig_mcn_0.05, na.rm = TRUE),
    first_sig_n    = if (any(sig_mcn_0.05, na.rm = TRUE))
      min(n[sig_mcn_0.05], na.rm = TRUE) else NA_integer_,
    last_sig_n     = if (any(sig_mcn_0.05, na.rm = TRUE))
      max(n[sig_mcn_0.05], na.rm = TRUE) else NA_integer_,
    max_diff_power = if (n_cells > 0) max(diff_power, na.rm = TRUE) else NA_real_,
    .groups        = "drop"
  )

H1_mcnemar_midp <- chi2_all %>%
  filter(grepl("H1", scenario)) %>%
  group_by(k_row, k_col, k_lab, n, theta) %>%
  summarise(
    total_sim   = n(),
    power_midp  = mean(sig_midp),
    power_boot  = mean(sig_boot),
    a           = sum(sig_boot &  sig_midp),
    b           = sum(sig_boot & !sig_midp),
    c           = sum(!sig_boot &  sig_midp),
    d           = sum(!sig_boot & !sig_midp),
    .groups     = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_boot_gt_midp  = mcnemar_binom(b, c),
    p_midp_gt_boot  = mcnemar_binom(c, b),
    diff_power      = power_boot - power_midp,
    sig_boot_gt_midp = is.finite(p_boot_gt_midp) & p_boot_gt_midp < 0.05,
    sig_midp_gt_boot = is.finite(p_midp_gt_boot) & p_midp_gt_boot < 0.05
  ) %>%
  ungroup() %>%
  arrange(k_row, k_col, n)

H1_mcnemar_midp_summary <- H1_mcnemar_midp %>%
  group_by(k_row, k_col, k_lab) %>%
  summarise(
    n_cells              = n(),
    n_sig_boot_gt_midp   = sum(sig_boot_gt_midp, na.rm = TRUE),
    n_sig_midp_gt_boot   = sum(sig_midp_gt_boot, na.rm = TRUE),
    first_n_boot_gt_midp = if (any(sig_boot_gt_midp, na.rm = TRUE))
      min(n[sig_boot_gt_midp], na.rm = TRUE) else NA_integer_,
    first_n_midp_gt_boot = if (any(sig_midp_gt_boot, na.rm = TRUE))
      min(n[sig_midp_gt_boot], na.rm = TRUE) else NA_integer_,
    max_diff_power       = if (n_cells > 0) max(diff_power, na.rm = TRUE) else NA_real_,
    .groups              = "drop"
  )

# -------------------- SIZE TESTS (H0): BH + RAW ----------------

h0_base <- chi2_summary_all %>%
  filter(grepl("H0", scenario)) %>%
  transmute(
    k_row = as.integer(k_row),
    k_col = as.integer(k_col),
    k_lab,
    n = as.integer(n),
    alpha,
    total_sim,
    rej_perm,
    rej_midp,
    rej_boot,
    rate_perm,
    rate_midp,
    rate_boot
  )

h0_long <- h0_base %>%
  pivot_longer(
    c(rej_perm, rej_midp, rej_boot, rate_perm, rate_midp, rate_boot),
    names_to  = c("what", "method_raw"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = what, values_from = value) %>%
  mutate(
    method = recode(method_raw, perm = "Perm", midp = "mid-p", boot = "BOOTEI")
  )

binom_tests_row <- function(rej, tot, alpha) {
  bt_two <- binom.test(rej, tot, p = alpha, alternative = "two.sided")
  tibble(
    p_two     = bt_two$p.value,
    ci_lo     = unname(bt_two$conf.int[1]),
    ci_hi     = unname(bt_two$conf.int[2]),
    p_greater = binom.test(rej, tot, p = alpha, alternative = "greater")$p.value,
    p_less    = binom.test(rej, tot, p = alpha, alternative = "less")$p.value
  )
}

H0_tests <- h0_long %>%
  rowwise() %>%
  mutate(tmp = list(binom_tests_row(rej, total_sim, alpha))) %>%
  unnest(cols = c(tmp)) %>%
  ungroup() %>%
  group_by(method) %>%
  mutate(
    p_two_BH     = p.adjust(p_two,     method = "BH"),
    p_greater_BH = p.adjust(p_greater, method = "BH"),
    p_less_BH    = p.adjust(p_less,    method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    size_flag_BH = case_when(
      p_greater_BH < 0.05 ~ "↑ (size > α, BH)",
      p_less_BH    < 0.05 ~ "↓ (size < α, BH)",
      TRUE                ~ "OK"
    )
  ) %>%
  arrange(k_row, k_col, n, method)

H0_summary <- H0_tests %>%
  group_by(method, k_row, k_col, k_lab) %>%
  summarise(
    n_cells       = n(),
    n_above_BH    = sum(p_greater_BH < 0.05, na.rm = TRUE),
    n_below_BH    = sum(p_less_BH    < 0.05, na.rm = TRUE),
    first_above_n = if (any(p_greater_BH < 0.05, na.rm = TRUE))
      min(n[p_greater_BH < 0.05], na.rm = TRUE) else NA_integer_,
    .groups       = "drop"
  )

H0_tests_noBH <- H0_tests %>%
  mutate(
    size_flag_noBH = case_when(
      p_greater < 0.05 ~ "↑ (size > α, raw)",
      p_less    < 0.05 ~ "↓ (size < α, raw)",
      TRUE             ~ "OK"
    )
  ) %>%
  arrange(k_row, k_col, n, method)

H0_summary_noBH <- H0_tests_noBH %>%
  group_by(method, k_row, k_col, k_lab) %>%
  summarise(
    n_cells        = n(),
    n_above_raw    = sum(p_greater < 0.05, na.rm = TRUE),
    n_below_raw    = sum(p_less    < 0.05, na.rm = TRUE),
    first_above_n  = if (any(p_greater < 0.05, na.rm = TRUE))
      min(n[p_greater < 0.05], na.rm = TRUE) else NA_integer_,
    first_below_n  = if (any(p_less < 0.05, na.rm = TRUE))
      min(n[p_less < 0.05], na.rm = TRUE) else NA_integer_,
    .groups        = "drop"
  )

# ----------------------------- PLOTS ---------------------------

plot_power_vs_n <- function(k_row_pick, k_col_pick) {
  df_k <- chi2_long_H1 %>%
    filter(k_row == k_row_pick, k_col == k_col_pick) %>%
    arrange(n)
  
  k_lab_here <- unique(df_k$k_lab)
  if (length(k_lab_here) == 0L) return(NULL)
  
  ggplot(df_k, aes(x = n, y = power, group = method, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous("n", breaks = sort(unique(df_k$n))) +
    scale_y_continuous("Power", limits = c(0, 1), labels = percent_format()) +
    labs(
      title    = glue("Chi-square – {k_lab_here} – Power vs n (H1)"),
      subtitle = "Permutation vs BOOTEI",
      color    = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_size_vs_n <- function(k_row_pick, k_col_pick) {
  df_k <- chi2_summary_all %>%
    filter(grepl("H0", scenario), k_row == k_row_pick, k_col == k_col_pick) %>%
    select(k_lab, n, rate_perm, rate_midp, rate_boot) %>%
    pivot_longer(c(rate_perm, rate_midp, rate_boot),
                 names_to = "method", values_to = "size") %>%
    mutate(
      method = recode(method, rate_perm = "Perm", rate_midp = "mid-p", rate_boot = "BOOTEI"),
      method = factor(method, levels = c("Perm", "mid-p", "BOOTEI"))
    ) %>%
    arrange(n)
  
  k_lab_here <- unique(df_k$k_lab)
  if (length(k_lab_here) == 0L) return(NULL)
  
  ggplot(df_k, aes(x = n, y = size, group = method, color = method)) +
    geom_hline(yintercept = alpha_nominal, linetype = "dashed") +
    geom_line() +
    geom_point() +
    scale_x_continuous("n", breaks = sort(unique(df_k$n))) +
    scale_y_continuous(
      "Type I error",
      limits = c(0, max(0.25, alpha_nominal * 2)),
      labels = percent_format()
    ) +
    labs(
      title    = glue("Chi-square – {k_lab_here} – Type I error vs n (H0)"),
      subtitle = glue("Dashed line = nominal alpha ({alpha_nominal})"),
      color    = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pH0_2 <- plot_size_vs_n(2L, 2L); pH1_2 <- plot_power_vs_n(2L, 2L)
pH0_3 <- plot_size_vs_n(3L, 3L); pH1_3 <- plot_power_vs_n(3L, 3L)
pH0_4 <- plot_size_vs_n(4L, 4L); pH1_4 <- plot_power_vs_n(4L, 4L)

if (!is.null(pH0_2)) print(pH0_2); if (!is.null(pH1_2)) print(pH1_2)
if (!is.null(pH0_3)) print(pH0_3); if (!is.null(pH1_3)) print(pH1_3)
if (!is.null(pH0_4)) print(pH0_4); if (!is.null(pH1_4)) print(pH1_4)

# ------------------ COMBINED POWER FIGURE (LOG10 X) ------------

col_perm   <- "#D55E00"
col_bootei <- "#0072B2"

theme_clean <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      plot.title         = element_text(face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1)
    )
}

plot_power_log10 <- function(k_row_pick, k_col_pick) {
  df_k <- chi2_long_H1 %>%
    filter(k_row == k_row_pick, k_col == k_col_pick) %>%
    mutate(Method = factor(method, levels = c("Perm", "BOOTEI"),
                           labels = c("Permutation", "BOOTEI"))) %>%
    arrange(n)
  
  k_lab_here <- unique(df_k$k_lab)
  n_vals     <- sort(unique(df_k$n))
  
  ggplot(df_k, aes(x = n, y = power, colour = Method, group = Method)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.8) +
    scale_x_log10(expression(log[10](n)), breaks = n_vals, labels = n_vals) +
    scale_y_continuous("Power", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_colour_manual(values = c("Permutation" = col_perm, "BOOTEI" = col_bootei)) +
    labs(title = k_lab_here, subtitle = expression("Power under " * H[1]), colour = NULL) +
    coord_cartesian(expand = FALSE) +
    theme_clean()
}

fig_power <- plot_power_log10(2, 2) + plot_power_log10(3, 3) + plot_power_log10(4, 4) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

print(fig_power)

# --------------------------- OUTPUT ----------------------------

print(H1_mcnemar_summary)
print(H1_mcnemar_midp_summary)
print(H0_summary)
print(H0_summary_noBH)
