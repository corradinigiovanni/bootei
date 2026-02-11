# ================================================================
# BOOTEI – Welch t-test: analysis of welch_DETAIL_GLOBAL_dgp2.csv
#   - Power (H1): Perm vs BOOTEI, by n_eff = n1 + n2 and k
#   - McNemar tests under H1: BOOTEI > Perm; BOOTEI vs mid-p
#   - Size under H0: binom.test vs alpha (BH-adjusted and raw)
#   - Summary power/size figures with log10(n_eff) x-axis
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

ROOT <- "_"
file_global <- file.path(ROOT, "welch_DETAIL_GLOBAL_dgp2.csv")
stopifnot("welch_DETAIL_GLOBAL_dgp2.csv not found!" = file.exists(file_global))

# --------------------- LOAD + CELL SUMMARIES -------------------

welch_detail_all <- read_csv(file_global, show_col_types = FALSE) %>%
  mutate(
    k      = as.integer(k),
    n1     = as.integer(n1),
    n2     = as.integer(n2),
    n_eff  = as.integer(n1 + n2),
    n_pair = glue("{n1}+{n2}"),
    k_lab  = glue("k = {k}")
  )

welch_summary_all <- welch_detail_all %>%
  group_by(scenario, k, k_lab, n1, n2, n_eff, n_pair, alpha, assoc) %>%
  summarise(
    total_sim = n(),
    rej_perm  = sum(sig_perm),
    rej_midp  = sum(sig_midp),
    rej_boot  = sum(sig_boot),
    rate_perm = mean(sig_perm),
    rate_midp = mean(sig_midp),
    rate_boot = mean(sig_boot),
    .groups   = "drop"
  )

alpha_nominal <- welch_summary_all$alpha %>%
  unique() %>%
  sort() %>%
  tail(1)

cat(glue("Nominal alpha: {alpha_nominal}\n"))

# ------------------- LONG FORMAT (POWER, H1) -------------------

welch_long_H1 <- welch_summary_all %>%
  filter(grepl("H1", scenario)) %>%
  select(scenario, k, k_lab, n1, n2, n_eff, n_pair, assoc, alpha,
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

welch_H1_mcnemar <- welch_detail_all %>%
  filter(grepl("H1", scenario)) %>%
  group_by(k, k_lab, n1, n2, n_eff, n_pair, assoc) %>%
  summarise(
    total_sim   = n(),
    power_perm  = mean(sig_perm),
    power_boot  = mean(sig_boot),
    a           = sum(sig_boot &  sig_perm),
    b           = sum(sig_boot & !sig_perm),
    c           = sum(!sig_boot &  sig_perm),
    d           = sum(!sig_boot & !sig_perm),
    .groups     = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_mcn        = mcnemar_binom(b, c),  # BOOTEI > Perm
    diff_power   = power_boot - power_perm,
    sig_mcn_0.05 = is.finite(p_mcn) & p_mcn < 0.05
  ) %>%
  ungroup() %>%
  arrange(k, n_eff, n1, n2)

welch_H1_mcnemar_summary <- welch_H1_mcnemar %>%
  group_by(k, k_lab) %>%
  summarise(
    n_cells        = n(),
    n_sig_raw      = sum(sig_mcn_0.05, na.rm = TRUE),
    first_sig_neff = if (any(sig_mcn_0.05, na.rm = TRUE))
      min(n_eff[sig_mcn_0.05], na.rm = TRUE) else NA_integer_,
    last_sig_neff  = if (any(sig_mcn_0.05, na.rm = TRUE))
      max(n_eff[sig_mcn_0.05], na.rm = TRUE) else NA_integer_,
    max_diff_power = if (n_cells > 0) max(diff_power, na.rm = TRUE) else NA_real_,
    .groups        = "drop"
  )

welch_H1_mcnemar_midp <- welch_detail_all %>%
  filter(grepl("H1", scenario)) %>%
  group_by(k, k_lab, n1, n2, n_eff, n_pair, assoc) %>%
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
    p_boot_gt_midp   = mcnemar_binom(b, c),
    p_midp_gt_boot   = mcnemar_binom(c, b),
    diff_power       = power_boot - power_midp,
    sig_boot_gt_midp = is.finite(p_boot_gt_midp) & p_boot_gt_midp < 0.05,
    sig_midp_gt_boot = is.finite(p_midp_gt_boot) & p_midp_gt_boot < 0.05
  ) %>%
  ungroup() %>%
  arrange(k, n_eff, n1, n2)

welch_H1_mcnemar_midp_summary <- welch_H1_mcnemar_midp %>%
  group_by(k, k_lab) %>%
  summarise(
    n_cells                 = n(),
    n_sig_boot_gt_midp      = sum(sig_boot_gt_midp, na.rm = TRUE),
    n_sig_midp_gt_boot      = sum(sig_midp_gt_boot, na.rm = TRUE),
    first_neff_boot_gt_midp = if (any(sig_boot_gt_midp, na.rm = TRUE))
      min(n_eff[sig_boot_gt_midp], na.rm = TRUE) else NA_integer_,
    first_neff_midp_gt_boot = if (any(sig_midp_gt_boot, na.rm = TRUE))
      min(n_eff[sig_midp_gt_boot], na.rm = TRUE) else NA_integer_,
    max_diff_power          = if (n_cells > 0) max(diff_power, na.rm = TRUE) else NA_real_,
    .groups                 = "drop"
  )

# -------------------- SIZE TESTS (H0): BH + RAW ----------------

welch_h0_base <- welch_summary_all %>%
  filter(grepl("H0", scenario)) %>%
  transmute(
    k      = as.integer(k),
    k_lab,
    n1     = as.integer(n1),
    n2     = as.integer(n2),
    n_eff  = as.integer(n_eff),
    n_pair,
    alpha,
    total_sim,
    rej_perm,
    rej_midp,
    rej_boot,
    rate_perm,
    rate_midp,
    rate_boot
  )

welch_h0_long <- welch_h0_base %>%
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

welch_H0_tests <- welch_h0_long %>%
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
  arrange(k, n_eff, n1, n2, method)

welch_H0_summary <- welch_H0_tests %>%
  group_by(method, k, k_lab) %>%
  summarise(
    n_cells          = n(),
    n_above_BH       = sum(p_greater_BH < 0.05, na.rm = TRUE),
    n_below_BH       = sum(p_less_BH    < 0.05, na.rm = TRUE),
    first_above_neff = if (any(p_greater_BH < 0.05, na.rm = TRUE))
      min(n_eff[p_greater_BH < 0.05], na.rm = TRUE) else NA_integer_,
    .groups          = "drop"
  )

welch_H0_tests_noBH <- welch_H0_tests %>%
  mutate(
    size_flag_noBH = case_when(
      p_greater < 0.05 ~ "↑ (size > α, raw)",
      p_less    < 0.05 ~ "↓ (size < α, raw)",
      TRUE             ~ "OK"
    )
  ) %>%
  arrange(k, n_eff, n1, n2, method)

welch_H0_summary_noBH <- welch_H0_tests_noBH %>%
  group_by(method, k, k_lab) %>%
  summarise(
    n_cells          = n(),
    n_above_raw      = sum(p_greater < 0.05, na.rm = TRUE),
    n_below_raw      = sum(p_less    < 0.05, na.rm = TRUE),
    first_above_neff = if (any(p_greater < 0.05, na.rm = TRUE))
      min(n_eff[p_greater < 0.05], na.rm = TRUE) else NA_integer_,
    first_below_neff = if (any(p_less < 0.05, na.rm = TRUE))
      min(n_eff[p_less < 0.05], na.rm = TRUE) else NA_integer_,
    .groups          = "drop"
  )

# ----------------------------- PLOTS ---------------------------

plot_power_vs_neff <- function(k_pick) {
  df_k <- welch_long_H1 %>%
    filter(k == k_pick) %>%
    arrange(n_eff)
  
  k_lab_here <- unique(df_k$k_lab)
  if (length(k_lab_here) == 0L) return(NULL)
  
  ggplot(df_k, aes(x = n_eff, y = power, group = method, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous("n_eff = n1 + n2", breaks = sort(unique(df_k$n_eff))) +
    scale_y_continuous("Power", limits = c(0, 1), labels = percent_format()) +
    labs(
      title    = glue("Welch t-test (DGP2) – {k_lab_here} – Power vs n_eff"),
      subtitle = "Permutation vs BOOTEI",
      color    = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

welch_long_H0 <- welch_summary_all %>%
  filter(grepl("H0", scenario)) %>%
  select(k, k_lab, n1, n2, n_eff, alpha, rate_perm, rate_midp, rate_boot) %>%
  pivot_longer(c(rate_perm, rate_midp, rate_boot),
               names_to = "method", values_to = "size") %>%
  mutate(
    method = recode(method, rate_perm = "Perm", rate_midp = "mid-p", rate_boot = "BOOTEI"),
    method = factor(method, levels = c("Perm", "mid-p", "BOOTEI"))
  )

plot_size_vs_neff <- function(k_pick) {
  df_k <- welch_long_H0 %>%
    filter(k == k_pick) %>%
    arrange(n_eff)
  
  k_lab_here <- unique(df_k$k_lab)
  if (length(k_lab_here) == 0L) return(NULL)
  
  ggplot(df_k, aes(x = n_eff, y = size, group = method, color = method)) +
    geom_hline(yintercept = alpha_nominal, linetype = "dashed") +
    geom_line() +
    geom_point() +
    scale_x_continuous("n_eff = n1 + n2", breaks = sort(unique(df_k$n_eff))) +
    scale_y_continuous(
      "Type I error",
      limits = c(0, max(0.25, alpha_nominal * 2)),
      labels = percent_format()
    ) +
    labs(
      title    = glue("Welch t-test (DGP2) – {k_lab_here} – Type I error vs n_eff (H0)"),
      subtitle = glue("Dashed line = nominal alpha ({alpha_nominal})"),
      color    = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pH0_3 <- plot_size_vs_neff(3L); pH1_3 <- plot_power_vs_neff(3L)
pH0_5 <- plot_size_vs_neff(5L); pH1_5 <- plot_power_vs_neff(5L)
pH0_7 <- plot_size_vs_neff(7L); pH1_7 <- plot_power_vs_neff(7L)

if (!is.null(pH0_3)) print(pH0_3); if (!is.null(pH1_3)) print(pH1_3)
if (!is.null(pH0_5)) print(pH0_5); if (!is.null(pH1_5)) print(pH1_5)
if (!is.null(pH0_7)) print(pH0_7); if (!is.null(pH1_7)) print(pH1_7)

# ------------------ COMBINED POWER FIGURE (LOG10 X) ------------

col_perm   <- "#D55E00"
col_bootei <- "#0072B2"

theme_clean_big <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      
      axis.text.x  = element_text(size = base_size + 0.5, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = base_size + 0.5),
      
      axis.title.x = element_text(size = base_size + 4, margin = margin(t = 10)),
      axis.title.y = element_text(size = base_size + 4, margin = margin(r = 10)),
      
      plot.title    = element_text(face = "bold", size = base_size + 6),
      plot.subtitle = element_text(size = base_size + 4, margin = margin(b = 8)),
      
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.text     = element_text(size = base_size + 2),
      
      plot.margin = margin(10, 12, 10, 12)
    )
}

plot_power_log10 <- function(k_pick, base_size = 14) {
  df_k <- welch_long_H1 %>%
    filter(k == k_pick) %>%
    mutate(Method = factor(method,
                           levels = c("Perm", "BOOTEI"),
                           labels = c("CLASSIC", "BOOTEI"))) %>%
    arrange(n_eff)
  
  k_lab_here <- unique(df_k$k_lab)
  n_vals     <- sort(unique(df_k$n_eff))
  
  ggplot(df_k, aes(x = n_eff, y = power, colour = Method, group = Method)) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2.4) +
    scale_x_log10(expression(log[10](n[1] + n[2])), breaks = n_vals, labels = n_vals) +
    scale_y_continuous("Power", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_colour_manual(values = c("CLASSIC" = col_perm, "BOOTEI" = col_bootei)) +
    labs(
      title    = k_lab_here,
      subtitle = expression("Power under " * H[1]),
      colour   = NULL
    ) +
    coord_cartesian(expand = FALSE) +
    theme_clean_big(base_size = base_size)
}

fig_power <- plot_power_log10(3L, base_size = 14) +
  plot_power_log10(5L, base_size = 14) +
  plot_power_log10(7L, base_size = 14) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

print(fig_power)

# ------------------ COMBINED SIZE FIGURE (LOG10 X) H0 ----------

plot_size_log10 <- function(k_pick, base_size = 14) {
  
  df_k <- welch_summary_all %>%
    filter(grepl("H0", scenario), k == k_pick) %>%
    select(k, k_lab, n_eff, rate_perm, rate_boot) %>%
    pivot_longer(
      c(rate_perm, rate_boot),
      names_to  = "method",
      values_to = "size"
    ) %>%
    mutate(
      Method = recode(method,
                      rate_perm = "CLASSIC",
                      rate_boot = "BOOTEI"),
      Method = factor(Method, levels = c("CLASSIC", "BOOTEI"))
    ) %>%
    arrange(n_eff)
  
  k_lab_here <- unique(df_k$k_lab)
  n_vals     <- sort(unique(df_k$n_eff))
  
  ggplot(df_k, aes(x = n_eff, y = size, colour = Method, group = Method)) +
    geom_hline(yintercept = alpha_nominal, linetype = "dashed", linewidth = 0.9) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2.4) +
    scale_x_log10(expression(log[10](n[1] + n[2])), breaks = n_vals, labels = n_vals) +
    scale_y_continuous(
      "Type I error",
      labels = percent_format(accuracy = 1),
      limits = c(0, max(0.25, alpha_nominal * 2))
    ) +
    scale_colour_manual(values = c("CLASSIC" = col_perm, "BOOTEI" = col_bootei)) +
    labs(
      title    = k_lab_here,
      subtitle = bquote("Size under " * H[0]),
      colour   = NULL
    ) +
    coord_cartesian(expand = FALSE) +
    theme_clean_big(base_size = base_size) +
    theme(plot.caption = element_text(size = base_size + 2, hjust = 0, margin = margin(t = 6)))
}

fig_size <- plot_size_log10(3L, base_size = 14) +
  plot_size_log10(5L, base_size = 14) +
  plot_size_log10(7L, base_size = 14) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

print(fig_size)

# --------------------------- OUTPUT ----------------------------

print(welch_H1_mcnemar_summary)
print(welch_H1_mcnemar_midp_summary)
print(welch_H0_summary)
print(welch_H0_summary_noBH)


# ------------------ ONE FIGURE: POWER (H1) + SIZE (H0), LOG10 X ----------

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(glue)
library(patchwork)

col_perm   <- "#D55E00"
col_bootei <- "#0072B2"

theme_clean_big <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      
      axis.text.x  = element_text(size = base_size + 0.5, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = base_size + 0.5),
      
      axis.title.x = element_text(size = base_size + 4, margin = margin(t = 10)),
      axis.title.y = element_text(size = base_size + 4, margin = margin(r = 10)),
      
      axis.ticks.length = grid::unit(0.25, "cm"),
      axis.ticks        = element_line(linewidth = 0.8),
      
      plot.title = element_text(face = "bold", size = base_size + 6),
      
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.text     = element_text(size = base_size + 2),
      
      # extra left margin so the external tick/label are not clipped
      plot.margin = margin(10, 12, 10, 26)
    )
}

plot_size_power_log10 <- function(k_pick, base_size = 14) {
  
  # -------------------- Power under H1 --------------------
  df_H1 <- welch_summary_all %>%
    filter(grepl("H1", scenario), k == k_pick) %>%
    select(k_lab, n_eff, rate_perm, rate_boot) %>%
    pivot_longer(
      c(rate_perm, rate_boot),
      names_to  = "method",
      values_to = "rate"
    ) %>%
    mutate(
      Method = recode(method,
                      rate_perm = "CLASSIC",
                      rate_boot = "BOOTEI"),
      Method = factor(Method, levels = c("CLASSIC", "BOOTEI")),
      Metric = "Power"
    )
  
  # -------------------- Size under H0 ---------------------
  df_H0 <- welch_summary_all %>%
    filter(grepl("H0", scenario), k == k_pick) %>%
    select(k_lab, n_eff, rate_perm, rate_boot) %>%
    pivot_longer(
      c(rate_perm, rate_boot),
      names_to  = "method",
      values_to = "rate"
    ) %>%
    mutate(
      Method = recode(method,
                      rate_perm = "CLASSIC",
                      rate_boot = "BOOTEI"),
      Method = factor(Method, levels = c("CLASSIC", "BOOTEI")),
      Metric = "Size"
    )
  
  # -------------------- bind + labels ---------------------
  df_k <- bind_rows(df_H1, df_H0) %>%
    mutate(Metric = factor(Metric, levels = c("Power", "Size"))) %>%
    arrange(n_eff, Method, Metric)
  
  k_lab_here <- unique(df_k$k_lab)
  n_vals     <- sort(unique(df_k$n_eff))
  
  # -------------------- external alpha tick + label (same visual length across panels) ----
  # fixed fraction of the panel width in log10-scale => k=3 tick matches k=5/7 visually
  tick_frac      <- 0.085
  lab_gap_log10  <- 0.04
  
  x_min <- min(n_vals)
  x_max <- max(n_vals)
  
  log_range     <- log10(x_max) - log10(x_min)
  tick_w_log10  <- tick_frac * log_range
  
  x_axis_left <- x_min
  x_tick_out  <- x_axis_left / (10^tick_w_log10)
  x_lab_out   <- x_axis_left / (10^(tick_w_log10 + lab_gap_log10))
  
  alpha_lab <- paste0(round(100 * alpha_nominal), "%")  # typically "5%"
  
  ggplot(df_k,
         aes(x = n_eff, y = rate,
             colour = Method,
             group  = interaction(Method, Metric))) +
    
    geom_hline(yintercept = alpha_nominal, linetype = "dashed", linewidth = 0.9) +
    
    geom_line(linewidth = 1.3) +
    geom_point(aes(shape = Metric), size = 2.4) +
    scale_shape_manual(values = c(Power = 16, Size = 17), guide = "none") +
    
    # ---- external tick + label
    annotate(
      "segment",
      x = x_tick_out, xend = x_axis_left,
      y = alpha_nominal, yend = alpha_nominal,
      linewidth = 1.6,
      colour = "black"
    ) +
    annotate(
      "text",
      x = x_lab_out, y = alpha_nominal,
      label = alpha_lab,
      hjust = 1, vjust = 0.5,
      size  = (base_size + 2) / ggplot2::.pt,
      colour = "black"
    ) +
    
    scale_x_log10(
      expression(log[10](n[1] + n[2])),
      breaks = n_vals, labels = n_vals,
      limits = range(n_vals),
      oob    = scales::oob_keep
    ) +
    scale_y_continuous(
      "Rate",
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_colour_manual(values = c("CLASSIC" = col_perm, "BOOTEI" = col_bootei)) +
    labs(title = k_lab_here[1], colour = NULL) +
    
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_clean_big(base_size = base_size)
}

p3 <- plot_size_power_log10(3L, base_size = 14)
p5 <- plot_size_power_log10(5L, base_size = 14)
p7 <- plot_size_power_log10(7L, base_size = 14)

fig_size_power <- (p3 + p5 + p7) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Power under H1 and Size under H0",
    theme = theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 10))
    )
  ) &
  theme(legend.position = "bottom")

print(fig_size_power)
