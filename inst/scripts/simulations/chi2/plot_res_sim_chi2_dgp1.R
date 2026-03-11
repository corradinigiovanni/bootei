# ================================================================
# BOOTEI – Pearson chi-square: analysis of chi2_DETAIL_GLOBAL_dgp1.csv
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
  library(cowplot)
})

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")

# ---------------------------- PATH -----------------------------

ROOT <- "."
file_global <- file.path(ROOT, "chi2_DETAIL_GLOBAL_dgp1.csv")
stopifnot("chi2_DETAIL_GLOBAL_dgp1.csv not found!" = file.exists(file_global))

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

    # --------------------------- OUTPUT ----------------------------

print(H1_mcnemar_summary)
print(H1_mcnemar_midp_summary)
print(H0_summary)
print(H0_summary_noBH)

# ----------------------------- PLOTS ---------------------------

# ------------------ Theme ------------------
theme_jrssb_fig_shared_axes <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line  = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = grid::unit(0.14, "cm"),
      
      axis.text.x = element_text(
        size  = base_size - 2,
        angle = 35,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = base_size - 1),
      
      # shared axis titles outside the panels
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      plot.title    = element_blank(),
      plot.subtitle = element_text(
        size = base_size,
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      
      # move panel tags slightly higher
      plot.tag = element_text(size = base_size + 1, face = "bold"),
      plot.tag.position = c(0.01, 1.03),
      
      legend.position = "none",
      
      # a bit more top margin so A/B/C have room
      plot.margin = margin(10, 4, 5, 4)
    )
}

# ------------------ Single-panel function ------------------
make_powersize_panel_chi2_dgp1 <- function(k_row_pick, k_col_pick, panel_tag, base_size = 11) {
  
  # ---------- Power under H1 ----------
  df_H1 <- chi2_summary_all %>%
    filter(grepl("H1", scenario), k_row == k_row_pick, k_col == k_col_pick) %>%
    select(k_lab, n, rate_perm, rate_boot) %>%
    pivot_longer(
      cols = c(rate_perm, rate_boot),
      names_to = "method",
      values_to = "rate"
    ) %>%
    mutate(
      Method = recode(
        method,
        rate_perm = "CLASSIC",
        rate_boot = "BOOTEI"
      ),
      Metric = "Power"
    )
  
  # ---------- Size under H0 ----------
  df_H0 <- chi2_summary_all %>%
    filter(grepl("H0", scenario), k_row == k_row_pick, k_col == k_col_pick) %>%
    select(k_lab, n, rate_perm, rate_boot) %>%
    pivot_longer(
      cols = c(rate_perm, rate_boot),
      names_to = "method",
      values_to = "rate"
    ) %>%
    mutate(
      Method = recode(
        method,
        rate_perm = "CLASSIC",
        rate_boot = "BOOTEI"
      ),
      Metric = "Size"
    )
  
  df_plot <- bind_rows(df_H1, df_H0) %>%
    mutate(
      Method = factor(Method, levels = c("CLASSIC", "BOOTEI")),
      Metric = factor(Metric, levels = c("Power", "Size"))
    ) %>%
    arrange(Method, Metric, n)
  
  n_vals <- sort(unique(df_plot$n))
  panel_label <- paste0(k_row_pick, "\u00D7", k_col_pick)
  
  # Fewer x ticks only in panel A
  x_breaks <- if (k_row_pick == 2L && k_col_pick == 2L) {
    c(10, 20, 50, 100, 500, 1000)
  } else {
    n_vals
  }
  
  ggplot(
    df_plot,
    aes(
      x = n,
      y = rate,
      group = interaction(Method, Metric)
    )
  ) +
    geom_line(
      aes(colour = Method, linetype = Method),
      linewidth = 0.75
    ) +
    geom_point(
      aes(shape = Metric, colour = Method),
      size = 1.45,
      stroke = 0.60,
      fill = "white"
    ) +
    geom_hline(
      yintercept = alpha_nominal,
      colour = "black",
      linewidth = 0.75,
      linetype = "dotted",
      lineend = "round"
    ) +
    scale_colour_manual(
      values = c(
        "CLASSIC" = "grey35",
        "BOOTEI"  = "black"
      ),
      guide = "none"
    ) +
    scale_linetype_manual(
      values = c(
        "CLASSIC" = "11",
        "BOOTEI"  = "solid"
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "Power" = 21,
        "Size"  = 22
      ),
      guide = "none"
    ) +
    scale_x_log10(
      breaks = x_breaks,
      labels = label_number()(x_breaks),
      limits = range(n_vals)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, alpha_nominal, 0.25, 0.50, 0.75, 1.00),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      subtitle = panel_label,
      tag = panel_tag
    ) +
    theme_jrssb_fig_shared_axes(base_size = base_size)
}

# ------------------ Build panels ------------------
powersize_panel_A_chi2_dgp1 <- make_powersize_panel_chi2_dgp1(
  k_row_pick = 2L, k_col_pick = 2L, panel_tag = "A", base_size = 11
)

powersize_panel_B_chi2_dgp1 <- make_powersize_panel_chi2_dgp1(
  k_row_pick = 3L, k_col_pick = 3L, panel_tag = "B", base_size = 11
)

powersize_panel_C_chi2_dgp1 <- make_powersize_panel_chi2_dgp1(
  k_row_pick = 4L, k_col_pick = 4L, panel_tag = "C", base_size = 11
)

# ------------------ Panels row ------------------
panels_row_chi2_dgp1 <-
  powersize_panel_A_chi2_dgp1 +
  powersize_panel_B_chi2_dgp1 +
  powersize_panel_C_chi2_dgp1 +
  plot_layout(
    nrow   = 1,
    widths = c(1.28, 1, 1)
  )

# ------------------ Final combined plot with shared axis titles ------------------
powersize_plot_chi2_dgp1 <-
  ggdraw() +
  draw_plot(
    panels_row_chi2_dgp1,
    x = 0.085, y = 0.12,
    width = 0.90, height = 0.82
  ) +
  draw_label(
    "Empirical rejection rate",
    x = 0.025, y = 0.53,
    angle = 90,
    size = 11
  ) +
  draw_label(
    "Sample size, n (log10 scale)",
    x = 0.54, y = 0.035,
    size = 11
  )

print(powersize_plot_chi2_dgp1)

# ------------------ Save outputs ------------------
ggsave(
  filename = "/media/data/corradini/BOOTEI/sim_power_size_chi2/powersize_plot_chi2_dgp1.pdf",
  plot     = powersize_plot_chi2_dgp1,
  width    = 7,
  height   = 6,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "/media/data/corradini/BOOTEI/sim_power_size_chi2/powersize_plot_chi2_dgp1.tiff",
  plot     = powersize_plot_chi2_dgp1,
  width    = 7,
  height   = 6,
  units    = "in",
  dpi      = 600,
  compression = "lzw",
  device   = "tiff"
)
