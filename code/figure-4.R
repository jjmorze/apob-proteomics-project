# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-26
# Last modified: 2025-12-07
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, forcats, ggplot2, grid, gridExtra, cowplot,
#                      scales, tidyr, purrr, patchwork, ggplotify
#
# Goal:
#   Generate a multi-panel figure combining:
#   Panel A: Forest plots showing protein-fraction associations (UKBB & MESA)
#            - Simplified version without shape differentiation
#   Panel B: Spearman correlation heatmaps (MPS vs conventional biomarkers)
#   Panel C: Forest plot of MPS hazard ratios across models
#
# Required input data (in global environment):
#   1) Panel A:
#      - combined_ukbb_mesa_coxph.rda
#   2) Panel B:
#      - dataset.rda
#   3) Panel C:
#      - results_mps_long.rda
#
# Outputs:
#   - figure-4.png: Combined multi-panel figure
#     (11.8 x 5.4 inches, 1200 dpi)
#
# ==============================================================================


# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(scales)
  library(tidyr)
  library(purrr)
  library(patchwork)
})


# ==============================================================================
# PANEL A: FOREST PLOTS (UKBB & MESA, 3 PANES)
# ==============================================================================

# ------------------------------------------------------------------------------
# Data Preparation
# ------------------------------------------------------------------------------

df0 <- combined_ukbb_mesa_coxph %>%
  mutate(fraction = as.factor(trimws(as.character(fraction))))

# Select top 10 proteins by UKBB p-value within each fraction
top10_by_fraction <- df0 %>%
  group_by(fraction) %>%
  arrange(p.value.protein, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()


# Create long plotting table: one row per (fraction, gene, Study)
plot_df <- bind_rows(
  # UKBB rows
  top10_by_fraction %>%
    transmute(
      fraction,
      gene,
      Study    = "UKBB",
      estimate = estimate.ukbb.protein,
      se       = std.ukbb.protein,
      sign_ukbb = direction.ukbb
    ),
  # MESA rows
  top10_by_fraction %>%
    transmute(
      fraction,
      gene,
      Study    = "MESA",
      estimate = estimate.mesa,
      se       = std.error.mesa
    )
) %>%
  mutate(
    # Calculate HR and confidence intervals
    HR     = exp(estimate),
    LB     = exp(estimate - 1.96 * se),
    UB     = exp(estimate + 1.96 * se),
    log_HR = log(HR),
    log_LB = log(LB),
    log_UB = log(UB)
  )

# ------------------------------------------------------------------------------
# Global HR Limits & Aesthetics
# ------------------------------------------------------------------------------

hr_min <- floor(min(plot_df$LB, na.rm = TRUE) * 100) / 100
hr_max <- ceiling(max(plot_df$UB, na.rm = TRUE) * 100) / 100
hr_min <- max(hr_min, 0.2)
hr_max <- min(hr_max, 5)

candidate_breaks <- c(0.25, 0.33, 0.5, 0.67, 0.75, 1, 1.25, 1.5, 2, 3, 4, 5)
hr_breaks <- candidate_breaks[candidate_breaks >= hr_min & candidate_breaks <= hr_max]
log_breaks <- log(hr_breaks)
lab_breaks <- sprintf("%.2g", hr_breaks)

# Color mapping
study_colors <- c(UKBB = "steelblue", MESA = "firebrick")

# ------------------------------------------------------------------------------
# Helper Function: Build One Forest Panel
# ------------------------------------------------------------------------------

make_panel <- function(frac_label) {
  df_panel <- plot_df %>% filter(fraction == frac_label)
  
  # Per-fraction gene ordering: top-10 by UKBB p (top at top)
  levels_this <- top10_by_fraction %>%
    filter(fraction == frac_label) %>%
    arrange(p.value.protein) %>%
    pull(gene) %>%
    rev()
  
  df_panel <- df_panel %>%
    mutate(gene_fac = forcats::fct_relevel(factor(gene), levels_this))
  
  # Background stripes aligned to factor order
  stripe_df <- df_panel %>%
    distinct(gene_fac) %>%
    arrange(gene_fac) %>%
    mutate(
      y_id    = as.integer(gene_fac),
      ymin    = y_id - 0.5,
      ymax    = y_id + 0.5,
      fill_bg = if_else(y_id %% 2 == 0, "grey92", "white")
    )
  
  pd <- position_dodge(width = 0.55)
  
  p <- ggplot() +
    # Background stripes
    geom_rect(
      data = stripe_df %>% filter(fill_bg == "white"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      inherit.aes = FALSE,
      fill = "white",
      color = NA
    ) +
    geom_rect(
      data = stripe_df %>% filter(fill_bg == "grey92"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      inherit.aes = FALSE,
      fill = "grey92",
      color = NA
    ) +
    # Reference line at HR = 1
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    # Confidence interval bars
    geom_errorbarh(
      data = df_panel,
      aes(y = gene_fac, xmin = log_LB, xmax = log_UB, color = Study),
      height = 0.3,
      linewidth = 0.6,
      position = pd
    ) +
    # Point estimates (filled squares for both studies)
    geom_point(
      data = df_panel,
      aes(y = gene_fac, x = log_HR, color = Study),
      shape = 15,
      size = 2,
      stroke = 0.6,
      position = pd
    ) +
    # Scale settings
    scale_x_continuous(
      limits = log(c(hr_min - 0.02, hr_max + 0.02)),
      breaks = log_breaks,
      labels = lab_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(
      values = study_colors,
      name = "Study"
    ) +
    # Labels and theme
    labs(
      x = "Hazard ratio (per 1-SD)",
      y = NULL,
      title = as.character(frac_label)
    ) +
    theme_linedraw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title.x       = element_text(size = 11),
      legend.position    = "none",
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      plot.margin        = margin(4, 6, 4, 6),
      panel.background   = element_rect(fill = NA, color = NA),
      plot.background    = element_rect(fill = NA, color = NA)
    )
  
  ggplotGrob(p)
}

# ------------------------------------------------------------------------------
# Build Panels with TRL in Middle
# ------------------------------------------------------------------------------

fractions_all <- levels(df0$fraction)

# Put TRL-containing fraction in the middle
trl_idx <- which(grepl("trl", fractions_all, ignore.case = TRUE))
if (length(trl_idx) == 1 && length(fractions_all) == 3) {
  other <- setdiff(seq_along(fractions_all), trl_idx)
  pane_order <- c(other[1], trl_idx, other[2])
} else {
  pane_order <- seq_along(fractions_all)
}

panes <- lapply(fractions_all[pane_order], make_panel)

# ------------------------------------------------------------------------------
# Create Shared Legend
# ------------------------------------------------------------------------------

legend_dummy <- ggplot(
  plot_df,
  aes(x = log_HR, y = factor(gene), color = Study)
) +
  geom_point() +
  scale_color_manual(values = study_colors, name = "Study") +
  guides(color = guide_legend(order = 1)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    legend.box       = "vertical",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10)
  )

legend_grob <- cowplot::get_legend(legend_dummy)

# ------------------------------------------------------------------------------
# Arrange Panel A: 3 Panels + Legend
# ------------------------------------------------------------------------------

three_col <- arrangeGrob(grobs = panes, ncol = 3)
panelA_grob <- arrangeGrob(
  three_col,
  legend_grob,
  ncol   = 2,
  widths = c(1, 0.15)
)


# ==============================================================================
# PANEL B: SPEARMAN CORRELATION HEATMAPS
# ==============================================================================

# ------------------------------------------------------------------------------
# Define Variables and Labels
# ------------------------------------------------------------------------------

row_vars <- c("LDL_ridge_full", "TRL_ridge_full", "LPA_ridge_full")
panelA_cols <- c("LDL_ridge_full", "TRL_ridge_full", "LPA_ridge_full")
panelB_cols <- c(
  "conven_apob_imp",
  "conven_ldlp",
  "conven_ldlp_pct",
  "conven_trlp",
  "conven_trlp_pct",
  "conven_lpa_corr",
  "conven_lpa_corr_pct"
)

pretty_row <- c(
  LDL_ridge_full = "LDL MPS",
  TRL_ridge_full = "TRL MPS",
  LPA_ridge_full = "Lp(a) MPS"
)

pretty_col <- c(
  LDL_ridge_full      = "LDL MPS",
  TRL_ridge_full      = "TRL MPS",
  LPA_ridge_full      = "Lp(a) MPS",
  conven_apob_imp     = "apoB",
  conven_ldlp         = "LDL",
  conven_ldlp_pct     = "LDL (%apoB)",
  conven_trlp         = "TRL",
  conven_trlp_pct     = "TRL (%apoB)",
  conven_lpa_corr     = "Lp(a)",
  conven_lpa_corr_pct = "Lp(a) (%apoB)"
)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

as_numeric_safely <- function(v) {
  if (is.null(v)) return(NULL)
  if (!is.numeric(v)) suppressWarnings(v <- as.numeric(v))
  v
}

safe_spearman <- function(x, y) {
  x <- as_numeric_safely(x)
  y <- as_numeric_safely(y)
  if (is.null(x) || is.null(y)) return(NA_real_)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

get_corr_long <- function(dat, row_vars, col_vars) {
  rows_present <- intersect(row_vars, names(dat))
  cols_present <- intersect(col_vars, names(dat))
  if (length(rows_present) == 0 || length(cols_present) == 0) {
    return(tibble())
  }
  
  tidyr::expand_grid(row_var = rows_present, col_var = cols_present) %>%
    mutate(
      rho = purrr::map2_dbl(
        row_var,
        col_var,
        ~ safe_spearman(dat[[.x]], dat[[.y]])
      ),
      row_lab = pretty_row[row_var],
      col_lab = pretty_col[col_var],
      lab = sprintf("%.2f", rho)
    )
}

# ------------------------------------------------------------------------------
# Main Plotting Function
# ------------------------------------------------------------------------------

plot_spearman_heatmap_merged <- function(df) {
  longA <- get_corr_long(df, row_vars, panelA_cols)
  longB <- get_corr_long(df, row_vars, panelB_cols)
  
  if (nrow(longA) == 0 && nrow(longB) == 0) {
    stop("None of the required columns were found in the data.")
  }
  
  # Define ordering
  row_order   <- c("LDL MPS", "TRL MPS", "Lp(a) MPS")
  col_order_A <- c("LDL MPS", "TRL MPS", "Lp(a) MPS")
  col_order_B <- c(
    "apoB",
    "LDL",
    "LDL (%apoB)",
    "TRL",
    "TRL (%apoB)",
    "Lp(a)",
    "Lp(a) (%apoB)"
  )
  
  fix_orders <- function(d, col_levels) {
    d %>%
      mutate(
        row_lab = factor(row_lab, levels = row_order),
        col_lab = factor(col_lab, levels = col_levels)
      )
  }
  
  longA <- fix_orders(longA, col_order_A)
  longB <- fix_orders(longB, col_order_B)
  
  # Base tile function
  base_tile <- function(data, show_legend = FALSE) {
    ggplot(data, aes(x = col_lab, y = row_lab, fill = rho)) +
      geom_tile(color = "grey90", linewidth = 0.3) +
      geom_text(aes(label = ifelse(is.na(rho), "", lab)), size = 3) +
      scale_fill_gradient2(
        low = "steelblue",
        mid = "white",
        high = "firebrick",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Correlation",
        guide = if (show_legend) "colourbar" else "none"
      ) +
      labs(x = NULL, y = NULL) +
      theme_linedraw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()
      )
  }
  
  # Create subplots
  pA <- base_tile(longA, show_legend = FALSE)
  pB <- base_tile(longB, show_legend = TRUE) +
    theme(axis.text.y = element_blank())
  
  # Combine with patchwork
  pA + patchwork::plot_spacer() + pB +
    patchwork::plot_layout(widths = c(0.8, 0.05, 1.8))
}


# ==============================================================================
# PANEL C: MPS FOREST PLOT
# ==============================================================================

make_forest_plot_mps <- function(results_mps_long, return_legend = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(cowplot)
  })
  
  # Filter to target models and exposures
  target_models <- c(
    "Model 1: Clinical covariates",
    "Model 2: Clinical cov. + conventional"
  )
  target_expos <- c("LDL_ridge_full", "TRL_ridge_full", "LPA_ridge_full")
  
  df <- results_mps_long %>%
    filter(model_name %in% target_models, Exposure %in% target_expos) %>%
    mutate(
      ExposureLabel = case_when(
        Exposure == "LDL_ridge_full" ~ "LDL MPS",
        Exposure == "TRL_ridge_full" ~ "TRL MPS",
        Exposure == "LPA_ridge_full" ~ "Lp(a) MPS",
        TRUE ~ ExposureLabel
      ),
      # Top-to-bottom: LDL → TRL → Lp(a)
      ExposureLabel = factor(
        ExposureLabel,
        levels = c("Lp(a) MPS", "TRL MPS", "LDL MPS")
      ),
      model_short = recode(
        model_name,
        "Model 1: Clinical covariates" = "M1: Clinical\n covariates",
        "Model 2: Clinical cov. + conventional" = "M2: M1 +\n conventional"
      ),
      # M1 above, M2 below in legend
      model_short = factor(
        model_short,
        levels = c("M2: M1 +\n conventional", "M1: Clinical\n covariates")
      ),
      log_HR = log(HR),
      log_LB = log(HR_LCL),
      log_UB = log(HR_UCL)
    ) %>%
    tidyr::drop_na(log_HR, log_LB, log_UB)
  
  # Background stripes
  stripe_df <- df %>%
    distinct(ExposureLabel) %>%
    arrange(desc(ExposureLabel)) %>%
    mutate(
      ymin = as.integer(ExposureLabel) - 0.5,
      ymax = as.integer(ExposureLabel) + 0.5,
      fill_bg = if_else(row_number() %% 2 == 0, "grey92", "white")
    )
  
  # HR scale settings
  hr_min <- floor(min(df$HR_LCL, na.rm = TRUE) * 100) / 100
  hr_max <- ceiling(max(df$HR_UCL, na.rm = TRUE) * 100) / 100
  hr_min <- max(hr_min, 0.2)
  hr_max <- min(hr_max, 5)
  
  candidate_breaks <- c(
    0.25, 0.33, 0.5, 0.67, 0.75, 0.8, 0.9, 1,
    1.1, 1.25, 1.4, 1.5, 2, 3, 4, 5
  )
  hr_breaks <- candidate_breaks[
    candidate_breaks >= hr_min & candidate_breaks <= hr_max
  ]
  log_breaks <- log(hr_breaks)
  lab_breaks <- sprintf("%.2g", hr_breaks)
  
  # Color scheme
  cols <- c(
    "M1: Clinical\n covariates" = "#6EA6CD",
    "M2: M1 +\n conventional" = "#2B5C8A"
  )
  
  pd <- position_dodge(width = 0.28)
  
  base_plot <- ggplot() +
    geom_rect(
      data = stripe_df,
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = fill_bg),
      inherit.aes = FALSE,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
      data = df,
      aes(
        y = ExposureLabel,
        xmin = log_LB,
        xmax = log_UB,
        color = model_short,
        group = model_short
      ),
      height = 0.3,
      linewidth = 0.7,
      position = pd
    ) +
    geom_point(
      data = df,
      aes(
        y = ExposureLabel,
        x = log_HR,
        color = model_short,
        group = model_short
      ),
      shape = 15,
      size = 2.3,
      stroke = 0.6,
      position = pd
    ) +
    scale_x_continuous(
      limits = log(c(hr_min - 0.02, hr_max + 0.02)),
      breaks = log_breaks,
      labels = lab_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(
      values = cols,
      breaks = c("M1: Clinical\n covariates", "M2: M1 +\n conventional"),
      name = "Model"
    ) +
    labs(x = "Hazard ratio (per 1-SD)", y = NULL) +
    theme_linedraw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title.x       = element_text(size = 11),
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.title       = element_text(size = 11, face = "plain"),
      legend.text        = element_text(size = 10),
      plot.margin        = margin(4, 18, 4, 18),
      panel.background   = element_rect(fill = NA, color = NA),
      plot.background    = element_rect(fill = NA, color = NA)
    )
  
  if (return_legend) {
    cowplot::get_legend(base_plot)
  } else {
    base_plot + theme(legend.position = "none")
  }
}


# ==============================================================================
# COMBINE ALL PANELS
# ==============================================================================

# ------------------------------------------------------------------------------
# Generate Individual Panels
# ------------------------------------------------------------------------------

# Panel B plot
pB <- plot_spearman_heatmap_merged(dataset)

# Panel C plot and legend
pC_no_legend <- make_forest_plot_mps(results_mps_long, return_legend = FALSE)
pC_legend    <- make_forest_plot_mps(results_mps_long, return_legend = TRUE)

# ------------------------------------------------------------------------------
# Add Corner Labels
# ------------------------------------------------------------------------------

label_plot <- function(plot_or_grob, label_text, size = 14) {
  is_gg <- inherits(plot_or_grob, "ggplot")
  
  if (is_gg) {
    p <- cowplot::ggdraw(plot_or_grob) +
      cowplot::draw_plot_label(
        label = label_text,
        x = 0.01,
        y = 0.99,
        hjust = 0,
        vjust = 1,
        size = size,
        fontface = "plain"
      )
    return(ggplotGrob(p))
  } else {
    # Handle grobs
    p <- cowplot::ggdraw() +
      cowplot::draw_plot(plot_or_grob, x = 0, y = 0, width = 1, height = 1) +
      cowplot::draw_plot_label(
        label = label_text,
        x = 0.01,
        y = 0.99,
        hjust = 0,
        vjust = 1,
        size = size,
        fontface = "plain"
      )
    return(ggplotGrob(p))
  }
}

gA <- label_plot(panelA_grob,  "A", size = 14)
gB <- label_plot(pB,           "B", size = 14)
gC <- label_plot(pC_no_legend, "C", size = 14)

# ------------------------------------------------------------------------------
# Create Layout Matrix
# ------------------------------------------------------------------------------

# Row 1: Panel A spans all 3 columns
# Row 2: Panel B | Panel C | Panel C legend
lay <- rbind(
  c(1, 1, 1),
  c(2, 3, 4)
)

combined_grob <- arrangeGrob(
  grobs = list(gA, gB, gC, pC_legend),
  layout_matrix = lay,
  widths  = c(1.3, 0.6, 0.24),
  heights = c(1.10, 1.00)
)

# ------------------------------------------------------------------------------
# Display Combined Figure
# ------------------------------------------------------------------------------

grid::grid.newpage()
grid::grid.draw(combined_grob)

# ------------------------------------------------------------------------------
# Save Figure
# ------------------------------------------------------------------------------

ggsave(
  "figure-4.png",
  ggplotify::as.ggplot(combined_grob),
  width = 11.8,
  height = 5.4,
  dpi = 1200,
  bg = "white"
)
