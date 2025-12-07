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
#   Required packages: dplyr, tidyr, ggplot2, scales, purrr, cowplot
#
# Goal:
#   Create publication-quality visualizations of lipoprotein distributions:
#   1) Half-boxplot showing absolute and relative concentrations of 
#      LDL, TRL, and Lp(a) particles
#   2) Correlation heatmap showing Spearman correlations between
#      ApoB, LDL, TRL, and Lp(a) measurements
#   3) Combine both visualizations into a two-panel figure
#
# Required input data:
#   - dataset.rda
#
# Main outputs:
#   - sfigure-lipoprotein.png (combined visualization)
#
# ==============================================================================

# ==============================================================================
# LIBRARIES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)
library(cowplot)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Define column names for absolute concentrations (nmol/L)
nmol_cols <- c("conven_ldlp", "conven_trlp", "conven_lpa_corr")

# Define column names for relative concentrations (% of apoB)
pct_cols <- paste0(nmol_cols, "_pct")

# All variables including apoB
vars_abs <- c("conven_apob", nmol_cols)

# Display names for variables
nice_names <- c(
  conven_apob     = "ApoB",
  conven_ldlp     = "LDL",
  conven_trlp     = "TRL",
  conven_lpa_corr = "Lp(a)"
)

# Color palette for visualizations
my_cols <- c(
  "LDL"   = "#FF9E9DFF",
  "TRL"   = "#3FB8AFFF",
  "Lp(a)" = "#7FA7C7FF"
)

# ==============================================================================
# DATA PREPARATION: HALF-BOXPLOT
# ==============================================================================

# Prepare absolute concentration data (nmol/L)
df_nm <- dataset %>%
  select(all_of(nmol_cols)) %>%
  pivot_longer(everything(), names_to = "analyte", values_to = "nmol")

# Prepare relative concentration data (% of apoB)
df_pct <- dataset %>%
  select(all_of(pct_cols)) %>%
  pivot_longer(everything(), names_to = "analyte_pct", values_to = "pct") %>%
  mutate(
    analyte = sub("_pct$", "", analyte_pct),
    pct = pmin(pct, 100)
  ) %>%
  select(-analyte_pct)

# Convert percentage to 0-100 scale if necessary
if (max(df_pct$pct, na.rm = TRUE) <= 1) {
  df_pct <- mutate(df_pct, pct = pct * 100)
}

# Calculate scaling factor for secondary axis
max_nmol <- max(df_nm$nmol, na.rm = TRUE)
s_factor <- max_nmol / 100

# Scale percentage data to nmol scale for plotting
df_pct <- df_pct %>%
  mutate(nmol_scaled = pct * s_factor)

# Apply factor labels to analyte names
df_nm <- df_nm %>%
  mutate(analyte = factor(analyte, levels = nmol_cols, labels = nice_names[nmol_cols]))

df_pct <- df_pct %>%
  mutate(analyte = factor(analyte, levels = nmol_cols, labels = nice_names[nmol_cols]))

# ==============================================================================
# COMPUTE STATISTICS FOR BOXPLOTS
# ==============================================================================

# Function to compute key percentiles (1st, 25th, 50th, 75th, 99th)
compute_stats <- function(x) {
  qs <- quantile(x, probs = c(0.01, 0.25, 0.5, 0.75, 0.99), na.rm = TRUE)
  data.frame(
    p1  = qs[1],
    q1  = qs[2],
    med = qs[3],
    q3  = qs[4],
    p99 = qs[5]
  )
}

# Compute statistics for absolute concentrations
stats_nm <- df_nm %>%
  group_by(analyte) %>%
  summarise(compute_stats(nmol), .groups = "drop") %>%
  mutate(type = "nmol")

# Compute statistics for relative concentrations (scaled)
stats_pct <- df_pct %>%
  group_by(analyte) %>%
  summarise(compute_stats(nmol_scaled), .groups = "drop") %>%
  mutate(type = "pct")

# Combine statistics
stats_all <- bind_rows(stats_nm, stats_pct) %>%
  mutate(x = as.numeric(analyte) * 0.35)

# ==============================================================================
# DEFINE BOXPLOT GEOMETRY
# ==============================================================================

# Geometric parameters for half-boxplots
half_width     <- 0.12
whisker_offset <- 0.02
cap_width      <- 0.01
gap            <- 0

# Calculate coordinates for boxplot elements
stats_all <- stats_all %>%
  mutate(
    # Whisker x-position
    x_w = ifelse(type == "nmol", x - whisker_offset, x + whisker_offset),
    
    # Preliminary box boundaries
    prelim_xmin = ifelse(type == "nmol", x - half_width, x),
    prelim_xmax = ifelse(type == "nmol", x, x + half_width),
    
    # Final box boundaries with gap adjustment
    xmin = ifelse(type == "nmol", prelim_xmin, pmax(prelim_xmin, x_w + gap)),
    xmax = ifelse(type == "nmol", pmin(prelim_xmax, x_w - gap), prelim_xmax),
    
    # Ensure valid box dimensions
    xmin = ifelse(xmin >= xmax, xmax - 1e-4, xmin),
    
    # Whisker cap boundaries
    xlc = x_w - cap_width,
    xrc = x_w + cap_width
  )

# Calculate y-axis upper limit with headroom
y_lim_top <- max(stats_all$p99, na.rm = TRUE) * 1.05

# ==============================================================================
# CREATE HALF-BOXPLOT VISUALIZATION (PANEL A)
# ==============================================================================

p1 <- ggplot() +
  # Draw boxes (IQR)
  geom_rect(
    data = stats_all,
    aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3,
        fill = analyte, colour = analyte),
    alpha = 0.6
  ) +
  
  # Draw median lines
  geom_segment(
    data = stats_all,
    aes(x = xmin, xend = xmax, y = med, yend = med, color = analyte),
    linewidth = 0.6
  ) +
  
  # Draw upper whiskers (median to 99th percentile)
  geom_segment(
    data = stats_all,
    aes(x = x_w, xend = x_w, y = med, yend = p99, color = analyte),
    linewidth = 0.5
  ) +
  
  # Draw upper whisker caps
  geom_segment(
    data = stats_all,
    aes(x = xlc, xend = xrc, y = p99, yend = p99, color = analyte),
    linewidth = 0.5
  ) +
  
  # Draw lower whiskers (1st percentile to median)
  geom_segment(
    data = stats_all,
    aes(x = x_w, xend = x_w, y = p1, yend = med, color = analyte),
    linewidth = 0.5
  ) +
  
  # Draw lower whisker caps
  geom_segment(
    data = stats_all,
    aes(x = xlc, xend = xrc, y = p1, yend = p1, color = analyte),
    linewidth = 0.5
  ) +
  
  # X-axis: analyte names
  scale_x_continuous(
    breaks = sort(unique(stats_all$x)),
    labels = levels(df_nm$analyte),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  # Y-axis: square root scale with dual axes
  scale_y_sqrt(
    name   = "Absolute concentration [nmol/L]",
    limits = c(0, y_lim_top),
    breaks = c(0, 100, 500, 1000, 2000, 4000),
    sec.axis = sec_axis(
      ~ . / s_factor,
      name   = "Relative concentration [% of apoB]",
      breaks = c(0.01, 1, 10, 25, 50, 75, 100),
      labels = c("0%", "1%", "10%", "25%", "50%", "75%", "100%")
    )
  ) +
  
  # Apply color schemes
  scale_fill_manual(values = my_cols) +
  scale_color_manual(values = my_cols) +
  
  # Theme settings
  theme_linedraw(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x       = element_blank(),
    plot.margin        = margin(t = 30, r = 10, b = 10, l = 15)
  )

# ==============================================================================
# DATA PREPARATION: CORRELATION HEATMAP
# ==============================================================================

# Helper function to apply nice names and factor ordering
factorize_nice <- function(x) {
  factor(nice_names[x], levels = nice_names[vars_abs])
}

# Create grid of all variable pairs
pair_grid <- expand_grid(row_var = vars_abs, col_var = vars_abs) %>%
  mutate(
    i = match(row_var, vars_abs),
    j = match(col_var, vars_abs),
    
    # Determine which columns to use for correlation
    # Upper triangle: use relative concentrations (% apoB)
    # Lower triangle: use absolute concentrations
    row_col_to_use = case_when(
      j > i ~ if_else(row_var == "conven_apob", "conven_apob", 
                      paste0(row_var, "_pct")),
      j < i ~ row_var,
      TRUE  ~ NA_character_
    ),
    col_col_to_use = case_when(
      j > i ~ if_else(col_var == "conven_apob", "conven_apob", 
                      paste0(col_var, "_pct")),
      j < i ~ col_var,
      TRUE  ~ NA_character_
    )
  )

# Function to compute Spearman correlation
spearman_corr <- function(x, y, data) {
  if (is.na(x) || is.na(y)) return(NA_real_)
  if (!all(c(x, y) %in% names(data))) return(NA_real_)
  
  v1 <- data[[x]]
  v2 <- data[[y]]
  
  if (!is.numeric(v1)) v1 <- suppressWarnings(as.numeric(v1))
  if (!is.numeric(v2)) v2 <- suppressWarnings(as.numeric(v2))
  
  suppressWarnings(cor(v1, v2, method = "spearman", use = "complete.obs"))
}

# Compute correlations for all pairs
pair_grid <- pair_grid %>%
  mutate(
    rho = pmap_dbl(
      list(row_col_to_use, col_col_to_use),
      ~ spearman_corr(..1, ..2, dataset)
    )
  ) %>%
  mutate(
    row_lab = factorize_nice(row_var),
    col_lab = factorize_nice(col_var),
    rho_lab = ifelse(is.na(rho), "", sprintf("%.2f", rho))
  )

# ==============================================================================
# CREATE CORRELATION HEATMAP (PANEL B)
# ==============================================================================

p2 <- ggplot(pair_grid, aes(x = col_lab, y = row_lab, fill = rho)) +
  # Draw heatmap tiles
  geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
  
  # Add correlation values as text
  geom_text(aes(label = rho_lab), size = 3.5) +
  
  # Color scale: blue (negative) to white (zero) to red (positive)
  scale_fill_gradient2(
    name     = "Spearman \u03C1",
    low      = "#2c7bb6",
    mid      = "white",
    high     = "#c14f4b",
    midpoint = 0,
    limits   = c(-1, 1),
    oob      = scales::squish
  ) +
  
  # Square aspect ratio
  coord_equal() +
  
  # Labels
  labs(x = NULL, y = NULL) +
  
  # Theme settings
  theme_linedraw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.margin    = margin(t = 2, r = 0, b = 2, l = 0),
    plot.margin      = margin(t = 5, r = 5, b = 0, l = 5)
  )

# ==============================================================================
# COMBINE PANELS INTO FINAL FIGURE
# ==============================================================================

# Combine panels A and B with labels
combined <- cowplot::plot_grid(
  p1, p2,
  ncol        = 2,
  rel_widths  = c(1, 0.8),
  align       = "v",
  axis        = "tb",
  labels      = c("A", "B"),
  label_size  = 14,
  label_fontface = "plain",
  label_x     = 0.02,
  label_y     = 0.98,
  hjust       = 0,
  vjust       = 1
)

# Final layout adjustment
combined <- cowplot::plot_grid(combined, ncol = 1, rel_heights = c(1))

# Display figure
print(combined)

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

# Save high-resolution figure
ggsave(
  "sfigure-lipoprotein.png",
  plot   = combined,
  width  = 8,
  height = 3.65,
  dpi    = 1200,
  bg     = "white"
)
