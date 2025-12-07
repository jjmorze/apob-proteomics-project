# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-27
# Last modified: 2025-12-07
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: ggplot2, cowplot, patchwork, dplyr
#
# Goal:
#   Generate scatter plots comparing beta coefficients across different
#   measurement methods (conventional vs NMR/Olink) with Spearman correlation
#   statistics. Creates a three-panel figure showing:
#   1) Conventional TRL vs NMR TRL
#   2) Conventional LDL vs NMR LDL
#   3) Conventional Lp(a) vs Olink Lp(a)
#
# Required input data:
#   - model_1.rda
#   - model_2.rda
#   - model_3.rda
#
# Main outputs:
#   - sfigure-scatter-nmr-olink.png 
#
# ==============================================================================


# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

# Install packages if needed
# install.packages(c("ggplot2", "cowplot", "patchwork", "dplyr"))

library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)


# ==============================================================================
#  FUNCTION DEFINITIONS
# ==============================================================================

#' Create Scatter Plot with Spearman Correlation
#'
#' Generates a scatter plot with linear trend line and annotated Spearman
#' correlation coefficient and p-value.
#'
#' @param df Data frame containing the variables to plot
#' @param x Character string naming the x-axis variable
#' @param y Character string naming the y-axis variable
#' @param point_alpha Numeric value for point transparency (default: 0.6)
#' @param xlab Custom label for x-axis (default: variable name)
#' @param ylab Custom label for y-axis (default: variable name)
#'
#' @return A ggplot object
#'
scatter_with_spearman <- function(df, x, y,
                                  point_alpha = 0.6,
                                  xlab = NULL,
                                  ylab = NULL) {
  # Validate inputs
  stopifnot(x %in% names(df), y %in% names(df))
  
  # Extract variables
  xv <- df[[x]]
  yv <- df[[y]]
  
  # Calculate Spearman correlation
  ct   <- suppressWarnings(cor.test(xv, yv, method = "spearman", exact = FALSE))
  rho  <- unname(ct$estimate)
  pval <- ct$p.value
  
  # Format p-value for display
  if (pval < 2.2e-16) {
    pval_text <- "p < 2.2e-16"
  } else {
    pval_text <- paste0("p = ", format(pval, scientific = TRUE, digits = 2))
  }
  
  # Create two-line correlation label
  label_text <- sprintf("Spearman \u03C1 = %.3f\n%s", rho, pval_text)
  
  # Generate plot
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = point_alpha) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "#E40428") +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.07,
      vjust = 1.25,
      size = 4
    ) +
    labs(
      x = ifelse(is.null(xlab), x, xlab),
      y = ifelse(is.null(ylab), y, ylab)
    ) +
    theme_linedraw(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA)
    )
}


# ==============================================================================
#  DATA PREPARATION
# ==============================================================================

# Merge model results
# Note: Assumes model_1, model_2, and model_3 are already loaded in environment
data_plot <- left_join(model_1, model_2, by = "gene")
data_plot <- left_join(data_plot, model_3, by = "gene")


# ==============================================================================
#  PLOT GENERATION
# ==============================================================================

# Panel A: Conventional TRL vs NMR TRL
p1 <- scatter_with_spearman(
  data_plot,
  x = "estimate.ukbb.trl",
  y = "estimate.ukbb.nmr.trl",
  xlab = "Beta - conventional TRL",
  ylab = "Beta - NMR TRL"
)

# Panel B: Conventional LDL vs NMR LDL
p2 <- scatter_with_spearman(
  data_plot,
  x = "estimate.ukbb.ldl",
  y = "estimate.ukbb.nmr.ldl",
  xlab = "Beta - conventional LDL",
  ylab = "Beta - NMR LDL"
)

# Panel C: Conventional Lp(a) vs Olink Lp(a)
p3 <- scatter_with_spearman(
  data_plot,
  x = "estimate.ukbb.lpa",
  y = "estimate.ukbb.olinklpa.lpa",
  xlab = "Beta - conventional Lp(a)",
  ylab = "Beta - Olink Lp(a)"
)


# ==============================================================================
#  COMBINE AND SAVE
# ==============================================================================

# Combine panels into three-column layout
combined <- p1 + p2 + p3 +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

# Add subplot labels (A, B, C)
combined <- combined + plot_annotation(tag_levels = "A")

# Save figure with high resolution and transparent background
ggsave(
  filename = "sfigure-scatter-nmr-olink.png",
  plot = combined,
  width = 12,
  height = 3.5,
  bg = "transparent",
  dpi = 1200
)