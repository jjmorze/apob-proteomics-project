# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-12-07
# Last modified: 2025-12-07
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.0
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, ggplot2, paletteer, patchwork
#
# Goal:
#   Create composition plots showing consistency of Mendelian Randomization (MR)
#   results with observational analyses across lipoprotein fractions.
#   Generates a 100% stacked bar chart with:
#   1) Agreement categories (++, +, 0, -, --) for protein associations
#   2) Share of proteins in each category per lipoprotein fraction
#   3) Total protein counts displayed alongside
#
# Required input data:
#   - combined_long.rda 
#
# Main outputs:
#   - sfigure-mr-summary.png
#
# ==============================================================================


# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(paletteer)
library(patchwork)  # For layout side-by-side and legend collection


# ==============================================================================
# PREPARE 100% COMPOSITION DATA
# ==============================================================================

comp_counts <- combined_long %>%
  mutate(
    agree    = factor(agree, levels = c("++", "+", "0", "-", "--")),
    fraction = factor(fraction, levels = c("Lp(a)", "TRL", "LDL"))
  ) %>%
  group_by(fraction, agree) %>%
  summarise(n_proteins = n_distinct(protein_name), .groups = "drop") %>%
  group_by(fraction) %>%
  mutate(
    total = sum(n_proteins),
    pct   = ifelse(total > 0, 100 * n_proteins / total, 0)
  ) %>%
  ungroup()


# ==============================================================================
# DEFINE COLOR PALETTE AND LEGEND LABELS
# ==============================================================================

# Miami2 color palette from NineteenEightyR
miami2_cols <- as.character(paletteer_d("NineteenEightyR::miami2", n = 5))
names(miami2_cols) <- c("++", "+", "0", "-", "--")

# Manual legend labels for agreement categories
legend_labs <- c(
  "++" = "Agreement",
  "+"  = "Nominal agreement",
  "0"  = "Not associated",
  "-"  = "Nominal disagreement",
  "--" = "Disagreement"
)

# Consistent y-axis order for both plots
frac_levels <- levels(comp_counts$fraction)


# ==============================================================================
# CREATE MAIN 100% STACKED BAR PLOT
# ==============================================================================

p_comp <- ggplot(comp_counts, aes(x = pct, y = fraction, fill = agree)) +
  geom_col(width = 0.6) +
  
  # Add text labels: line 1 = percent, line 2 = count (hide if < 2%)
  geom_text(
    aes(label = ifelse(
      pct >= 2,
      paste0(round(pct, 1), "%\n", n_proteins),
      ""
    )),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    lineheight = 0.9,
    size = 3.6
  ) +
  
  # Color scale with custom labels
  scale_fill_manual(
    values = miami2_cols,
    breaks = names(legend_labs),
    labels = legend_labs,
    name   = "MR outcome"
  ) +
  
  # X-axis: percentage scale
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  # Y-axis: lipoprotein fractions
  scale_y_discrete(limits = frac_levels) +
  
  # Labels and theme
  labs(
    x = "Share of protein associations",
    y = "Lipoprotein fraction"
  ) +
  theme_linedraw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )


# ==============================================================================
# CREATE TOTALS COLUMN PANEL
# ==============================================================================

# Extract total protein counts per fraction
totals_df <- comp_counts %>% 
  distinct(fraction, total)

# Create slim totals panel (between plot and legend)
p_totals <- ggplot(totals_df, aes(y = fraction)) +
  geom_text(aes(x = 1, label = total), size = 4.2, fontface = "bold") +
  scale_y_discrete(limits = frac_levels) +
  coord_cartesian(xlim = c(0.5, 1.5), clip = "off") +
  labs(title = "N of\nproteins") +
  theme_linedraw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = -2)
  )


# ==============================================================================
# COMBINE PLOTS WITH PATCHWORK
# ==============================================================================

# Combine: totals panel sits between main plot and the collected legend
# guides='collect' moves a single legend to the right of the combined layout
p_final <- (p_comp + p_totals + 
              plot_layout(widths = c(4, 0.9), guides = "collect")) &
  theme(legend.position = "right")

# Display the final plot
p_final


# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

ggsave(
  filename = "sfigure-mr-summary.png",
  plot = p_final,
  width = 10.5,
  height = 4.2,
  dpi = 1200
)