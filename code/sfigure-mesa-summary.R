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
#   Required packages: dplyr, ggplot2, paletteer, patchwork
#
# Goal:
#   Create validation visualizations comparing UK Biobank MPS (Multi-Protein 
#   Scores) associations with MESA (Multi-Ethnic Study of Atherosclerosis) 
#   validation results. Generate 100% stacked bar plots showing agreement,
#   disagreement, and non-associated proteins across lipoprotein fractions.
#
# Required input data:
#   - mesa_validation_with_ukbb.rda  
#
# Main outputs:
#   - sfigure-mesa-summary.png (validation summary plot)
#
# ==============================================================================


# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(paletteer)
library(patchwork)


# ==============================================================================
# PREPARE COMPOSITION DATA (100% STACKED)
# ==============================================================================

# Calculate agreement counts and percentages for each lipoprotein fraction
mesa_comp_counts <- mesa_validation_with_ukbb %>%
  mutate(
    # Ensure consistent fraction naming and ordering
    fraction = factor(fraction, levels = c("Lp(a)", "TRL", "LDL")),
    
    # Classify agreement status (function defined elsewhere: agree_mesa)
    agree = agree_mesa(direction.ukbb, estimate.mesa, p.value.mesa),
    agree = factor(agree, levels = c("-", "0", "+"))
  ) %>%
  
  # Count proteins by fraction and agreement status
  group_by(fraction, agree) %>%
  summarise(n_proteins = n_distinct(protein_name), .groups = "drop") %>%
  
  # Calculate percentages within each fraction
  group_by(fraction) %>%
  mutate(
    total = sum(n_proteins, na.rm = TRUE),
    pct   = ifelse(total > 0, 100 * n_proteins / total, 0)
  ) %>%
  ungroup()


# ==============================================================================
# DEFINE COLORS AND LABELS
# ==============================================================================

# Extract 3 colors from miami2 palette for agreement categories
miami2_cols3 <- as.character(paletteer_d("NineteenEightyR::miami2", n = 5))[c(2, 3, 4)]
names(miami2_cols3) <- c("+", "0", "-")

# Legend labels for agreement categories
legend_labs_mesa <- c(
  "-" = "Disagreement",
  "0" = "Not associated",
  "+" = "Agreement"
)

# Consistent y-axis ordering
frac_levels <- levels(mesa_comp_counts$fraction)


# ==============================================================================
# CREATE MAIN 100% STACKED BAR PLOT
# ==============================================================================

p_comp_mesa <- ggplot(mesa_comp_counts, aes(x = pct, y = fraction, fill = agree)) +
  
  # Stacked horizontal bars
  geom_col(width = 0.6) +
  
  # Add labels showing percentage and count (hide if < 2%)
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
  
  # Apply color scheme
  scale_fill_manual(
    values = miami2_cols3,
    breaks = names(legend_labs_mesa),
    labels = legend_labs_mesa,
    name   = "MESA validation"
  ) +
  
  # X-axis: 0-100% with proper formatting
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  # Y-axis: lipoprotein fractions
  scale_y_discrete(limits = frac_levels) +
  
  # Labels
  labs(
    x = "Share of protein associations",
    y = "Lipoprotein fraction"
  ) +
  
  # Theme customization
  theme_linedraw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )


# ==============================================================================
# CREATE PROTEIN COUNT COLUMN
# ==============================================================================

# Prepare data frame with total protein counts per fraction
mesa_totals_df <- mesa_comp_counts %>% 
  distinct(fraction, total)

# Create slim column showing total N
p_totals_mesa <- ggplot(mesa_totals_df, aes(y = fraction)) +
  
  # Display total count
  geom_text(aes(x = 1, label = total), size = 4.2, fontface = "bold") +
  
  # Match y-axis with main plot
  scale_y_discrete(limits = frac_levels) +
  
  # Minimal x-axis space
  coord_cartesian(xlim = c(0.5, 1.5), clip = "off") +
  
  # Column header
  labs(title = "N of\nproteins") +
  
  # Minimal theme (remove borders, grids, axes)
  theme_linedraw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = -4)
  )


# ==============================================================================
# COMBINE PLOTS
# ==============================================================================

# Combine main plot with totals column using patchwork
p_final_mesa <- (p_comp_mesa + p_totals_mesa + 
                   plot_layout(widths = c(4, 0.9), guides = "collect")) &
  theme(legend.position = "right")

# Display combined plot
p_final_mesa


# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

ggsave(
  filename = "sfigure-mesa-summary.png",
  plot = p_final_mesa,
  width = 10.5,
  height = 4.2,
  dpi = 1200
)