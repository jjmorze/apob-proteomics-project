# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-26
# Last modified: 2025-11-27
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, stringr, ggplot2, grid, gridExtra, ggrepel,
#                      ggvenn, forcats, tibble, cowplot
#
# Goal:
#   Generate a multi-panel figure (Figure 3) combining:
#   - Panel A: Venn diagram showing protein signature overlap across LDL, TRL, 
#              and Lp(a) fractions
#   - Panel B: Scatter plot comparing Lp(a) vs TRL protein coefficients with 
#              Pearson correlation
#   - Panel C: Three barplots showing top enriched pathways for each fraction
#              with shared legend
#
# Required input data:
#   - combined_long_signif.rda
#   - combined.rda
#   - ea_ldl.rda
#   - ea_trl.rda
#   - ea_lpa.rda
#
# Output:
#   - figure-3.png (12.3" × 6.6", 1200 DPI)
#
# ==============================================================================


# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(ggrepel)
  library(ggvenn)
  library(forcats)
  library(tibble)
  library(cowplot)
})


# ==============================================================================
# PANEL A: VENN DIAGRAM (PROTEIN SIGNATURES)
# ==============================================================================

# Harmonize and filter protein signatures
df_sig <- combined_long_signif %>%
  filter(agree == "++") %>%
  mutate(
    # Standardize fraction labels
    fraction_std = case_when(
      str_to_upper(fraction) == "LDL" ~ "LDL",
      str_to_upper(fraction) == "TRL" ~ "TRL",
      str_to_upper(fraction) %in% c("LPA", "LP(A)") ~ "Lp(a)",
      TRUE ~ fraction
    ),
    # Map direction to +/- signs
    sign_chr = case_when(
      str_detect(direction.ukbb, regex("^pos", ignore_case = TRUE)) ~ "+",
      str_detect(direction.ukbb, regex("^neg", ignore_case = TRUE)) ~ "-",
      direction.ukbb %in% c("+", "-") ~ direction.ukbb,
      TRUE ~ direction.ukbb
    ),
    # Create gene_sign signature
    signature = paste0(gene, "_", sign_chr)
  ) %>%
  select(fraction_std, signature) %>%
  distinct()

# Build sets for each fraction
set_LDL <- df_sig %>% filter(fraction_std == "LDL") %>% pull(signature)
set_TRL <- df_sig %>% filter(fraction_std == "TRL") %>% pull(signature)
set_LPA <- df_sig %>% filter(fraction_std == "Lp(a)") %>% pull(signature)

venn_list <- list(
  LDL     = unique(set_LDL),
  TRL     = unique(set_TRL),
  "Lp(a)" = unique(set_LPA)
)

# Create Venn diagram
pA <- ggvenn(
  venn_list,
  fill_color   = c("#3F7DB8FF", "#3FB8AFFF", "#FF9E9DFF"),
  fill_alpha   = 0.65,
  stroke_color = "grey60",
  stroke_size  = 0.5,
  text_size    = 5.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))


# ==============================================================================
# PANEL B: SCATTER PLOT (Lp(a) vs TRL COEFFICIENTS)
# ==============================================================================

# Helper function to find p-value column
find_pcol <- function(df, frac_pattern) {
  p_cols <- names(df)[grepl(
    "\\bp\\b|p.value|p_value|pval|p-val|pvalue",
    names(df),
    ignore.case = TRUE
  )]
  
  if (length(p_cols) == 0) return(NA_character_)
  
  frac_cols <- p_cols[grepl(frac_pattern, p_cols, ignore.case = TRUE)]
  if (length(frac_cols) > 0) return(frac_cols[1])
  
  return(p_cols[1])
}

# Prepare data for scatter plot
df_plot <- combined %>%
  filter(agree.lpa == "++") %>%
  select(any_of(c("gene", "estimate.ukbb.lpa", "estimate.ukbb.trl")), everything()) %>%
  rename(x = estimate.ukbb.lpa, y = estimate.ukbb.trl) %>%
  filter(!is.na(x) & !is.na(y))

# Find and add p-value ranking
pcol_lpa <- find_pcol(df_plot, "lpa|lp\\(?a\\)?")
pcol_trl <- find_pcol(df_plot, "trl")

if (!is.na(pcol_lpa) && pcol_lpa %in% names(df_plot) &&
    !is.na(pcol_trl) && pcol_trl %in% names(df_plot)) {
  df_plot <- df_plot %>%
    mutate(p_rank = pmin(
      !!rlang::sym(pcol_lpa),
      !!rlang::sym(pcol_trl),
      na.rm = TRUE
    ))
} else if (!is.na(pcol_lpa) && pcol_lpa %in% names(df_plot)) {
  df_plot <- df_plot %>% mutate(p_rank = !!rlang::sym(pcol_lpa))
} else if (!is.na(pcol_trl) && pcol_trl %in% names(df_plot)) {
  df_plot <- df_plot %>% mutate(p_rank = !!rlang::sym(pcol_trl))
} else {
  stop("No p-value column found. Include a column like p.value or p_value.")
}

# Label top 10 and bottom 10 proteins by p-value
df_plot <- df_plot %>%
  arrange(p_rank) %>%
  mutate(label_flag = NA_character_) %>%
  mutate(label_flag = if_else(row_number() <= 10, gene, label_flag)) %>%
  mutate(label_flag = if_else(row_number() > n() - 10, gene, label_flag))

# Calculate Pearson correlation
cor_test <- cor.test(df_plot$x, df_plot$y, method = "pearson")
r <- round(unname(cor_test$estimate), 2)
pval <- ifelse(
  cor_test$p.value < 2.2e-16,
  "P < 2.2 × 10⁻¹⁶",
  paste0("P = ", signif(cor_test$p.value, 2))
)
label_text <- paste0("R = ", r, ", ", pval)

# Create scatter plot
pB <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(color = "grey40", size = 2.0, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.9) +
  geom_text_repel(
    aes(label = label_flag),
    size = 3.0,
    na.rm = TRUE,
    max.overlaps = Inf,
    segment.size = 0.2
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label_text,
    hjust = -0.05, vjust = 1.3,
    size = 4
  ) +
  labs(
    x = expression("Protein coefficients for Lp(a)"),
    y = expression("Protein coefficients for TRL")
  ) +
  theme_linedraw(base_size = 12) +
  theme(
    plot.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 140, 5.5, 75)
  )


# ==============================================================================
# PANEL C: PATHWAY ENRICHMENT BARPLOTS
# ==============================================================================

# Manual pathway mapping to labels and themes
mapping <- tribble(
  ~Term, ~Label, ~Theme,
  # TRL pathways
  "Chemokine-Mediated Signaling Pathway (GO:0070098)",
  "Chemokine-mediated signaling", "Inflammation",
  "Cellular Response To Chemokine (GO:1990869)",
  "Cellular response to chemokine", "Inflammation",
  "Lymphocyte Chemotaxis (GO:0048247)",
  "Lymphocyte chemotaxis", "Inflammation",
  "Positive Regulation Of ERK1 And ERK2 Cascade (GO:0070374)",
  "↑regulation of ERK cascade", "Signalling",
  "Neutrophil Chemotaxis (GO:0030593)",
  "Neutrophil chemotaxis", "Inflammation",
  "Regulation Of ERK1 And ERK2 Cascade (GO:0070372)",
  "Regulation of ERK cascade", "Signalling",
  "Granulocyte Chemotaxis (GO:0071621)",
  "Granulocyte chemotaxis", "Inflammation",
  "Monocyte Chemotaxis (GO:0002548)",
  "Monocyte chemotaxis", "Inflammation",
  "Neutrophil Migration (GO:1990266)",
  "Neutrophil migration", "Inflammation",
  "Positive Regulation Of MAPK Cascade (GO:0043410)",
  "↑regulation of MAPK cascade", "Signalling",
  # LDL pathways
  "Cholesterol Homeostasis (GO:0042632)",
  "Cholesterol homeostasis", "Lipid/Lipoprotein",
  "Sterol Homeostasis (GO:0055092)",
  "Sterol homeostasis", "Lipid/Lipoprotein",
  "Regulation Of Inflammatory Response (GO:0050727)",
  "Regulation of inflammatory response", "Inflammation",
  "Platelet Activating Factor Metabolic Process (GO:0046469)",
  "PAF metabolic process", "Coagulation/Hemostasis",
  "High-Density Lipoprotein Particle Assembly (GO:0034380)",
  "HDL particle assembly", "Lipid/Lipoprotein",
  "Phosphatidylcholine Catabolic Process (GO:0034638)",
  "PC catabolism", "Lipid/Lipoprotein",
  "Negative Regulation Of Hemostasis (GO:1900047)",
  "↓regulation of hemostasis", "Coagulation/Hemostasis",
  "Low-Density Lipoprotein Particle Remodeling (GO:0034374)",
  "LDL particle remodeling", "Lipid/Lipoprotein",
  "Positive Regulation Of T-helper 1 Type Immune Response (GO:0002827)",
  "↑regulation of Th1 response", "Inflammation",
  "Negative Regulation Of Coagulation (GO:0050819)",
  "↓regulation of coagulation", "Coagulation/Hemostasis",
  # Lp(a) pathways
  "Neuropeptide Signaling Pathway (GO:0007218)",
  "Neuropeptide signaling", "Signalling",
  "Inflammatory Response (GO:0006954)",
  "Inflammatory response", "Inflammation",
  "Response To Type II Interferon (GO:0034341)",
  "Response to type II interferon", "Inflammation",
  "Response To Interleukin-1 (GO:0070555)",
  "Response to IL-1", "Inflammation"
)

# Define theme colors
theme_colors <- c(
  "Inflammation"           = "#CFAF8BCC",
  "Lipid/Lipoprotein"      = "#DAD8A7FF",
  "Signalling"             = "#B772B0CC",
  "Coagulation/Hemostasis" = "#AEC8A8CC"
)
theme_levels <- names(theme_colors)

# Helper function to format p-values
fmt_p <- function(p) {
  ifelse(p < 1e-4, formatC(p, format = "e", digits = 2), signif(p, 2))
}

# Function to prepare data for each panel
prep_panel <- function(df_full, panel_title) {
  stopifnot(all(c("Term", "P.value") %in% names(df_full)))
  
  df_full %>%
    mutate(P.value = as.numeric(P.value)) %>%
    arrange(P.value) %>%
    slice_head(n = 10) %>%
    left_join(mapping, by = "Term") %>%
    mutate(
      Theme = factor(
        if_else(is.na(Theme), "Other", Theme),
        levels = c(theme_levels, "Other")
      ),
      Label = if_else(
        is.na(Label),
        str_replace(Term, "\\s*\\(GO:\\d+\\)\\s*$", ""),
        Label
      ),
      neglog10P = -log10(P.value),
      Label = fct_reorder(Label, neglog10P),
      p_lab = fmt_p(P.value),
      Panel = panel_title
    )
}

# Function to create barplot
make_plot <- function(df_prepped, title_txt, show_legend = TRUE) {
  max_y <- if (nrow(df_prepped) > 0) {
    max(df_prepped$neglog10P, na.rm = TRUE)
  } else {
    1
  }
  
  ggplot(df_prepped, aes(x = Label, y = neglog10P, fill = Theme)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = p_lab, y = neglog10P), hjust = 1.05, size = 3) +
    coord_flip(clip = "off") +
    expand_limits(y = max_y * 1.10) +
    scale_fill_manual(
      values = theme_colors,
      breaks = theme_levels,
      drop = TRUE,
      name = "Meta-pathway"
    ) +
    labs(
      x = NULL,
      y = expression(-log[10](P)),
      title = title_txt
    ) +
    theme_linedraw(base_size = 11) +
    theme(
      panel.border = element_rect(linewidth = 0.3),
      plot.title = element_text(hjust = 0.5, face = "plain"),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(5.5, 18, 5.5, 5.5)
    )
}

# Prepare data for each fraction
df_ldl <- prep_panel(ea_ldl$full, "LDL")
df_trl <- prep_panel(ea_trl$full, "TRL")
df_lpa <- prep_panel(ea_lpa$full, "Lp(a)")

# Create shared legend
df_all <- bind_rows(df_ldl, df_trl, df_lpa)
legend_breaks <- intersect(theme_levels, unique(as.character(df_all$Theme)))

legend_plot <- ggplot(df_all, aes(Label, neglog10P, fill = Theme)) +
  geom_col() +
  scale_fill_manual(
    values = theme_colors,
    breaks = legend_breaks,
    drop = TRUE,
    name = "Meta-pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

shared_legend <- cowplot::get_legend(legend_plot)

# Create individual plots without legends
p_ldl <- make_plot(df_ldl, "LDL", show_legend = FALSE)
p_trl <- make_plot(df_trl, "TRL", show_legend = FALSE)
p_lpa <- make_plot(df_lpa, "Lp(a)", show_legend = FALSE)

# Combine plots with shared legend
row_plots <- cowplot::plot_grid(
  p_ldl, p_trl, p_lpa,
  ncol = 3,
  rel_widths = c(1.1, 1, 1),
  align = "h"
)

pC <- cowplot::plot_grid(
  row_plots, shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.15)
)


# ==============================================================================
# COMBINE ALL PANELS
# ==============================================================================

# Add panel labels (A, B, C)
pA_lab <- ggdraw(pA) +
  draw_plot_label(
    "A",
    x = 0.01, y = 0.99,
    hjust = 0, vjust = 1,
    size = 14,
    fontface = "plain"
  )

pB_lab <- ggdraw(pB) +
  draw_plot_label(
    "C",
    x = 0.01, y = 0.99,
    hjust = 0, vjust = 1,
    size = 14,
    fontface = "plain"
  )

pC_lab <- ggdraw(pC) +
  draw_plot_label(
    "B",
    x = 0.01, y = 0.99,
    hjust = 0, vjust = 1,
    size = 14,
    fontface = "plain"
  )

# Define layout: A and B in first row, C spanning second row
layout_mat <- rbind(
  c(1, 2),
  c(3, 3)
)

# Define relative dimensions
first_row_widths <- c(2, 2.54)  # Width ratio for A:B
row_heights <- c(1, 1)          # Height ratio for rows

# Combine all panels
final_plot <- grid.arrange(
  grobs = list(pA_lab, pB_lab, pC_lab),
  layout_matrix = layout_mat,
  widths = first_row_widths,
  heights = row_heights
)


# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

ggsave(
  "figure-3.png",
  final_plot,
  width = 12.3,
  height = 6.6,
  dpi = 1200,
  bg = "white"
)