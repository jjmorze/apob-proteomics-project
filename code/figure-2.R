# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2024-09-26
# Last modified: 2025-11-27
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, tidyr, ggplot2, ggrepel, cowplot, RColorBrewer,
#                      stringr, tibble, scales
#
# Goal:
#   Generate Figure 2: Multi-trait volcano and forest plots for lipoprotein 
#   particle-protein associations.
#   - Panel A: Pathway annotation bar (biological class)
#   - Panel B: Volcano plots (LDL, TRL, Lp(a))
#   - Panel C: Forest plots comparing UKBB and MESA cohorts
#   Includes MR-confirmed associations and pathway-based coloring.
#
# Required input data:
#   - combined.rda
#   - mesa_validation_with_ukbb.rda
# Outputs:
#   - figure-2.png: Final composite figure (10 x 12.5 inches, 600 DPI)
#
# ==============================================================================


# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(RColorBrewer)
  library(stringr)
  library(tibble)
  library(scales)
}) 

# ==============================================================================
#  LOAD DATA
# ==============================================================================

load(combined.rda)
load(mesa_validation_with_ukbb.rda)

# ==============================================================================
#  CONSTANTS AND CONFIGURATION
# ==============================================================================

# Statistical thresholds
SMALL_P_FLOOR <- 1e-300
BONF_THRESH   <- 0.05 / 1417  # Bonferroni correction for 1417 tests

# Cohort color scheme
COHORT_COLORS <- c(
  UKBB = "steelblue",
  MESA = "#8B2D25"
)

# Pathway color palette (biological classes)
PATHWAY_COLORS <- c(
  "Cell Cycle and Apoptosis Regulators"          = "#5C93C2FF",
  "Gene Expression and Regulatory Proteins"      = "#6E9FC5CC",
  "Hormones and Growth Factors"                  = "#7FA7C7FF",
  "Immune and Defense Proteins"                  = "#5FB2BBCC",
  "Metabolic Enzymes"                            = "#3FB8AFFF",
  "Others"                                       = "#5FBEAECC",
  "Plasma and Extracellular Proteins"            = "#7FC7AFFF",
  "Protein Synthesis, Folding, and Degradation"  = "#AEC8A8CC",
  "Receptors and Signal Transducers"             = "#DAD8A7FF",
  "Structural and Cytoskeletal Proteins"         = "#F5BBA2CC",
  "Transporters and Ion Channels"                = "#FF9E9DFF"
)


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

# Annotate statistical significance
# 
# Categorizes associations into three groups:
#   - ns: not significant (p_obs >= Bonferroni threshold)
#   - obs_sig: observational significance only
#   - mr_confirmed: both observational and MR-confirmed
#
# @param df Data frame with p_obs and p_mr columns
# @param bonf_thresh Bonferroni significance threshold
# @return Data frame with sig_group column added
annotate_significance <- function(df, bonf_thresh = BONF_THRESH) {
  df %>%
    mutate(
      # Floor zero p-values to avoid -Inf
      p_obs = if_else(p_obs == 0, SMALL_P_FLOOR, p_obs),
      p_mr  = if_else(is.na(p_mr), 1, p_mr),
      p_mr  = if_else(p_mr == 0, SMALL_P_FLOOR, p_mr),
      
      # Classify significance
      sig_group = case_when(
        p_obs >= bonf_thresh                      ~ "ns",
        p_obs < bonf_thresh & p_mr >= bonf_thresh ~ "obs_sig",
        p_obs < bonf_thresh & p_mr < bonf_thresh  ~ "mr_confirmed"
      )
    )
}


# Create pathway color palette
# 
# Generates a color palette for pathways, using manual colors where available
# and automatically generating colors for unmapped pathways.
#
# @param pathways Vector of unique pathway names
# @param manual Named vector of manual color assignments (optional)
# @return Named vector of colors for all pathways
make_pathway_palette <- function(pathways, manual = NULL) {
  pathways <- sort(unique(pathways))
  
  if (!is.null(manual)) {
    # Remove empty or NA names from manual palette
    manual <- manual[!is.na(names(manual)) & nchar(names(manual)) > 0]
    
    # Separate known and unknown pathways
    known   <- intersect(pathways, names(manual))
    unknown <- setdiff(pathways, names(manual))
    
    # Auto-generate colors for unknown pathways
    auto <- if (length(unknown) <= 12) {
      RColorBrewer::brewer.pal(
        n = max(3, length(unknown)), 
        "Paired"
      )[seq_along(unknown)]
    } else {
      scales::hue_pal()(length(unknown))
    }
    
    # Combine manual and auto colors
    pal <- c(
      setNames(manual[known], known),
      setNames(auto, unknown)
    )
    
    return(pal[pathways])
    
  } else {
    # No manual palette: auto-generate all colors
    if (length(pathways) <= 12) {
      base <- RColorBrewer::brewer.pal(
        n = max(3, length(pathways)), 
        "Paired"
      )[1:length(pathways)]
    } else {
      base <- scales::hue_pal()(length(pathways))
    }
    
    return(setNames(base, pathways))
  }
}


# Create volcano plot
# 
# Generates a bidirectional volcano plot with:
#   - Y-axis split at zero (positive/negative associations)
#   - Compressed scale for highly significant p-values
#   - Pathway-based coloring for MR-confirmed associations
#   - Automatic labeling of top associations per pathway
#
# @param data Data frame with columns: pathway, gene, beta_obs, p_obs, p_mr
# @param bonf_thresh Bonferroni significance threshold
# @param small_p_divisor Scale compression factor for p < bonf_thresh
# @param cap_at_p Optional p-value cap for y-axis limits
# @param fixed_max_x Fixed maximum x-axis value for alignment
# @param pathway_palette Manual pathway color palette (optional)
# @return ggplot object
plot_volcano <- function(data,
                         bonf_thresh = BONF_THRESH,
                         small_p_divisor = 10,
                         cap_at_p = NULL,
                         fixed_max_x = NULL,
                         pathway_palette = NULL) {
  
  # Prepare data
  df <- data %>%
    arrange(pathway, gene) %>%
    mutate(x_id = row_number()) %>%
    annotate_significance(bonf_thresh)
  
  # Custom p-value transformation (compressed scale for p < bonf_thresh)
  half_log_bonf <- -log10(bonf_thresh)
  custom_p_transform <- function(p) {
    if_else(
      p >= bonf_thresh,
      -log10(p),
      -log10(p) / small_p_divisor + half_log_bonf
    )
  }
  
  # Calculate plot coordinates
  df <- df %>%
    mutate(
      p_trans = custom_p_transform(p_obs),
      y_plot  = if_else(beta_obs >= 0, 2 + p_trans, -2 - p_trans),
      
      # Visual attributes
      alpha = case_when(
        sig_group == "ns"           ~ 0.1,
        sig_group == "obs_sig"      ~ 0.5,
        sig_group == "mr_confirmed" ~ 1
      ),
      shape_group = sig_group,
      border_group = case_when(
        sig_group == "ns"           ~ "ns",
        sig_group == "obs_sig"      ~ pathway,
        sig_group == "mr_confirmed" ~ "special"
      ),
      fill_group = if_else(sig_group == "mr_confirmed", pathway, NA_character_)
    )
  
  # Color mapping
  pathways <- sort(unique(df$pathway))
  pathway_colors <- make_pathway_palette(pathways, manual = pathway_palette)
  border_colors  <- c(ns = "grey70", special = "grey20", pathway_colors)
  fill_colors    <- pathway_colors
  shape_map      <- c(ns = 21, obs_sig = 21, mr_confirmed = 23)
  
  # Identify top associations for labeling (one positive + one negative per pathway)
  top_pos <- df %>%
    filter(sig_group == "mr_confirmed", beta_obs > 0) %>%
    group_by(pathway) %>%
    slice_max(beta_obs, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  top_neg <- df %>%
    filter(sig_group == "mr_confirmed", beta_obs < 0) %>%
    group_by(pathway) %>%
    slice_min(beta_obs, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  top_labeled <- bind_rows(top_pos, top_neg)
  
  # Reference p-value lines
  ref_pvals <- if (is.null(cap_at_p)) {
    c(1e-50, 1e-100, 1e-200, 1e-300)
  } else {
    c(cap_at_p, 1e-25, 1e-10)
  }
  ref_labels <- parse(text = paste0("10^", log10(ref_pvals)))
  ref_yvals  <- sapply(ref_pvals, custom_p_transform)
  bonf_y     <- custom_p_transform(bonf_thresh)
  
  # X-axis limits
  max_x <- if (!is.null(fixed_max_x)) fixed_max_x else max(df$x_id)
  
  # Build plot
  p <- ggplot(df, aes(x = x_id, y = y_plot)) +
    
    # Points
    geom_point(
      aes(color = border_group, fill = fill_group, alpha = alpha, shape = shape_group),
      size = 2,
      stroke = 0.4,
      show.legend = FALSE
    )
  
  # Add labels for top associations
  if (nrow(top_labeled) > 0) {
    p <- p +
      geom_text_repel(
        data = top_labeled,
        aes(label = gene, x = x_id, y = y_plot),
        size = 2.8,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.3
      )
  }
  
  # Reference lines and styling
  p <- p +
    # Zero lines
    geom_hline(yintercept = c(-2, 2), color = "black") +
    
    # P-value reference lines
    geom_hline(
      yintercept = 2 + ref_yvals,
      linetype = "solid",
      color = "grey70",
      size = 0.1
    ) +
    geom_hline(
      yintercept = -2 - ref_yvals,
      linetype = "solid",
      color = "grey70",
      size = 0.1
    ) +
    
    # Bonferroni threshold lines
    geom_hline(
      yintercept = 2 + bonf_y,
      linetype = "dotted",
      color = "red",
      size = 0.4
    ) +
    geom_hline(
      yintercept = -2 - bonf_y,
      linetype = "dotted",
      color = "red",
      size = 0.4
    ) +
    
    # Scales
    scale_shape_manual(values = shape_map, guide = "none") +
    scale_color_manual(values = border_colors, guide = "none") +
    scale_fill_manual(values = fill_colors, na.value = NA, guide = "none") +
    scale_alpha_identity(guide = "none") +
    
    scale_x_continuous(
      limits = c(0.5, max_x + 0.5),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      breaks = c(-2, -2 - ref_yvals, 2, 2 + ref_yvals),
      labels = c("0", ref_labels, "0", ref_labels),
      expand = expansion(mult = c(0.05, 0.05)),
      limits = if (!is.null(cap_at_p)) {
        c(-(2 + custom_p_transform(cap_at_p)),
          (2 + custom_p_transform(cap_at_p)))
      } else {
        NULL
      }
    ) +
    
    # Labels and theme
    labs(
      y = expression("\u2190" ~ Negative ~ phantom("XX") ~ 
                       -log[10](P) ~ phantom("XX") ~ Positive ~ "\u2192")
    ) +
    
    theme_linedraw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_blank(),
      plot.margin = margin(2, 6, 2, 2),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )
  
  return(p)
}


# Create multi-cohort forest plot
# 
# Generates a forest plot comparing UKBB and MESA cohorts for top associations.
# Shows top N positive and top N negative MR-confirmed associations from UKBB.
#
# @param data Data frame with columns: cohort, gene, beta_obs, se, p_obs, p_mr
# @param trait_name Trait name (used for LPA filtering)
# @param bonf_thresh Bonferroni significance threshold
# @param top_n Number of top positive and negative associations to show
# @param x_label X-axis label
# @param show_legend Whether to show cohort legend
# @return ggplot object
make_forest_plot_posneg_cohorts <- function(data,
                                            trait_name,
                                            bonf_thresh = BONF_THRESH,
                                            top_n = 5,
                                            x_label = "",
                                            show_legend = FALSE) {
  
  # Get UKBB MR-confirmed associations
  ukbb_df  <- data %>% filter(cohort == "UKBB") %>% annotate_significance(bonf_thresh)
  ukbb_sig <- ukbb_df %>% filter(sig_group == "mr_confirmed")
  
  # Special handling: exclude LPA gene for Lp(a) trait
  if (toupper(trait_name) == "LPA") {
    ukbb_sig <- ukbb_sig %>% filter(gene != "LPA")
  }
  
  # Return empty plot if no significant associations
  if (nrow(ukbb_sig) == 0) {
    return(ggplot() + theme_void())
  }
  
  # Select top positive and negative associations
  top_pos <- ukbb_sig %>%
    filter(beta_obs > 0) %>%
    arrange(desc(beta_obs)) %>%
    slice_head(n = top_n)
  
  top_neg <- ukbb_sig %>%
    filter(beta_obs < 0) %>%
    arrange(beta_obs) %>%
    slice_head(n = top_n)
  
  selected_genes <- unique(c(top_pos$gene, top_neg$gene))
  
  # Prepare plotting data
  sel <- data %>%
    filter(gene %in% selected_genes) %>%
    annotate_significance(bonf_thresh) %>%
    mutate(
      lower = beta_obs - 1.96 * se,
      upper = beta_obs + 1.96 * se,
      gene  = as.character(gene)
    ) %>%
    mutate(
      gene = factor(
        gene,
        levels = ukbb_sig %>%
          filter(gene %in% selected_genes) %>%
          arrange(desc(beta_obs)) %>%
          pull(gene) %>%
          unique()
      )
    ) %>%
    arrange(desc(beta_obs), cohort)
  
  # Create alternating row backgrounds
  gene_levels <- levels(sel$gene)
  gene_df <- tibble(
    gene = gene_levels,
    row_idx = seq_along(gene_levels)
  ) %>%
    mutate(
      ymin = row_idx - 0.5,
      ymax = row_idx + 0.5,
      fill_bg = if_else(row_idx %% 2 == 0, "grey90", "white")
    )
  
  position_dodge_val <- position_dodge(width = 0.6)
  
  # Build plot
  plot <- ggplot() +
    
    # Alternating backgrounds
    geom_rect(
      data = gene_df %>% filter(fill_bg == "white"),
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "white",
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = gene_df %>% filter(fill_bg == "grey90"),
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "grey90",
      inherit.aes = FALSE
    ) +
    
    # Null effect line
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    
    # Error bars
    geom_errorbarh(
      data = sel,
      aes(xmin = lower, xmax = upper, y = gene, color = cohort),
      height = 0.2,
      position = position_dodge_val,
      stroke = 0.6,
      show.legend = show_legend
    ) +
    
    # Point estimates
    geom_point(
      data = sel,
      aes(x = beta_obs, y = gene, color = cohort),
      shape = 15,
      size = 2.5,
      stroke = 0,
      position = position_dodge_val,
      show.legend = show_legend
    ) +
    
    # Scales and theme
    scale_color_manual(
      values = COHORT_COLORS,
      name = NULL,
      guide = if (show_legend) guide_legend(nrow = 1, byrow = TRUE) else "none"
    ) +
    
    labs(x = x_label, y = NULL) +
    
    theme_linedraw(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = NA, colour = NA),
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(2, 6, 2, 2),
      axis.title.x = if (x_label == "") {
        element_text(color = "transparent")
      } else {
        element_text()
      },
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  return(plot)
}


# Create pathway annotation bar
# 
# Generates a horizontal bar showing pathway groupings with protein counts.
#
# @param volcano_data Data frame from volcano plot
# @param max_x Maximum x-axis value for alignment
# @param height Bar height
# @param offset X-axis offset
# @param pathway_palette Manual pathway color palette (optional)
# @return ggplot object
make_pathway_bar <- function(volcano_data,
                             max_x,
                             height = 1.5,
                             offset = 0,
                             pathway_palette = NULL) {
  
  # Assign x positions
  df <- volcano_data %>%
    arrange(pathway, gene) %>%
    mutate(x_id = row_number())
  
  # Calculate pathway blocks
  counts <- df %>%
    group_by(pathway) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(pathway) %>%
    mutate(
      xmin = lag(cumsum(n), default = 0) + 0.5 + offset,
      xmax = cumsum(n) + 0.5 + offset,
      xmid = (xmin + xmax) / 2,
      label_wrapped = paste0(n, " proteins")
    )
  
  # Color mapping
  fill_map <- make_pathway_palette(counts$pathway, manual = pathway_palette)
  
  # Build plot
  ggplot(counts) +
    
    # Pathway blocks
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = height, fill = pathway),
      color = "black",
      show.legend = FALSE
    ) +
    
    # Labels (rotated 90 degrees)
    geom_text(
      aes(x = xmid, y = height / 2, label = label_wrapped),
      angle = 90,
      size = 3,
      color = "black",
      lineheight = 0.95,
      hjust = 0.5,
      vjust = 0.5
    ) +
    
    # Scales and theme
    scale_fill_manual(values = fill_map, guide = "none") +
    scale_x_continuous(
      limits = c(0.5, max_x + 0.5) + offset,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    
    coord_cartesian(clip = "off") +
    
    theme_void(base_size = 13) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}


# Create side trait label
# 
# Generates a vertical text label for trait names.
#
# @param text Label text
# @return ggplot object
make_trait_side_label <- function(text) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = text,
      angle = 90,
      size = 4,
      fontface = "bold"
    ) +
    theme_void()
}


# Create column labels
# 
# Generates panel labels "B" and "C" for volcano and forest plots.
#
# @return List with volcano and forest label plots
make_column_labels <- function() {
  volcano_label <- ggdraw() +
    draw_label("B", fontface = "plain", size = 14, x = -0.07, hjust = 0, vjust = 0.75)
  
  forest_label <- ggdraw() +
    draw_label("C", fontface = "plain", size = 14, x = 0, hjust = 0, vjust = 0.75)
  
  list(volcano = volcano_label, forest = forest_label)
}


# ==============================================================================
#  DATA PREPARATION
# ==============================================================================

# Trait-to-column mapping for the 'combined' data frame
trait_map <- list(
  LDL = list(
    beta = "estimate.ukbb.ldl",
    se   = "std.error.ukbb.ldl",
    p    = "p.value.ukbb.ldl",
    pmr  = "p.value.mr.ldl"
  ),
  TRL = list(
    beta = "estimate.ukbb.trl",
    se   = "std.error.ukbb.trl",
    p    = "p.value.ukbb.trl",
    pmr  = "p.value.mr.trl"
  ),
  LPA = list(
    beta = "estimate.ukbb.lpa",
    se   = "std.error.ukbb.lpa",
    p    = "p.value.ukbb.lpa",
    pmr  = "p.value.mr.lpa"
  )
)


# Extract UKBB data for a specific trait
# 
# @param combined_df Wide format data frame with all traits
# @param trait Trait name: "LDL", "TRL", or "LPA"
# @return Data frame with standardized column names
prepare_ukbb_from_combined <- function(combined_df, trait = c("LDL", "TRL", "LPA")) {
  trait <- match.arg(trait)
  stems <- trait_map[[trait]]
  
  combined_df %>%
    dplyr::select(
      pathway,
      gene,
      beta_obs = all_of(stems$beta),
      se       = all_of(stems$se),
      p_obs    = all_of(stems$p),
      p_mr     = all_of(stems$pmr)
    ) %>%
    arrange(pathway, gene) %>%
    mutate(cohort = "UKBB")
}


# Extract MESA data for a specific trait
# 
# @param mesa_df Long format MESA validation data
# @param trait Trait name: "LDL", "TRL", or "LPA"
# @return Data frame with standardized column names
prepare_mesa_from_long <- function(mesa_df, trait = c("LDL", "TRL", "LPA")) {
  trait <- match.arg(trait)
  mesa_fraction <- if (trait == "LPA") "Lp(a)" else trait
  
  mesa_df %>%
    filter(fraction == mesa_fraction) %>%
    transmute(
      gene,
      beta_obs = estimate.mesa,
      se       = std.error.mesa,
      p_obs    = p.value.mesa,
      p_mr     = NA_real_,  # No MR data for MESA
      cohort   = "MESA"
    )
}


# ==============================================================================
#  BUILD TRAIT-SPECIFIC DATA AND PLOTS
# ==============================================================================

# Manual pathway palette (uses global PATHWAY_COLORS if defined)
.manual_palette <- if (exists("PATHWAY_COLORS", inherits = TRUE)) {
  get("PATHWAY_COLORS", inherits = TRUE)
} else {
  NULL
}


# -----------------------------------------------------------------------------
#  LDL
# -----------------------------------------------------------------------------

ldl_ukbb <- prepare_ukbb_from_combined(combined, trait = "LDL")

ldl_mesa <- prepare_mesa_from_long(mesa_validation_with_ukbb, trait = "LDL") %>%
  left_join(
    combined %>% dplyr::select(gene, pathway) %>% distinct(),
    by = "gene"
  ) %>%
  mutate(pathway = coalesce(pathway, "unknown")) %>%
  dplyr::select(pathway, gene, beta_obs, se, p_obs, p_mr, cohort)

ldl_combined <- bind_rows(ldl_ukbb, ldl_mesa) %>%
  arrange(pathway, gene) %>%
  mutate(x_id = row_number())

ldl_data  <- ldl_ukbb %>% mutate(x_id = row_number())
max_x_ldl <- max(ldl_data$x_id)

p_ldl_volcano <- plot_volcano(
  ldl_data,
  fixed_max_x = max_x_ldl,
  bonf_thresh = BONF_THRESH,
  pathway_palette = .manual_palette
)

f_ldl <- make_forest_plot_posneg_cohorts(
  ldl_combined,
  trait_name = "LDL",
  x_label = "",
  bonf_thresh = BONF_THRESH
)


# -----------------------------------------------------------------------------
#  TRL
# -----------------------------------------------------------------------------

trl_ukbb <- prepare_ukbb_from_combined(combined, trait = "TRL")

trl_mesa <- prepare_mesa_from_long(mesa_validation_with_ukbb, trait = "TRL") %>%
  left_join(
    combined %>% dplyr::select(gene, pathway) %>% distinct(),
    by = "gene"
  ) %>%
  mutate(pathway = coalesce(pathway, "unknown")) %>%
  dplyr::select(pathway, gene, beta_obs, se, p_obs, p_mr, cohort)

trl_combined <- bind_rows(trl_ukbb, trl_mesa) %>%
  arrange(pathway, gene) %>%
  mutate(x_id = row_number())

trl_data <- trl_ukbb %>% mutate(x_id = row_number())

p_trl_volcano <- plot_volcano(
  trl_data,
  bonf_thresh = BONF_THRESH,
  pathway_palette = .manual_palette
)

f_trl <- make_forest_plot_posneg_cohorts(
  trl_combined,
  trait_name = "TRL",
  x_label = "",
  bonf_thresh = BONF_THRESH
)


# -----------------------------------------------------------------------------
#  Lp(a)
# -----------------------------------------------------------------------------

lpa_ukbb <- prepare_ukbb_from_combined(combined, trait = "LPA")

lpa_mesa <- prepare_mesa_from_long(mesa_validation_with_ukbb, trait = "LPA") %>%
  left_join(
    combined %>% dplyr::select(gene, pathway) %>% distinct(),
    by = "gene"
  ) %>%
  mutate(pathway = coalesce(pathway, "unknown")) %>%
  dplyr::select(pathway, gene, beta_obs, se, p_obs, p_mr, cohort)

lpa_combined <- bind_rows(lpa_ukbb, lpa_mesa) %>%
  arrange(pathway, gene) %>%
  mutate(x_id = row_number())

lpa_data <- lpa_ukbb %>% mutate(x_id = row_number())

p_lpa_volcano <- plot_volcano(
  lpa_data,
  cap_at_p = 1e-50,
  small_p_divisor = 2,
  bonf_thresh = BONF_THRESH,
  pathway_palette = .manual_palette
)

f_lpa <- make_forest_plot_posneg_cohorts(
  lpa_combined,
  trait_name = "LPA",
  x_label = "Change in NPX per 1-SD particle increase",
  bonf_thresh = BONF_THRESH
)


# ==============================================================================
#  BUILD FIGURE COMPONENTS
# ==============================================================================

# -----------------------------------------------------------------------------
#  Row 1: Pathway annotation bar (Panel A) + Pathway legend
# -----------------------------------------------------------------------------

# Build pathway legend using LDL palette
palette_ldl <- make_pathway_palette(
  sort(unique(ldl_data$pathway)),
  manual = .manual_palette
)

pathway_legend_plot_raw <- (function(palette) {
  legend_df_pathway <- tibble(
    pathway = factor(sort(names(palette)), levels = sort(names(palette))),
    x = seq_along(palette),
    y = 1
  )
  
  ggplot(legend_df_pathway, aes(x = x, y = y, fill = pathway)) +
    geom_point(shape = 21, color = "black", size = 3.5, show.legend = TRUE) +
    scale_fill_manual(
      name = "Biological class",
      values = palette,
      guide = guide_legend(
        ncol = 2,
        byrow = TRUE,
        title.position = "top",
        label.position = "right"
      )
    ) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 2, 0, 2),
      legend.spacing = unit(2, "pt"),
      legend.box.spacing = unit(2, "pt"),
      legend.title = element_text(size = 8, face = "bold", margin = margin(b = 2)),
      legend.text = element_text(size = 7),
      plot.margin = margin(0, 4, 0, 50)
    )
})(palette_ldl)

pathway_legend_grob <- cowplot::get_legend(pathway_legend_plot_raw)

# Create pathway bar (Panel A) with label
pathway_bar_ldl <- make_pathway_bar(
  ldl_data,
  max_x = max_x_ldl,
  offset = 10000,
  pathway_palette = .manual_palette
) +
  theme(plot.margin = margin(0, 10, 0, 5))

pathway_bar_ldl_labeled <- ggdraw(pathway_bar_ldl) +
  draw_label(
    "A",
    x = -0.28,
    y = 1,
    hjust = 0,
    vjust = 0,
    fontface = "plain",
    size = 14
  )

# Spacer plots
spacer_left  <- ggplot() + theme_void()
spacer_right <- ggplot() + theme_void()

# Assemble Row 1
row1 <- cowplot::plot_grid(
  spacer_left,
  pathway_bar_ldl_labeled,
  pathway_legend_grob,
  spacer_right,
  ncol = 4,
  rel_widths = c(0.4, 1.03, 1.1, 0.1),
  align = "h",
  axis = "tb"
)


# -----------------------------------------------------------------------------
#  Trait rows (volcano + forest plots)
# -----------------------------------------------------------------------------

# Build a trait row with optional column labels
# 
# @param volcano_plot Volcano plot object
# @param forest_plot Forest plot object
# @param trait_name Trait name for side label
# @param include_bc Whether to include "B" and "C" column labels
# @return Combined plot grid
build_trait_row <- function(volcano_plot,
                            forest_plot,
                            trait_name,
                            include_bc = FALSE) {
  
  side_label <- make_trait_side_label(trait_name)
  
  if (include_bc) {
    # Include column labels B and C
    cols <- make_column_labels()
    header <- plot_grid(
      cols$volcano,
      cols$forest,
      ncol = 2,
      rel_widths = c(1, 1)
    )
    body <- plot_grid(
      volcano_plot,
      forest_plot,
      ncol = 2,
      rel_widths = c(1, 1),
      align = "h"
    )
    combined <- plot_grid(
      header,
      body,
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
  } else {
    # No column labels
    combined <- plot_grid(
      volcano_plot,
      forest_plot,
      ncol = 2,
      rel_widths = c(1, 1),
      align = "h"
    )
  }
  
  plot_grid(
    side_label,
    combined,
    ncol = 2,
    rel_widths = c(0.08, 1)
  )
}

row_ldl <- build_trait_row(p_ldl_volcano, f_ldl, "LDL", include_bc = TRUE)
row_trl <- build_trait_row(p_trl_volcano, f_trl, "TRL", include_bc = FALSE)
row_lpa <- build_trait_row(p_lpa_volcano, f_lpa, "Lp(a)", include_bc = FALSE)


# -----------------------------------------------------------------------------
#  Cohort legend row (horizontal, bottom of figure)
# -----------------------------------------------------------------------------

# Build forest plot with legend enabled
f_ldl_with_legend <- make_forest_plot_posneg_cohorts(
  ldl_combined,
  trait_name = "LDL",
  x_label = "",
  bonf_thresh = BONF_THRESH,
  show_legend = TRUE
)

legend_grob <- cowplot::get_legend(f_ldl_with_legend)

# Two-column layout: spacer + legend (for alignment with forest plots)
legend_row <- cowplot::plot_grid(
  ggplot() + theme_void(),  # Transparent spacer
  legend_grob,
  ncol = 2,
  rel_widths = c(0.6, 0.42),
  align = "h"
)


# -----------------------------------------------------------------------------
#  Spacers
# -----------------------------------------------------------------------------

spacer_top <- ggplot() + theme_void() +
  theme(panel.background = element_rect(fill = "transparent", color = NA))

spacer_mid <- ggplot() + theme_void() +
  theme(panel.background = element_rect(fill = "transparent", color = NA))

spacer_bottom <- ggplot() + theme_void()


# ==============================================================================
#  ASSEMBLE FINAL FIGURE
# ==============================================================================

final_plot <- cowplot::plot_grid(
  spacer_top,
  row1,
  spacer_mid,
  row_ldl,
  row_trl,
  row_lpa,
  legend_row,
  spacer_bottom,
  ncol = 1,
  rel_heights = c(0.1, 0.4, 0.05, 1, 1, 1, 0.12, 0.03)
)


# ==============================================================================
#  SAVE OUTPUT
# ==============================================================================

print(final_plot)

ggsave(
  "figure-2.png",
  final_plot,
  width = 10,
  height = 12.5,
  dpi = 600
)