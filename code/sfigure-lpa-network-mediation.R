# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-28
# Last modified: 2025-12-07
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, stringr, forcats, ggplot2, grid, gridExtra,
#                      cowplot, scales, igraph, ggraph, graphlayouts, ggrepel,
#                      DiagrammeR
#
# Goal:
#   Generate a three-panel combined figure for Lp(a)-protein network analysis:
#   Panel A: LPA network visualization (clusters 1, 3, 6) using ggraph
#   Panel B: Forest plot of cluster-protein hazard ratios
#   Panel C: Mediation flowchart showing Lp(a)-P → {PC1 clusters 1/3/6} → CAD
#
# Required input objects:
#   - g_lpa.rda
#   - m_lpa.rda
#   - cox_clusters_long.rda
#   - fit_lpa.rda
#
# Main outputs:
#   - sfigure-lpa-network-mediation.png
#
# ==============================================================================


# ==============================================================================
# LOAD REQUIRED PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(scales)
  library(igraph)
  library(ggraph)
  library(graphlayouts)
  library(ggrepel)
  library(DiagrammeR)
})


# ==============================================================================
# DEFINE COLOR PALETTES
# ==============================================================================

# Cluster color palette
my_cols_named2 <- c(
  `1` = "#3F7DB8FF",
  `2` = "#7FA7C7FF",
  `3` = "#3FB8AFFF",
  `4` = "#7FC7AFFF",
  `5` = "#DAD8A7FF",
  `6` = "#FF9E9DFF"
)

# Mediation diagram colors
col_pale_blue <- "#D6ECFA"
col_c1 <- "#7FA7C7FF"
col_c3 <- "#3FB8AFFF"
col_c6 <- "#FF9E9DFF"


# ==============================================================================
# PANEL A: LP(A) PROTEIN NETWORK VISUALIZATION
# ==============================================================================

pA <- viz_protein_network(
  g_lpa,
  membership      = m_lpa,
  edge_cut        = 0.05,
  layout          = "fr",
  node_size       = "strength",
  label_mode      = "proteins",
  cluster_palette = my_cols_named2,
  seed            = 1,
  hubs_per_comm   = 5
)


# ==============================================================================
# PANEL B: CLUSTER FOREST PLOT (LP(A))
# ==============================================================================

#' Generate forest plot for Lp(a) cluster hazard ratios
#'
#' @param df_all Data frame containing Cox regression results
#' @param highlight_clusters Numeric vector of clusters to highlight
#' @return ggplot object
make_panel_B_lpa_plot <- function(df_all,
                                  highlight_clusters = c(1, 3, 6)) {
  
  # Clean column names (remove backticks)
  nm <- names(df_all)
  nm <- gsub("[`]", "", nm)
  names(df_all) <- nm
  
  # Identify correlation column
  corr_col <- names(df_all)[grepl(
    "Correlation\\s*with\\s*lipoprotein\\s*fraction",
    names(df_all),
    ignore.case = TRUE
  )]
  
  if (length(corr_col) == 0) {
    corr_col <- "Correlation.with.lipoprotein.fraction"
  }
  
  names(df_all)[names(df_all) == corr_col] <- "corr_frac"
  
  # Filter and prepare Lp(a) data
  df_lpa <- df_all %>%
    mutate(Fraction = as.character(Fraction)) %>%
    filter(grepl("^\\s*(LPA|Lp\\(a\\))\\s*$", Fraction, ignore.case = TRUE)) %>%
    mutate(
      Cluster = as.integer(Cluster),
      HR      = exp(logHR),
      LB      = exp(logHR - 1.96 * SE),
      UB      = exp(logHR + 1.96 * SE),
      log_HR  = log(HR),
      log_LB  = log(LB),
      log_UB  = log(UB),
      highlight = if_else(
        Cluster %in% highlight_clusters,
        "highlight",
        "other"
      )
    )
  
  if (nrow(df_lpa) == 0) {
    stop("No LPA rows found in the input table.")
  }
  
  # Order clusters by HR
  df_lpa <- df_lpa %>%
    arrange(HR, Cluster) %>%
    mutate(cluster_fac = fct_inorder(factor(Cluster)))
  
  # Create alternating background stripes
  stripe_df <- df_lpa %>%
    distinct(cluster_fac) %>%
    arrange(cluster_fac) %>%
    mutate(
      y_id   = as.integer(cluster_fac),
      ymin   = y_id - 0.5,
      ymax   = y_id + 0.5,
      fill_bg = if_else(y_id %% 2 == 0, "grey92", "white")
    )
  
  # Define HR axis breaks
  hr_min <- floor(min(df_lpa$LB, na.rm = TRUE) * 100) / 100
  hr_max <- ceiling(max(df_lpa$UB, na.rm = TRUE) * 100) / 100
  hr_min <- max(hr_min, 0.2)
  hr_max <- min(hr_max, 5)
  
  candidate_breaks <- unique(sort(c(
    0.25, 0.33, 0.5, 0.67, 0.75, 1, 1.1, 1.25, 1.4, 1.5, 2, 3, 4, 5
  )))
  
  hr_breaks <- candidate_breaks[
    candidate_breaks >= hr_min & candidate_breaks <= hr_max
  ]
  
  log_breaks <- log(hr_breaks)
  lab_breaks <- sprintf("%.2g", hr_breaks)
  
  # Define colors for highlighting
  cols <- c(
    other     = "steelblue",
    highlight = "firebrick"
  )
  
  # Create forest plot
  ggplot() +
    # Background stripes (white)
    geom_rect(
      data = dplyr::filter(stripe_df, fill_bg == "white"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      fill = "white",
      color = NA
    ) +
    # Background stripes (grey)
    geom_rect(
      data = dplyr::filter(stripe_df, fill_bg == "grey92"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      fill = "grey92",
      color = NA
    ) +
    # Null effect reference line
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    # Confidence intervals
    geom_errorbarh(
      data = df_lpa,
      aes(
        y = cluster_fac,
        xmin = log_LB,
        xmax = log_UB,
        color = highlight
      ),
      height = 0.3,
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    # Point estimates
    geom_point(
      data = df_lpa,
      aes(y = cluster_fac, x = log_HR, color = highlight),
      shape = 15,
      size = 2,
      stroke = 0.6,
      show.legend = FALSE
    ) +
    # Scales and labels
    scale_x_continuous(
      limits = log(c(hr_min - 0.02, hr_max + 0.02)),
      breaks = log_breaks,
      labels = lab_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = cols, guide = "none") +
    labs(
      x = "Hazard ratio (per 1-SD)",
      y = "Protein cluster no."
    ) +
    # Theme customization
    theme_linedraw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10),
      axis.title.x       = element_text(size = 11),
      axis.title.y       = element_text(size = 10),
      legend.position    = "none",
      plot.title         = element_blank(),
      plot.margin        = margin(4, 6, 4, 6),
      panel.background   = element_rect(fill = NA, color = NA),
      plot.background    = element_rect(fill = NA, color = NA)
    )
}

# Generate Panel B
pB <- make_panel_B_lpa_plot(
  cox_clusters_long,
  highlight_clusters = c(1, 3, 6)
)


# ==============================================================================
# PANEL C: CAUSAL MEDIATION FLOWCHART
# ==============================================================================

# Define box dimensions and positions
box_w <- 2.2
box_h <- 1.1

x_lpa <- 0
x_cad <- 8
x_mid <- (x_lpa + x_cad) / 2

y_base <- 0
y_c1 <- 2.5
y_c3 <- 4.5
y_c6 <- 6.5
y_top <- y_base + box_h / 2

# Define nodes
nodes <- data.frame(
  name  = c("Lp(a)", "CAD", "C1", "C3", "C6"),
  x     = c(x_lpa, x_cad, x_mid, x_mid, x_mid),
  y     = c(y_base, y_base, y_c1, y_c3, y_c6),
  fill  = c(col_pale_blue, col_pale_blue, col_c1, col_c3, col_c6),
  label = c("Lp(a)", "CAD", "Cluster 1", "Cluster 3", "Cluster 6")
)

# Calculate rectangle coordinates
rects <- within(nodes, {
  xmin <- x - box_w / 2
  xmax <- x + box_w / 2
  ymin <- y - box_h / 2
  ymax <- y + box_h / 2
})

# Calculate edge positions
lpa_right_x  <- x_lpa + box_w / 2
cad_left_x   <- x_cad - box_w / 2
x_left_edge  <- x_mid - box_w / 2
x_right_edge <- x_mid + box_w / 2

# Define arrows from Lp(a) to clusters (upward)
arrows_up <- data.frame(
  x = x_lpa,
  y = y_top,
  xend = x_left_edge,
  yend = c(y_c1, y_c3, y_c6),
  arrow = TRUE
)

# Define arrows from clusters to CAD (downward)
arrows_down <- data.frame(
  x = x_right_edge,
  y = c(y_c1, y_c3, y_c6),
  xend = x_cad,
  yend = y_top,
  arrow = c(FALSE, TRUE, FALSE)
)

# Define direct arrow from Lp(a) to CAD
direct_arrow <- data.frame(
  x = lpa_right_x,
  y = y_base,
  xend = cad_left_x,
  yend = y_base
)

#' Add curved arrows to a ggplot object
#'
#' @param p ggplot object
#' @param df Data frame with arrow coordinates
#' @return Modified ggplot object
add_curve_rows <- function(p, df) {
  for (i in seq_len(nrow(df))) {
    p <- p + geom_curve(
      data = df[i, ],
      aes(x = x, y = y, xend = xend, yend = yend),
      curvature = 0,
      linewidth = 0.8,
      color = "black",
      arrow = if (isTRUE(df$arrow[i])) {
        arrow(type = "closed", length = unit(6, "pt"))
      } else {
        NULL
      }
    )
  }
  return(p)
}

# Initialize plot with boxes
pC <- ggplot() +
  # Lp(a) and CAD boxes
  geom_rect(
    data = rects[rects$name %in% c("Lp(a)", "CAD"), ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = col_pale_blue,
    color = NA
  ) +
  # Cluster boxes
  geom_rect(
    data = rects[rects$name %in% c("C1", "C3", "C6"), ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
    color = NA
  ) +
  scale_fill_manual(
    values = c(C1 = col_c1, C3 = col_c3, C6 = col_c6),
    guide = "none"
  ) +
  # Box labels
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label),
    size = 4
  ) +
  # Direct effect arrow
  geom_segment(
    data = direct_arrow,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    linewidth = 0.8,
    arrow = arrow(type = "closed", length = unit(4, "pt"))
  ) +
  # Coordinate system
  coord_equal(
    xlim = c(-2, 10),
    ylim = c(-1, 8),
    expand = FALSE
  ) +
  theme_void()

# Add curved arrows
pC <- add_curve_rows(pC, arrows_up)
pC <- add_curve_rows(pC, arrows_down)

#' Build effect estimate labels for mediation diagram
#'
#' @param fit Mediation analysis results object
#' @return List of formatted label strings
build_effect_labels <- function(fit) {
  pe <- fit$effect.pe
  lo <- fit$effect.ci.low
  hi <- fit$effect.ci.high
  
  fnum <- function(x) sprintf("%.2f", x)
  fpct <- function(x) sprintf("%.0f", x * 100)
  
  list(
    total_label = paste0(
      "Total effect: ", fnum(pe["Rte"]),
      " (95%CI: ", fnum(lo["Rte"]), "-", fnum(hi["Rte"]), ")"
    ),
    direct_label = paste0(
      "Direct effect: ", fnum(pe["Rcde"]),
      " (95%CI: ", fnum(lo["Rcde"]), "-", fnum(hi["Rcde"]), ")"
    ),
    indirect_label = paste0(
      "Indirect effect: ", fnum(pe["Rpnie"]),
      " (95%CI: ", fnum(lo["Rpnie"]), "-", fnum(hi["Rpnie"]), ")"
    ),
    prop_label = paste0(
      "Proportion mediated: ", fpct(pe["pm"]),
      "% (95%CI: ", fpct(lo["pm"]), "-", fpct(hi["pm"]), "%)"
    )
  )
}

# Generate effect labels
labs <- build_effect_labels(fit_lpa)

# Define label positions
x_center <- (x_lpa + x_cad) / 2
y_total <- y_base - 0.30
y_direct <- y_base + 0.30
y_top_cluster <- max(y_c1, y_c3, y_c6) + box_h / 2
y_indirect_1 <- y_top_cluster + 0.70
y_indirect_2 <- y_top_cluster + 0.25

# Add effect annotations
pC <- pC +
  annotate(
    "text",
    x = x_center,
    y = y_total,
    label = labs$total_label,
    fontface = "bold",
    size = 3.5,
    vjust = 1
  ) +
  annotate(
    "text",
    x = x_center,
    y = y_direct,
    label = labs$direct_label,
    size = 3.5,
    vjust = 0
  ) +
  annotate(
    "text",
    x = x_center,
    y = y_indirect_1,
    label = labs$indirect_label,
    fontface = "bold",
    size = 3.5,
    vjust = 0
  ) +
  annotate(
    "text",
    x = x_center,
    y = y_indirect_2,
    label = labs$prop_label,
    size = 3.5,
    vjust = 0
  )


# ==============================================================================
# COMBINE PANELS INTO FINAL FIGURE
# ==============================================================================

# Convert plots to grobs
pA_grob <- ggplotGrob(pA$plot)
pB_grob <- ggplotGrob(pB)
pC_grob <- ggplotGrob(pC)

#' Scale grob within grid cell
#'
#' @param grob Grid object to scale
#' @param width,height Relative dimensions
#' @param x,y Position within cell
#' @param just Justification
#' @return Modified grob tree
scale_inside_cell <- function(grob,
                              width = 1,
                              height = 1,
                              x = 0.5,
                              y = 0.5,
                              just = c("center", "center")) {
  grid::grobTree(
    grid::editGrob(
      grob,
      vp = grid::viewport(
        width = width,
        height = height,
        x = x,
        y = y,
        just = just
      )
    )
  )
}

# Scale individual panels
pA_scaled <- scale_inside_cell(pA_grob, width = 1.0, height = 1.0)
pB_scaled <- scale_inside_cell(pB_grob, width = 0.8, height = 0.8)
pC_scaled <- scale_inside_cell(pC_grob, width = 1.0, height = 1.0)

# Arrange top row (Panels A and B)
top_row <- gridExtra::arrangeGrob(
  pA_scaled,
  pB_scaled,
  ncol = 2,
  widths = c(1.0, 1.1),
  heights = 1
)

# Arrange bottom row (Panel C with spacer)
spacer <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
bottom_row <- gridExtra::arrangeGrob(
  pC_scaled,
  spacer,
  ncol = 2,
  widths = c(1.0, 1.1),
  heights = 1
)

# Combine rows into final figure
final_fig <- gridExtra::arrangeGrob(
  top_row,
  bottom_row,
  ncol = 1,
  heights = c(1.0, 0.8)
)

# Add panel labels (A, B, C)
labeled_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(final_fig) +
  cowplot::draw_plot_label(
    label = c("A", "B", "C"),
    x = c(0.02, 0.55, 0.02),
    y = c(0.98, 0.98, 0.46),
    size = 14,
    fontface = "plain"
  )

# Display figure
print(labeled_plot)


# ==============================================================================
# SAVE FIGURE
# ==============================================================================

# Save high-resolution figure
ggsave(
  filename = "sfigure-lpa-network-mediation.png",
  plot = labeled_plot,
  width = 10,
  height = 8,
  dpi = 1200,
  bg = "white"
)
