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
#   Required packages: dplyr, stringr, forcats, ggplot2, grid, gridExtra,
#                      cowplot, scales, igraph, ggraph, graphlayouts, 
#                      ggrepel, DiagrammeR
#
# Goal:
#   Generate a three-panel composite figure combining:
#   Panel A: TRL protein network visualization (clusters 8, 13, 34)
#   Panel B: Forest plot of cluster-protein hazard ratios (UKBB + MESA)
#   Panel C: Causal mediation flowchart (TRL → PC1 clusters → CAD)
#
# Required input objects:
#   - g_trl.rda
#   - m_trl.rda
#   - cox_clusters_long: data.frame (Cox regression results)
#   - fit_trlp.rda
#
# Outputs:
#   - figure-5.png 
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
# PANEL A: TRL PROTEIN NETWORK VISUALIZATION
# ==============================================================================

# Define cluster color palette
my_cols_named3 <- c(
  `8`  = "#7FA7C7FF",
  `13` = "#3FB8AFFF",
  `34` = "#FF9E9DFF"
)

# Set correlation threshold
rho_cut <- 0.1

# Build TRL subgraph plot
# Note: Requires pre-defined helper functions:
#   - viz_protein_network()
#   - plot_selected_clusters()
#   - save_transparent()
p_trl_sel <- plot_selected_clusters(
  g_trl, 
  m_trl, 
  clusters = c(8, 13, 34),
  edge_cut = rho_cut,
  layout = "stress",
  node_size = "strength",
  label_mode = "proteins",
  cluster_palette = my_cols_named3,
  seed = 1
)

# Optional manual adjustment for CCL7_rz and CCL8_rz positions
idx_ccl7 <- which(p_trl_sel$data$name == "CCL7_rz")
idx_ccl8 <- which(p_trl_sel$data$name == "CCL8_rz")
p_trl_sel$data$x[idx_ccl7] <- p_trl_sel$data$x[idx_ccl7] + 7
p_trl_sel$data$y[idx_ccl7] <- p_trl_sel$data$y[idx_ccl7] + 6
p_trl_sel$data$x[idx_ccl8] <- p_trl_sel$data$x[idx_ccl8] + 7
p_trl_sel$data$y[idx_ccl8] <- p_trl_sel$data$y[idx_ccl8] + 6

# Finalize Panel A
pA <- p_trl_sel + 
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 12))


# ==============================================================================
# PANEL B: FOREST PLOT OF CLUSTER HAZARD RATIOS
# ==============================================================================

#' Create TRL-specific forest plot
#'
#' @param df_all Data frame containing Cox regression results
#' @param highlight_clusters Numeric vector of cluster IDs to highlight
#' @return ggplot object
make_panel_B_trl_plot <- function(df_all,
                                  highlight_clusters = c(8, 13, 34)) {
  
  # Clean column names
  nm <- names(df_all)
  nm <- gsub("[`]", "", nm)
  names(df_all) <- nm
  
  # Identify correlation column
  corr_col <- names(df_all)[
    grepl("Correlation\\s*with\\s*lipoprotein\\s*fraction",
          names(df_all), 
          ignore.case = TRUE)
  ]
  if (length(corr_col) == 0) {
    corr_col <- "Correlation.with.lipoprotein.fraction"
  }
  names(df_all)[names(df_all) == corr_col] <- "corr_frac"
  
  # Filter and prepare TRL data
  df_trl <- df_all %>%
    mutate(Fraction = as.character(Fraction)) %>%
    filter(grepl("^\\s*TRL\\s*$", Fraction, ignore.case = TRUE)) %>%
    mutate(
      Cluster = as.integer(Cluster),
      HR  = exp(logHR),
      LB  = exp(logHR - 1.96 * SE),
      UB  = exp(logHR + 1.96 * SE),
      log_HR = log(HR), 
      log_LB = log(LB), 
      log_UB = log(UB),
      highlight = if_else(Cluster %in% highlight_clusters, 
                          "highlight", 
                          "other")
    )
  
  if (nrow(df_trl) == 0) {
    stop("No TRL rows found in the input table.")
  }
  
  # Order by increasing hazard ratio
  df_trl <- df_trl %>% 
    arrange(HR, Cluster) %>% 
    mutate(cluster_fac = fct_inorder(factor(Cluster)))
  
  # Create alternating background stripes
  stripe_df <- df_trl %>%
    distinct(cluster_fac) %>%
    arrange(cluster_fac) %>%
    mutate(
      y_id = as.integer(cluster_fac),
      ymin = y_id - 0.5, 
      ymax = y_id + 0.5,
      fill_bg = if_else(y_id %% 2 == 0, "grey92", "white")
    )
  
  # Define x-axis breaks
  hr_min <- floor(min(df_trl$LB, na.rm = TRUE) * 100) / 100
  hr_max <- ceiling(max(df_trl$UB, na.rm = TRUE) * 100) / 100
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
  
  # Define colors
  cols <- c(
    other = "steelblue", 
    highlight = "firebrick"
  )
  
  # Build plot
  ggplot() +
    # Background stripes
    geom_rect(
      data = filter(stripe_df, fill_bg == "white"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      inherit.aes = FALSE, 
      fill = "white", 
      color = NA
    ) +
    geom_rect(
      data = filter(stripe_df, fill_bg == "grey92"),
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      inherit.aes = FALSE, 
      fill = "grey92", 
      color = NA
    ) +
    # Null line
    geom_vline(
      xintercept = 0, 
      linetype = "dashed", 
      color = "grey50"
    ) +
    # Confidence intervals
    geom_errorbarh(
      data = df_trl,
      aes(y = cluster_fac, xmin = log_LB, xmax = log_UB, color = highlight),
      height = 0.3, 
      linewidth = 0.6, 
      show.legend = FALSE
    ) +
    # Point estimates
    geom_point(
      data = df_trl,
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
    # Theme
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

# Build Panel B
pB <- make_panel_B_trl_plot(cox_clusters_long)


# ==============================================================================
# PANEL C: CAUSAL MEDIATION FLOWCHART
# ==============================================================================

# --- Layout constants ---
box_w <- 2.2
box_h <- 1.1
x_trl <- 0
x_cad <- 8
x_mid <- (x_trl + x_cad) / 2
y_base <- 0
y_c8 <- 2.5
y_c13 <- 4.5
y_c34 <- 6.5
y_top <- y_base + box_h / 2
x_top_trl <- x_trl
x_top_cad <- x_cad

# --- Color palette ---
col_pale_blue <- "#D6ECFA"
col_c8  <- "#7FA7C7FF"
col_c13 <- "#3FB8AFFF"
col_c34 <- "#FF9E9DFF"

# --- Define nodes ---
nodes <- data.frame(
  name  = c("TRL", "CAD", "C8", "C13", "C34"),
  x     = c(x_trl, x_cad, x_mid, x_mid, x_mid),
  y     = c(y_base, y_base, y_c8, y_c13, y_c34),
  fill  = c(col_pale_blue, col_pale_blue, col_c8, col_c13, col_c34),
  label = c("TRL", "CAD", "Cluster 8", "Cluster 13", "Cluster 34")
)

# Define rectangles
rects <- within(nodes, {
  xmin <- x - box_w / 2
  xmax <- x + box_w / 2
  ymin <- y - box_h / 2
  ymax <- y + box_h / 2
})

# Edge coordinates
trl_right_x  <- x_trl + box_w / 2
cad_left_x   <- x_cad - box_w / 2
x_left_edge  <- x_mid - box_w / 2
x_right_edge <- x_mid + box_w / 2

# --- Define arrows ---
# Arrows from TRL to clusters (upward)
arrows_up <- data.frame(
  x = x_top_trl, 
  y = y_top,
  xend = x_left_edge, 
  yend = c(y_c8, y_c13, y_c34),
  curvature = c(0, 0, 0), 
  arrow = TRUE
)

# Arrows from clusters to CAD (downward)
arrows_down <- data.frame(
  x = x_right_edge, 
  y = c(y_c8, y_c13, y_c34),
  xend = x_top_cad, 
  yend = y_top,
  curvature = c(0, 0, 0),
  arrow = c(FALSE, TRUE, FALSE)
)

# Direct arrow from TRL to CAD
direct_arrow <- data.frame(
  x = trl_right_x, 
  y = y_base, 
  xend = cad_left_x, 
  yend = y_base
)

# --- Helper function to add curved arrows ---
add_curve_rows <- function(p, df, angle = 70, ncp = 12, lw = 0.8) {
  for (i in seq_len(nrow(df))) {
    p <- p + geom_curve(
      data = df[i, ],
      aes(x = x, y = y, xend = xend, yend = yend),
      curvature = df$curvature[i],
      angle = angle, 
      ncp = ncp,
      linewidth = lw, 
      color = "black",
      arrow = if (isTRUE(df$arrow[i])) {
        arrow(type = "closed", length = unit(6, "pt"))
      } else {
        NULL
      }
    )
  }
  p
}

# --- Build base structure ---
pC <- ggplot() +
  # TRL and CAD boxes
  geom_rect(
    data = rects[rects$name %in% c("TRL", "CAD"), ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = col_pale_blue, 
    color = NA
  ) +
  # Cluster boxes
  geom_rect(
    data = rects[rects$name %in% c("C8", "C13", "C34"), ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
    color = NA
  ) +
  scale_fill_manual(
    values = c(C8 = col_c8, C13 = col_c13, C34 = col_c34), 
    guide = "none"
  ) +
  # Box labels
  geom_text(
    data = nodes, 
    aes(x = x, y = y, label = label), 
    size = 4
  ) +
  # Direct arrow
  geom_segment(
    data = direct_arrow,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black", 
    linewidth = 0.8,
    arrow = arrow(type = "closed", length = unit(4, "pt"))
  ) +
  coord_equal(xlim = c(-2, 10), ylim = c(-1, 8), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Add curved arrows
pC <- add_curve_rows(pC, arrows_up)
pC <- add_curve_rows(pC, arrows_down)

# --- Add mediation effect labels ---

#' Build effect labels from mediation analysis results
#'
#' @param fit Mediation analysis results object
#' @param total_key Key name for total effect
#' @param direct_key Key name for direct effect
#' @param indirect_key Key name for indirect effect
#' @param digits Number of decimal places for effects
#' @param pct_digits Number of decimal places for percentages
#' @return List of formatted label strings
build_effect_labels <- function(fit,
                                total_key   = "Rte",
                                direct_key  = "Rcde",
                                indirect_key = "Rpnie",
                                digits = 2, 
                                pct_digits = 0) {
  
  pe  <- fit$effect.pe
  lo  <- fit$effect.ci.low
  hi  <- fit$effect.ci.high
  
  fnum <- function(x) sprintf(paste0("%.", digits, "f"), x)
  fpct <- function(x) sprintf(paste0("%.", pct_digits, "f"), x * 100)
  
  getv <- function(v, key) {
    if (!is.null(names(v)) && key %in% names(v)) {
      v[[key]]
    } else {
      NA_real_
    }
  }
  
  te  <- getv(pe, total_key)
  te_lo  <- getv(lo, total_key)
  te_hi  <- getv(hi, total_key)
  
  de  <- getv(pe, direct_key)
  de_lo  <- getv(lo, direct_key)
  de_hi  <- getv(hi, direct_key)
  
  ie  <- getv(pe, indirect_key)
  ie_lo <- getv(lo, indirect_key)
  ie_hi <- getv(hi, indirect_key)
  
  pm  <- getv(pe, "pm")
  pm_lo <- getv(lo, "pm")
  pm_hi <- getv(hi, "pm")
  
  list(
    total_label = paste0(
      "Total effect: ", fnum(te), 
      " (95%CI: ", fnum(te_lo), "-", fnum(te_hi), ")"
    ),
    direct_label = paste0(
      "Direct effect: ", fnum(de), 
      " (95%CI: ", fnum(de_lo), "-", fnum(de_hi), ")"
    ),
    indirect_label = paste0(
      "Indirect effect: ", fnum(ie), 
      " (95%CI: ", fnum(ie_lo), "-", fnum(ie_hi), ")"
    ),
    prop_label = paste0(
      "Proportion mediated: ", fpct(pm), 
      "% (95%CI: ", fpct(pm_lo), "-", fpct(pm_hi), "%)"
    )
  )
}

# Generate labels
labs <- build_effect_labels(fit_trlp)

# Define label positions
x_center <- (x_trl + x_cad) / 2
y_total <- y_base - 0.30
y_direct <- y_base + 0.30
y_top_cluster <- max(y_c8, y_c13, y_c34) + box_h / 2
y_indirect_1 <- y_top_cluster + 0.70
y_indirect_2 <- y_top_cluster + 0.25

# Add labels to plot
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
pA_grob <- ggplotGrob(pA)
pB_grob <- ggplotGrob(pB)
pC_grob <- ggplotGrob(pC)

# Shift Panel A upward (optional aesthetic adjustment)
pA_grob_shifted <- grid::editGrob(
  grob = pA_grob,
  vp   = grid::viewport(y = 0.20, height = 1.10, just = "bottom")
)

# --- Helper function to scale grobs within cells ---
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

# Scale Panel A and Panel C within their cells
pA_scaled <- scale_inside_cell(
  pA_grob_shifted, 
  width = 1.00, 
  height = 0.95, 
  x = 0.5, 
  y = 0.52
)

pC_scaled <- scale_inside_cell(
  pC_grob, 
  width = 1.5, 
  height = 1.1, 
  x = 0.5, 
  y = 0.55
)

# Arrange left column (A + C stacked)
left_col <- gridExtra::arrangeGrob(
  pA_scaled, 
  pC_scaled, 
  ncol = 1, 
  heights = c(1, 1)
)

# Right column (B)
right_col <- pB_grob

# Combine columns
final_two <- gridExtra::arrangeGrob(
  left_col, 
  right_col, 
  ncol = 2, 
  widths = c(1.8, 1)
)

# Add panel labels (A, B, C)
labeled_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(final_two) +
  cowplot::draw_plot_label(
    label = c("A", "C", "B"),
    x = c(0.02, 0.02, 0.62),
    y = c(0.99, 0.56, 0.99),
    size = 14,
    fontface = "plain"
  )

# Display plot
print(labeled_plot)


# ==============================================================================
# SAVE FINAL FIGURE
# ==============================================================================

# Save figure
ggsave(
  filename = "figure-5.png",
  plot = labeled_plot,
  width = 10,
  height = 8,
  dpi = 1200,
  bg = "white"
)

