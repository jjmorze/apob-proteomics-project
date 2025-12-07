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
#   Required packages: huge, igraph, tidygraph, ggraph, dplyr, tibble, stringr,
#                      ggplot2, cowplot, svglite (optional)
#
# Goal:
#   Perform protein network analysis using graphical LASSO (GLASSO) to derive
#   partial correlation networks, identify protein communities via walktrap
#   clustering, compute multi-protein scores (MPS) as PC1 per community, and
#   visualize networks with sign-colored edges and community-based layouts.
#
#   Key steps:
#   1) Fit GLASSO to obtain precision matrix and partial correlations
#   2) Build graph from partial correlations with |ρ| threshold
#   3) Cluster proteins using walktrap algorithm
#   4) Derive PC1 scores per cluster (MPS)
#   5) Align scores to anchor lipoproteins (TRL, LDL, Lp(a))
#   6) Visualize networks with customizable layouts and palettes
#
# Required input data:
#   - trl_protein_data.rda
#   - lpa_protein_data.rda
#   - dataset.rda
#
# Main outputs:
#   - fit_trl, fit_lpa: GLASSO fitted models
#   - rev_scores_trl.rda
#   - rev_scores_lpa.rda
#   - membership_trl.rda
#   - membership_lpa.rda
#   - g_trl, g_lpa: igraph network objects
#   - m_trl, m_lpa: Cluster membership vectors
#   - scores_trl, scores_lpa: MPS (multi-protein scores)
#   - Network visualization plots (PNG format)
#
# ==============================================================================


# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(huge)         # GLASSO fitting
  library(igraph)       # Graph operations
  library(tidygraph)    # Tidy graph manipulation
  library(ggraph)       # Graph visualization
  library(dplyr)        # Data manipulation
  library(tibble)       # Modern data frames
  library(stringr)      # String operations
  library(ggplot2)      # Plotting
  library(cowplot)      # Legend extraction
})


# ==============================================================================
#  UTILITY FUNCTIONS
# ==============================================================================

# Null-coalescing operator
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# ==============================================================================
#  STEP 1: FIT GLASSO AND COMPUTE PARTIAL CORRELATIONS
# ==============================================================================

#' Fit Graphical LASSO and Compute Partial Correlation Matrix
#'
#' @param X Numeric matrix of protein data (samples x proteins)
#' @param lambda Regularization parameter for GLASSO
#' @return List containing:
#'   - omega: Precision matrix (inverse covariance)
#'   - pcor: Partial correlation matrix
fit_glasso_pcor <- function(X, lambda) {
  X <- as.matrix(X)
  
  # Fit GLASSO model
  fit <- huge(
    X,
    method = "glasso",
    lambda = lambda,
    scr = FALSE,
    verbose = FALSE
  )
  
  # Extract precision matrix
  omega <- fit$icov[[1]]
  
  # Convert to partial correlations
  d <- sqrt(diag(omega))
  pcor <- -omega / (outer(d, d, `*`))
  diag(pcor) <- 1
  
  list(omega = omega, pcor = pcor)
}


# ==============================================================================
#  STEP 2: BUILD GRAPH FROM PARTIAL CORRELATIONS
# ==============================================================================

#' Build igraph Network from Partial Correlation Matrix
#'
#' @param pcor Partial correlation matrix
#' @param rho_cut Threshold for |ρ| to include edges (default: 0.10)
#' @param node_names Optional character vector of node names
#' @return igraph object with edge attributes (pcor, weight, sign)
graph_from_pcor <- function(pcor, rho_cut = 0.10, node_names = NULL) {
  P <- as.matrix(pcor)
  storage.mode(P) <- "numeric"
  stopifnot(nrow(P) == ncol(P))
  diag(P) <- 0
  
  # Assign node names
  vnames <- if (is.null(node_names)) {
    as.character(seq_len(nrow(P)))
  } else {
    as.character(node_names)
  }
  
  # Filter edges by |ρ| threshold
  keep <- abs(P) >= rho_cut
  keep[is.na(keep)] <- FALSE
  diag(keep) <- FALSE
  
  # Extract upper triangle indices
  idx <- which(keep, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  
  # Handle empty graph
  if (nrow(idx) == 0) {
    g <- make_empty_graph(n = length(vnames))
    V(g)$name <- vnames
    return(g)
  }
  
  # Create edge data frame
  df <- data.frame(
    from = vnames[idx[, 1]],
    to = vnames[idx[, 2]],
    pcor = P[idx],
    abs_weight = abs(P[idx]),
    sign = ifelse(P[idx] < 0, "negative", "positive"),
    stringsAsFactors = FALSE
  )
  
  # Build graph
  g <- graph_from_data_frame(
    df,
    directed = FALSE,
    vertices = tibble(name = vnames)
  )
  
  # Add edge attributes
  E(g)$pcor <- df$pcor
  E(g)$weight <- pmax(.Machine$double.eps, df$abs_weight)  # Safe for layouts
  E(g)$sign <- factor(df$sign, levels = c("negative", "positive"))
  
  g
}


# ==============================================================================
#  STEP 3: CLUSTER PROTEINS USING WALKTRAP ALGORITHM
# ==============================================================================

#' Perform Walktrap Clustering on Protein Network
#'
#' @param g igraph object
#' @param steps Number of random walk steps (default: 4)
#' @return Named integer vector of cluster membership
walktrap_clusters <- function(g, steps = 4L) {
  if (vcount(g) == 0) {
    return(setNames(integer(0), character(0)))
  }
  if (ecount(g) == 0) {
    return(setNames(rep(1L, vcount(g)), V(g)$name))
  }
  
  cw <- cluster_walktrap(g, weights = E(g)$weight, steps = steps)
  setNames(membership(cw), V(g)$name)
}


# ==============================================================================
#  STEP 4: COMPUTE PC1 SCORES PER CLUSTER (MPS)
# ==============================================================================

#' Compute First Principal Component per Cluster (Multi-Protein Scores)
#'
#' @param X Numeric matrix of protein data (samples x proteins)
#' @param membership Named integer vector of cluster assignments
#' @param ids Optional vector of sample IDs
#' @param prefix Prefix for score column names (default: "C")
#' @return Tibble with id column and cluster PC1 scores (C1_PC1, C2_PC1, ...)
cluster_pc1_scores <- function(X, membership, ids = NULL, prefix = "C") {
  X <- as.matrix(X)
  ids <- ids %||% rownames(X) %||% seq_len(nrow(X))
  stopifnot(nrow(X) == length(ids))
  
  # Map membership to columns
  cl_df <- tibble(
    node = colnames(X) %||% paste0("V", seq_len(ncol(X))),
    cluster = as.integer(membership[match(colnames(X), names(membership))])
  )
  clusters <- sort(unique(cl_df$cluster))
  
  # Compute PC1 per cluster
  scores <- lapply(clusters, function(cl) {
    cols <- cl_df$node[cl_df$cluster == cl]
    if (length(cols) == 0L) return(NULL)
    if (length(cols) == 1L) {
      as.numeric(scale(X[, cols, drop = TRUE]))
    } else {
      prcomp(X[, cols, drop = FALSE], center = FALSE, scale. = FALSE)$x[, 1]
    }
  })
  
  # Filter out NULL results
  keep <- !vapply(scores, is.null, logical(1))
  scores <- scores[keep]
  clusters <- clusters[keep]
  
  # Combine into matrix
  scores_mat <- if (length(scores)) {
    do.call(cbind, scores)
  } else {
    matrix(nrow = nrow(X), ncol = 0)
  }
  colnames(scores_mat) <- paste0(prefix, clusters, "_PC1")
  
  tibble(id = ids) %>% bind_cols(as_tibble(scores_mat))
}


# ==============================================================================
#  STEP 5: VISUALIZE PROTEIN NETWORK
# ==============================================================================

#' Visualize Protein Network with Community Detection
#'
#' @param g igraph object
#' @param membership Named integer vector of cluster assignments (optional)
#' @param edge_cut Minimum |ρ| for plotting edges (default: 0.10)
#' @param layout Layout algorithm: "stress", "drl", "fr", "kk", "lgl", "circle"
#' @param node_size Node sizing: "strength", "degree", "constant"
#' @param point_size Base point size for nodes (default: 2)
#' @param label_mode Label display: "proteins", "hubs", "none"
#' @param hubs_per_comm Number of hub labels per community (default: 3)
#' @param label_strip_regex Regex to strip from labels (default: "_rz$")
#' @param max_overlaps Max label overlaps allowed (default: 50)
#' @param cluster_palette Manual cluster colors (vector or named vector)
#' @param seed Random seed for layout reproducibility (default: 1)
#' @return List with tbl_graph, layout_df, and ggplot object
viz_protein_network <- function(
    g,
    membership = NULL,
    edge_cut = 0.10,
    layout = c("stress", "drl", "fr", "kk", "lgl", "circle"),
    node_size = c("strength", "degree", "constant"),
    point_size = 2,
    label_mode = c("proteins", "hubs", "none"),
    hubs_per_comm = 3,
    label_strip_regex = "_rz$",
    max_overlaps = 50,
    cluster_palette = NULL,
    seed = 1
) {
  
  # Match arguments
  layout <- match.arg(layout)
  node_size <- match.arg(node_size)
  label_mode <- match.arg(label_mode)
  set.seed(seed)
  
  # Validate inputs
  stopifnot(inherits(g, "igraph"))
  if (vcount(g) == 0) stop("Graph has no nodes.")
  if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
  if (is.null(E(g)$weight)) stop("Graph must have numeric edge weights (|ρ|).")
  if (is.null(E(g)$sign)) stop("Edges must have 'sign' attribute.")
  
  # Prune edges and remove isolates for plotting
  g2 <- delete_edges(g, which(E(g)$weight < edge_cut))
  g2 <- delete_vertices(g2, which(degree(g2) == 0))
  
  if (gorder(g2) == 0L || gsize(g2) == 0L) {
    stop("Graph empty after pruning; lower 'edge_cut'.")
  }
  
  # Map membership to nodes (by name)
  if (!is.null(membership)) {
    memb_vec <- as.integer(membership)
    names(memb_vec) <- names(membership)
    common <- intersect(names(memb_vec), V(g2)$name)
    
    if (!length(common)) {
      stop("No overlap between membership names and vertex names (after pruning).")
    }
    
    V(g2)$community <- NA_integer_
    V(g2)$community[match(common, V(g2)$name)] <-
      memb_vec[match(common, names(memb_vec))]
  }
  
  # Assign connected components to NA communities
  if (is.null(V(g2)$community) || any(is.na(V(g2)$community))) {
    comps <- components(g2)$membership
    if (is.null(V(g2)$community)) {
      V(g2)$community <- comps
    } else {
      na_nodes <- which(is.na(V(g2)$community))
      if (length(na_nodes)) V(g2)$community[na_nodes] <- comps[na_nodes]
    }
  }
  
  # Compute node metrics
  V(g2)$degree <- igraph::degree(g2, mode = "all")
  V(g2)$strength <- igraph::strength(g2, mode = "all", weights = E(g2)$weight)
  
  # Select layout
  tg <- as_tbl_graph(g2)
  lay_name <- switch(
    layout,
    stress = if (requireNamespace("graphlayouts", quietly = TRUE)) "stress" else "drl",
    drl = "drl",
    fr = "fr",
    kk = "kk",
    lgl = "lgl",
    circle = "linear"
  )
  glay <- ggraph::create_layout(tg, layout = lay_name, circular = (layout == "circle"))
  
  # Scale node size
  size_src <- switch(
    node_size,
    strength = glay$strength,
    degree = glay$degree,
    constant = rep(point_size, nrow(glay))
  )
  if (!is.numeric(size_src)) size_src <- rep(point_size, nrow(glay))
  
  rng <- range(size_src, na.rm = TRUE)
  glay$._size <- if (diff(rng) > 0) {
    scales::rescale(size_src, to = c(1.8, 6), from = rng)
  } else {
    rep(3, nrow(glay))
  }
  
  # Strip regex from labels for display
  disp_names <- glay$name
  if (!is.null(label_strip_regex) && nzchar(label_strip_regex)) {
    disp_names <- gsub(label_strip_regex, "", disp_names, perl = TRUE)
  }
  
  # Identify hub labels (top hubs per community)
  hub_labels <- rep(NA_character_, nrow(glay))
  if (label_mode == "hubs") {
    s <- if (!is.null(glay$strength)) glay$strength else glay$degree
    tmp <- data.frame(
      idx = seq_len(nrow(glay)),
      name = disp_names,
      community = glay$community,
      score = s
    )
    lab_idx <- unlist(
      by(tmp, tmp$community, function(d) {
        d$idx[order(-d$score)][seq_len(min(hubs_per_comm, nrow(d)))]
      }),
      use.names = FALSE
    )
    if (length(lab_idx)) hub_labels[lab_idx] <- disp_names[lab_idx]
  }
  
  # Define cluster color scale
  node_fill_scale <- {
    levs <- sort(unique(glay$community))
    if (!is.null(cluster_palette)) {
      if (is.null(names(cluster_palette))) {
        vals <- rep(cluster_palette, length.out = length(levs))
        ggplot2::scale_fill_manual(values = setNames(vals, levs), guide = "none")
      } else {
        ggplot2::scale_fill_manual(values = cluster_palette, guide = "none")
      }
    } else {
      ggplot2::scale_fill_brewer(palette = "Set2", guide = "none")
    }
  }
  
  # Build plot
  p <- ggraph(glay) +
    geom_edge_link(
      aes(colour = sign, width = weight, alpha = weight),
      show.legend = FALSE
    ) +
    ggraph::scale_edge_colour_manual(
      values = c(negative = "steelblue", positive = "firebrick"),
      guide = "none"
    ) +
    ggraph::scale_edge_width_continuous(range = c(0.3, 2.8), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
    geom_node_point(
      aes(fill = factor(community), size = ._size),
      shape = 21,
      stroke = 0.3,
      colour = "grey25",
      show.legend = FALSE
    ) +
    {
      if (label_mode == "proteins") {
        geom_node_text(
          aes(label = disp_names),
          repel = TRUE,
          max.overlaps = max_overlaps,
          size = 3
        )
      } else if (label_mode == "hubs") {
        geom_node_text(
          aes(label = hub_labels),
          repel = TRUE,
          max.overlaps = max_overlaps,
          na.rm = TRUE,
          size = 4
        )
      } else {
        NULL
      }
    } +
    ggplot2::scale_size_identity(guide = "none") +
    node_fill_scale +
    theme_graph() +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::guides(
      fill = "none",
      colour = "none",
      size = "none",
      edge_colour = "none",
      edge_alpha = "none",
      edge_width = "none"
    )
  
  invisible(list(tbl_graph = tg, layout_df = glay, plot = p))
}


# ==============================================================================
#  LIGHTWEIGHT WRAPPER: RUN COMPLETE SIGNATURE NETWORK ANALYSIS
# ==============================================================================

#' Run Complete Protein Network Analysis Pipeline
#'
#' Orchestrates steps 1-5: GLASSO fit, graph building, clustering, scoring, viz
#'
#' @param X Numeric matrix of protein data
#' @param lambda GLASSO regularization parameter
#' @param rho_cut Threshold for |ρ| to include edges (default: 0.10)
#' @param ids Optional vector of sample IDs
#' @param layout Layout algorithm (default: "stress")
#' @param node_size Node sizing method (default: "strength")
#' @param label_mode Label display mode (default: "proteins")
#' @param hubs_per_comm Number of hub labels per community (default: 3)
#' @param cluster_palette Manual cluster colors
#' @param walktrap_steps Number of random walk steps (default: 4)
#' @param seed Random seed (default: 1)
#' @return List containing omega, pcor, graph, membership_tbl, scores, plot, layout_df
run_signature_network <- function(
    X,
    lambda,
    rho_cut = 0.10,
    ids = NULL,
    layout = c("stress", "drl", "fr", "kk", "lgl", "circle"),
    node_size = c("strength", "degree", "constant"),
    label_mode = c("proteins", "hubs", "none"),
    hubs_per_comm = 3,
    cluster_palette = NULL,
    walktrap_steps = 4L,
    seed = 1
) {
  
  # Match arguments
  layout <- match.arg(layout)
  node_size <- match.arg(node_size)
  label_mode <- match.arg(label_mode)
  
  # Step 1: Fit GLASSO
  fit <- fit_glasso_pcor(X, lambda)
  
  # Step 2: Build graph
  g <- graph_from_pcor(fit$pcor, rho_cut = rho_cut, node_names = colnames(X))
  
  # Step 3: Cluster
  memb <- walktrap_clusters(g, steps = walktrap_steps)
  
  # Create membership table
  membership_tbl <- tibble(
    protein = V(g)$name,
    cluster = as.integer(memb[V(g)$name])
  ) %>%
    arrange(cluster, protein)
  
  # Step 4: Compute scores
  ids_vec <- ids %||% rownames(X) %||% seq_len(nrow(X))
  scores <- cluster_pc1_scores(X, memb, ids = ids_vec, prefix = "C")
  
  # Step 5: Visualize
  viz <- viz_protein_network(
    g,
    membership = memb,
    edge_cut = rho_cut,
    layout = layout,
    node_size = node_size,
    label_mode = label_mode,
    hubs_per_comm = hubs_per_comm,
    cluster_palette = cluster_palette,
    seed = seed
  )
  
  list(
    omega = fit$omega,
    pcor = fit$pcor,
    graph = g,
    membership_tbl = membership_tbl,
    scores = scores,
    plot = viz$plot,
    layout_df = viz$layout_df
  )
}


# ==============================================================================
#  HELPER: SUBSET GRAPH BY SELECTED CLUSTERS
# ==============================================================================

#' Subset Graph to Selected Clusters
#'
#' @param g igraph object
#' @param membership Named integer vector of cluster assignments
#' @param clusters Integer vector of cluster IDs to keep
#' @return igraph subgraph containing only specified clusters
subset_graph_by_clusters <- function(g, membership, clusters) {
  stopifnot(inherits(g, "igraph"))
  if (vcount(g) == 0) return(g)
  
  # Normalize inputs
  clusters <- as.integer(clusters)
  memb <- as.integer(membership)
  names(memb) <- names(membership)
  
  # Ensure vertices have names
  if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
  
  # Map membership onto vertices, keep only requested clusters
  v_cl <- memb[match(V(g)$name, names(memb))]
  keep_vids <- which(v_cl %in% clusters)
  
  if (!length(keep_vids)) {
    warning("No vertices belong to the selected clusters: ", paste(clusters, collapse = ", "))
    return(make_empty_graph())
  }
  
  induced_subgraph(g, vids = keep_vids)
}


# ==============================================================================
#  HELPER: PLOT SELECTED CLUSTERS ONLY
# ==============================================================================

#' Plot Only Selected Clusters from Network
#'
#' @param g igraph object
#' @param membership Named integer vector of cluster assignments
#' @param clusters Integer vector of cluster IDs to plot
#' @param edge_cut Minimum |ρ| for plotting edges (default: 0.10)
#' @param layout Layout algorithm (default: "fr")
#' @param node_size Node sizing method (default: "strength")
#' @param label_mode Label display mode (default: "proteins")
#' @param hubs_per_comm Number of hub labels per community (default: 3)
#' @param cluster_palette Manual cluster colors
#' @param seed Random seed (default: 1)
#' @return ggplot object or NULL if subgraph is empty
plot_selected_clusters <- function(
    g,
    membership,
    clusters,
    edge_cut = 0.10,
    layout = "fr",
    node_size = "strength",
    label_mode = "proteins",
    hubs_per_comm = 3,
    cluster_palette = NULL,
    seed = 1
) {
  
  # Subset graph
  g_sub <- subset_graph_by_clusters(g, membership, clusters)
  
  # Subset membership to nodes present in subgraph
  memb_sub <- membership[names(membership) %in% V(g_sub)$name]
  
  # Attempt to visualize (may fail if pruning creates empty graph)
  out <- try(
    viz_protein_network(
      g_sub,
      membership = memb_sub,
      edge_cut = edge_cut,
      layout = layout,
      node_size = node_size,
      label_mode = label_mode,
      hubs_per_comm = hubs_per_comm,
      cluster_palette = cluster_palette,
      seed = seed
    ),
    silent = TRUE
  )
  
  if (inherits(out, "try-error")) {
    message("Selected clusters produced an empty graph at edge_cut = ", edge_cut)
    return(NULL)
  }
  
  out$plot
}


# ==============================================================================
#  HELPER: ALIGN CLUSTER SCORES TO ANCHOR LIPOPROTEIN
# ==============================================================================

#' Flip (Sign-Align) Cluster PC1 Scores to Anchor Variable
#'
#' Ensures MPS are positively associated with anchor lipoproteins by
#' multiplying negatively correlated scores by -1.
#'
#' @param scores_df Data frame with cluster score columns (row-aligned to dataset)
#' @param dataset Modeling dataset containing anchor column
#' @param anchor_col Name of anchor column (e.g., "conven_trlp")
#' @param pattern Regex to detect score columns (default: "C\\d+_PC1")
#' @param method Correlation method (default: "pearson")
#' @param verbose Print flipped columns (default: TRUE)
#' @return Revised scores_df with sign-aligned scores
align_scores_to_anchor <- function(
    scores_df,
    dataset,
    anchor_col,
    pattern = "C\\d+_PC1",
    method = "pearson",
    verbose = TRUE
) {
  
  # Validate inputs
  if (!is.data.frame(scores_df)) stop("scores_df must be a data frame.")
  if (!is.data.frame(dataset)) stop("dataset must be a data frame.")
  if (!anchor_col %in% names(dataset)) {
    stop("Anchor column '", anchor_col, "' not found in dataset.")
  }
  if (nrow(scores_df) != nrow(dataset)) {
    stop(
      "Row count mismatch between scores_df (", nrow(scores_df),
      ") and dataset (", nrow(dataset), "). They must be row-aligned."
    )
  }
  
  # Detect cluster PC1 columns
  score_cols <- names(scores_df)[stringr::str_detect(names(scores_df), pattern)]
  if (!length(score_cols)) {
    if (verbose) {
      message("No score columns matching pattern '", pattern, "'. Returning input unchanged.")
    }
    return(scores_df)
  }
  
  # Extract anchor variable
  anchor <- dataset[[anchor_col]]
  out <- scores_df
  flipped <- character(0)
  
  # Check correlation and flip if negative
  for (sc in score_cols) {
    x <- out[[sc]]
    ok <- stats::complete.cases(x, anchor)
    
    if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(anchor[ok]) == 0) {
      if (verbose) {
        message("Skipping '", sc, "': insufficient variability or pairs.")
      }
      next
    }
    
    r <- suppressWarnings(stats::cor(x[ok], anchor[ok], method = method))
    
    if (is.finite(r) && r < 0) {
      out[[sc]] <- -1 * x
      flipped <- c(flipped, sc)
    }
  }
  
  # Report flipped columns
  if (verbose) {
    if (length(flipped)) {
      message(
        "Flipped (multiplied by -1) due to inverse correlation to '",
        anchor_col, "': ", paste(flipped, collapse = ", ")
      )
    } else {
      message("No flips needed for anchor '", anchor_col, "'.")
    }
  }
  
  out
}


# ==============================================================================
#  HELPER: SAVE PLOT WITH TRANSPARENT BACKGROUND
# ==============================================================================

#' Save ggplot/ggraph with Transparent Background
#'
#' @param plot ggplot or ggraph object
#' @param file Output filename (extension determines format: .png, .svg, .pdf)
#' @param width Plot width (default: 10)
#' @param height Plot height (default: 8)
#' @param units Units for width/height (default: "in")
#' @param dpi Resolution for raster formats (default: 300)
#' @param format File format (auto-detected from extension)
save_transparent <- function(
    plot,
    file,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    format = tools::file_ext(file)
) {
  
  # Add transparent background to plot
  plot <- plot +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)
    )
  
  # Choose device based on format
  if (tolower(format) %in% c("png", "")) {
    # PNG with alpha channel (Cairo backend)
    ggplot2::ggsave(
      filename = file,
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      bg = "transparent",
      device = "png",
      type = "cairo"
    )
  } else if (tolower(format) == "svg") {
    # SVG with transparent background
    svglite::svglite(file, width = width, height = height, bg = "transparent")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
  } else if (tolower(format) %in% c("pdf")) {
    # PDF with transparent page (Cairo-based)
    grDevices::cairo_pdf(file, width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(ggplotGrob(plot))
  } else {
    stop("Unsupported format: ", format, " (use .png, .svg, or .pdf)")
  }
}


# ==============================================================================
#  CUSTOM COLOR PALETTE (49 CLUSTERS)
# ==============================================================================

# Hand-curated palette for up to 49 distinct clusters
pal_49 <- c(
  `1`  = "#AEC8A8FF", `2`  = "#5FB2BBFF", `3`  = "#52959CFF", `4`  = "#5C93C2FF",
  `5`  = "#D98281FF", `6`  = "#FF9E9DFF", `7`  = "#BF7549FF", `8`  = "#FF7BA4CC",
  `9`  = "#458A75FF", `10` = "#6CA693FF", `11` = "#7FC7AFFF", `12` = "#456D8BFF",
  `13` = "#F5BBA2FF", `14` = "#8DD5BDFF", `15` = "#4D738FFF", `16` = "#6C8B65FF",
  `17` = "#83B2D9FF", `18` = "#FFC6AEFF", `19` = "#6E9FC5FF", `20` = "#6C8CA5FF",
  `21` = "#DAD8A7FF", `22` = "#3FB8AFFF", `23` = "#75A9D8FF", `24` = "#DC066ECC",
  `25` = "#227E78FF", `26` = "#7FA7C7FF", `27` = "#D20078CC", `28` = "#378276FF",
  `29` = "#92B9D9FF", `30` = "#FF99B6CC", `31` = "#92A78DFF", `32` = "#387A81FF",
  `33` = "#5FBEAEFF", `34` = "#B6B48BFF", `35` = "#73C3CCFF", `36` = "#E2E0AFFF",
  `37` = "#DD4946FF", `38` = "#71CEBEFF", `39` = "#CE9C86FF", `40` = "#409992FF",
  `41` = "#529F92FF", `42` = "#D76789CC", `43` = "#5F85A4FF", `44` = "#3F6586FF",
  `45` = "#547B9FFF", `46` = "#BAD4B4FF", `47` = "#FFB3B2FF", `48` = "#58C9C0FF",
  `49` = "#959351FF"
)


# ==============================================================================
#  EXECUTION: APPLY NETWORK ANALYSIS TO TRL AND LP(A) PROTEIN DATA
# ==============================================================================

# Define lambda regularization parameters and sample IDs
lambda_trl <- 0.06822683
lambda_lpa <- 0.2208841
id_vec <- rownames(zscore_protein)
if (is.null(id_vec)) id_vec <- seq_len(nrow(zscore_protein))

# ------------------------------------------------------------------------------
# STEP 1: Fit GLASSO models (slow step; save results to avoid refitting)
# ------------------------------------------------------------------------------
fit_trl <- fit_glasso_pcor(trl_protein_data, lambda_trl)
fit_lpa <- fit_glasso_pcor(lpa_protein_data, lambda_lpa)

# Optional: Save fitted models to disk
# saveRDS(fit_trl, "fit_trl.rds")
# saveRDS(fit_lpa, "fit_lpa.rds")
# fit_trl <- readRDS("fit_trl.rds")

# ------------------------------------------------------------------------------
# STEP 2: Build graphs from partial correlations
# ------------------------------------------------------------------------------
rho_cut1 <- 0.05  # Threshold for Lp(a) network
rho_cut2 <- 0.10  # Threshold for TRL network

g_trl <- graph_from_pcor(fit_trl$pcor, rho_cut2, node_names = colnames(trl_protein_data))
g_lpa <- graph_from_pcor(fit_lpa$pcor, rho_cut1, node_names = colnames(lpa_protein_data))

# ------------------------------------------------------------------------------
# STEP 3: Cluster proteins using walktrap algorithm
# ------------------------------------------------------------------------------
m_trl <- walktrap_clusters(g_trl, steps = 4)
m_lpa <- walktrap_clusters(g_lpa, steps = 4)

# Inspect membership tables
membership_trl <- tibble(
  protein = V(g_trl)$name,
  cluster = as.integer(m_trl[V(g_trl)$name])
) %>%
  arrange(cluster, protein)

membership_lpa <- tibble(
  protein = V(g_lpa)$name,
  cluster = as.integer(m_lpa[V(g_lpa)$name])
) %>%
  arrange(cluster, protein)

# ------------------------------------------------------------------------------
# STEP 4: Compute multi-protein scores (MPS) as PC1 per cluster
# ------------------------------------------------------------------------------
scores_trl <- cluster_pc1_scores(trl_protein_data, m_trl, ids = id_vec, prefix = "C")
scores_lpa <- cluster_pc1_scores(lpa_protein_data, m_lpa, ids = id_vec, prefix = "C")

# Align scores to anchor lipoproteins
rev_scores_trl <- align_scores_to_anchor(scores_trl, dataset, anchor_col = "conven_trlp")
rev_scores_lpa <- align_scores_to_anchor(scores_lpa, dataset, anchor_col = "conven_lpa_corr")

# Note: LDL alignment requires scores_ldl (not shown in original code)
# rev_scores_ldl <- align_scores_to_anchor(scores_ldl, dataset, anchor_col = "conven_ldlp")

# ------------------------------------------------------------------------------
# STEP 5: Visualize networks
# ------------------------------------------------------------------------------

# Lp(a) network visualization
plt_lpa <- viz_protein_network(
  g_lpa,
  membership = m_lpa,
  edge_cut = rho_cut1,
  layout = "fr",
  node_size = "strength",
  label_mode = "proteins",
  cluster_palette = my_cols_named2,  # Define this palette as needed
  seed = 1,
  hubs_per_comm = 5
)

# TRL network visualization
plt_trl <- viz_protein_network(
  g_trl,
  membership = m_trl,
  edge_cut = rho_cut2,
  layout = "stress",
  node_size = "strength",
  label_mode = "proteins",
  cluster_palette = pal_49_mixed,  # Define this palette as needed
  seed = 1,
  hubs_per_comm = 1
)

# Display plots
print(plt_lpa$plot)
print(plt_trl$plot)

# ------------------------------------------------------------------------------
# STEP 6: Save high-resolution network plots
# ------------------------------------------------------------------------------
save_transparent(
  plt_trl$plot,
  "sfigure-trl-network.png",
  dpi = 1200,
  height = 12,
  width = 10,
  units = "in"
)

save_transparent(
  plt_lpa$plot,
  "sfigure_lpa_network.png",
  dpi = 1200,
  height = 6,
  width = 6,
  units = "in"
)


# ==============================================================================
#  OPTIONAL: CREATE STANDALONE LEGEND FOR 49-CLUSTER PALETTE
# ==============================================================================

# Create dummy data for legend generation
df_legend <- data.frame(cluster = factor(1:49))

p_legend <- ggplot(df_legend, aes(x = cluster, y = 1, color = cluster)) +
  geom_point() +
  scale_color_manual(values = pal_49) +
  guides(
    color = guide_legend(
      nrow = 3,
      byrow = TRUE,
      title = "Cluster"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.7, "lines")
  )

# Extract and save legend only
legend_only <- cowplot::get_legend(p_legend)
cowplot::ggdraw(legend_only)

# Save legend with transparent background
ggplot2::ggsave(
  legend_only,
  file = "sfigure_trl-legend.png",
  dpi = 800,
  height = 1.5,
  width = 9,
  units = "in",
  bg = "transparent",
  device = ragg::agg_png
)
