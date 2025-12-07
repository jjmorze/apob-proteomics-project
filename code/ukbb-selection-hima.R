# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-08
# Last modified: 2025-12-07
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: HIMA, survival, dplyr, tibble, ggplot2
#
# Goal:
#   Perform HIMA (High-dimensional Mediation Analysis) survival mediation
#   analysis for lipoprotein fractions (TRL and Lp(a)) to identify protein
#   cluster mediators of the exposure-CAD relationship. Analyzes clusters
#   with n > 1 proteins using PC1 scores as mediators in a survival framework.
#
# Required input objects:
#   - df_model_trl           # Modeling dataset for TRL fraction
#   - df_model_lpa           # Modeling dataset for Lp(a) fraction
#   - membership_trl         # Cluster membership for TRL proteins
#   - membership_lpa         # Cluster membership for Lp(a) proteins
#
# Main outputs:
#   - hima_results           # List containing full HIMA results per fraction
#   - hima_summary_long      # Combined summary table (all fractions)
#   - hima_selected_all      # Significant mediators (FDR < 0.05)
#
# Optional outputs (commented):
#   - hima_by_fraction_results.rda
#
# ==============================================================================


# ==============================================================================
# DEPENDENCIES
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("HIMA", quietly = TRUE)) install.packages("HIMA")
  library(HIMA)
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# Set seed for reproducibility
set.seed(20251003)


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Prepare cluster metadata and filter clusters with > 1 protein
#'
#' @param membership_tbl Data frame with columns: protein, cluster
#' @return Tibble with cluster, n_prot, and proteins (comma-separated)
prep_cluster_meta <- function(membership_tbl) {
  stopifnot(
    is.data.frame(membership_tbl),
    all(c("protein", "cluster") %in% names(membership_tbl))
  )
  
  membership_tbl %>%
    mutate(protein_clean = sub("_rz$", "", protein)) %>%
    group_by(cluster) %>%
    summarise(
      n_prot   = n(),
      proteins = paste(sort(protein_clean), collapse = ", "),
      .groups  = "drop"
    ) %>%
    filter(n_prot > 1) %>%
    arrange(cluster)
}


#' Robust HIMA survival function caller
#'
#' Supports both hima_survival (newer) and hima.survival (older) naming
.hima_surv_fn <- local({
  f <- get0("hima_survival", envir = asNamespace("HIMA"), mode = "function")
  if (is.null(f)) {
    f <- get0("hima.survival", envir = asNamespace("HIMA"), mode = "function")
  }
  if (is.null(f)) {
    stop(
      "Couldn't find HIMA::hima_survival or HIMA::hima.survival. ",
      "Please update the HIMA package."
    )
  }
  f
})


#' Validate presence of required columns in data frame
#'
#' @param df Data frame to validate
#' @param cols Character vector of required column names
#' @param ctx Context description for error message
assert_has_cols <- function(df, cols, ctx = "data") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) {
    stop(
      sprintf("Missing columns in %s: %s", ctx, paste(miss, collapse = ", "))
    )
  }
}


# ==============================================================================
# CORE HIMA RUNNER
# ==============================================================================

#' Run HIMA survival mediation analysis for a single lipoprotein fraction
#'
#' @param frac_label Fraction name (e.g., "TRL", "Lp(a)")
#' @param df_fraction Modeling dataset containing outcome, exposure, mediators
#' @param membership_tbl Cluster membership table (protein, cluster)
#' @param exposure_var Name of exposure variable column
#' @param covariates Character vector of covariate column names
#' @param topN Maximum number of top mediators to screen (default: 5)
#' @param FDRcut FDR cutoff for significance (default: 0.05)
#' @param verbose Logical; print progress messages (default: TRUE)
#'
#' @return List containing:
#'   - Fraction: fraction label
#'   - exposure: exposure variable name
#'   - mediators: mediator column names used
#'   - n_mediators: count of mediators
#'   - hima: raw HIMA fit object
#'   - summary: full summary table
#'   - selected: significant mediators (FDR < FDRcut)
#'   - plot: ggplot object showing top mediators
run_hima_for_fraction <- function(
    frac_label,
    df_fraction,
    membership_tbl,
    exposure_var,
    covariates = c(
      "age", "sex", "bmi_imp", "race_imp",
      "smoking_status_imp", "sbp_mean_imp",
      "conven_hba1c_imp", "med_lipid"
    ),
    topN = 5,
    FDRcut = 0.05,
    verbose = TRUE
) {
  
  # --------------------------------------------------------------------------
  # 1. Identify eligible mediators (clusters with > 1 protein)
  # --------------------------------------------------------------------------
  meta <- prep_cluster_meta(membership_tbl)
  mediator_cols <- paste0("C", meta$cluster, "_PC1")
  mediator_cols <- intersect(mediator_cols, colnames(df_fraction))
  
  if (length(mediator_cols) == 0L) {
    if (verbose) {
      message(
        sprintf("[%s] No eligible mediator columns (n>1) found.", frac_label)
      )
    }
    return(list(
      Fraction    = frac_label,
      exposure    = exposure_var,
      mediators   = character(0),
      n_mediators = 0L,
      hima        = NULL,
      summary     = tibble(),
      selected    = tibble(),
      plot        = NULL
    ))
  }
  
  # --------------------------------------------------------------------------
  # 2. Build analysis frame and validate required columns
  # --------------------------------------------------------------------------
  needed <- c(
    "cad2_ukb_survtime", "cad2_ukb_ind",
    exposure_var, covariates, mediator_cols
  )
  assert_has_cols(
    df_fraction, needed,
    ctx = sprintf("df_fraction (%s)", frac_label)
  )
  
  df_use <- df_fraction %>% dplyr::select(all_of(needed))
  
  # --------------------------------------------------------------------------
  # 3. Construct HIMA input components
  # --------------------------------------------------------------------------
  X_vec      <- df_use[[exposure_var]]
  M_mat      <- as.matrix(df_use[, mediator_cols, drop = FALSE])
  OT_vec     <- df_use$cad2_ukb_survtime
  status_vec <- df_use$cad2_ukb_ind
  
  # Build covariate model matrix (without intercept)
  covariate_formula <- reformulate(covariates)
  COV_mat <- model.matrix(covariate_formula, data = df_use)[, -1, drop = FALSE]
  
  # --------------------------------------------------------------------------
  # 4. Filter to complete cases (HIMA requires no NAs)
  # --------------------------------------------------------------------------
  cc <- complete.cases(X_vec, OT_vec, status_vec, COV_mat, M_mat)
  if (!any(cc)) {
    stop(
      sprintf(
        "[%s] No complete cases for X, M, OT, status, and covariates.",
        frac_label
      )
    )
  }
  
  X_vec_cc  <- X_vec[cc]
  M_mat_cc  <- M_mat[cc, , drop = FALSE]
  OT_vec_cc <- OT_vec[cc]
  status_cc <- status_vec[cc]
  COV_cc    <- COV_mat[cc, , drop = FALSE]
  
  # --------------------------------------------------------------------------
  # 5. Limit topN to available mediators
  # --------------------------------------------------------------------------
  topN_eff <- min(topN, ncol(M_mat_cc))
  
  # --------------------------------------------------------------------------
  # 6. Run HIMA survival mediation analysis
  # --------------------------------------------------------------------------
  if (verbose) {
    message(
      sprintf(
        "[%s] Running HIMA with %d mediators and %d complete cases.",
        frac_label, ncol(M_mat_cc), sum(cc)
      )
    )
  }
  
  hima_fit <- .hima_surv_fn(
    X       = X_vec_cc,
    M       = M_mat_cc,
    OT      = OT_vec_cc,
    status  = status_cc,
    COV     = COV_cc,
    topN    = topN_eff,
    scale   = TRUE,
    FDRcut  = FDRcut,
    verbose = verbose
  )
  
  # --------------------------------------------------------------------------
  # 7. Generate summaries and visualization
  # --------------------------------------------------------------------------
  summ <- tibble()
  sel  <- tibble()
  plt  <- NULL
  
  if (!is.null(hima_fit$summary) && NROW(hima_fit$summary) > 0) {
    summ <- as_tibble(hima_fit$summary) %>%
      mutate(Fraction = frac_label, .before = 1)
    
    sel <- summ %>%
      filter(FDR < FDRcut) %>%
      arrange(FDR)
    
    # Create barplot for top mediators (up to 10)
    topn <- min(10L, nrow(summ))
    plot_df <- summ %>%
      arrange(FDR) %>%
      slice_head(n = topn) %>%
      mutate(neglog10FDR = -log10(pmax(FDR, .Machine$double.eps)))
    
    if (nrow(plot_df) > 0) {
      plt <- ggplot(
        plot_df,
        aes(x = reorder(mediator, neglog10FDR), y = neglog10FDR)
      ) +
        geom_col() +
        coord_flip() +
        labs(
          x = "Mediator (cluster PC1)",
          y = expression(-log[10](FDR)),
          title = sprintf("Top %d mediators (%s)", topn, frac_label),
          subtitle = sprintf("Exposure: %s; FDRcut = %.2f", exposure_var, FDRcut)
        ) +
        theme_minimal()
    }
  } else if (verbose) {
    message(sprintf("[%s] HIMA returned no rows in $summary.", frac_label))
  }
  
  # --------------------------------------------------------------------------
  # 8. Return results bundle
  # --------------------------------------------------------------------------
  list(
    Fraction    = frac_label,
    exposure    = exposure_var,
    mediators   = mediator_cols,
    n_mediators = length(mediator_cols),
    hima        = hima_fit,
    summary     = summ,
    selected    = sel,
    plot        = plt
  )
}


# ==============================================================================
# ANALYSIS EXECUTION
# ==============================================================================

# ------------------------------------------------------------------------------
# Run HIMA for TRL fraction
# ------------------------------------------------------------------------------
hima_TRL <- run_hima_for_fraction(
  frac_label     = "TRL",
  df_fraction    = df_model_trl,
  membership_tbl = membership_trl,
  exposure_var   = "conven_trlp",
  topN           = 49,
  FDRcut         = 0.05,
  verbose        = TRUE
)

# ------------------------------------------------------------------------------
# Run HIMA for Lp(a) fraction
# ------------------------------------------------------------------------------
hima_LPA <- run_hima_for_fraction(
  frac_label     = "Lp(a)",
  df_fraction    = df_model_lpa,
  membership_tbl = membership_lpa,
  exposure_var   = "conven_lpa_corr",
  topN           = 6,
  FDRcut         = 0.05,
  verbose        = TRUE
)


# ==============================================================================
# COMBINE AND SUMMARIZE RESULTS
# ==============================================================================

# Bundle all results
hima_results <- list(
  TRL   = hima_TRL,
  `Lp(a)` = hima_LPA
)

# Combined summary table (long format)
hima_summary_long <- bind_rows(
  hima_TRL$summary,
  hima_LPA$summary
)

# Extract significant mediators across all fractions
hima_selected_all <- bind_rows(
  hima_TRL$selected %>% mutate(Fraction = "TRL", .before = 1),
  hima_LPA$selected %>% mutate(Fraction = "Lp(a)", .before = 1)
)

# Print summary statistics
cat(
  sprintf(
    "HIMA summaries - TRL: %d mediators, Lp(a): %d mediators\n",
    nrow(hima_TRL$summary),
    nrow(hima_LPA$summary)
  )
)

cat("\nSignificant mediators (FDR < 0.05):\n")
print(hima_selected_all)


# ==============================================================================
# OPTIONAL: VISUALIZE AND SAVE RESULTS
# ==============================================================================

# Display plots (uncomment to view)
# print(hima_TRL$plot)
# print(hima_LPA$plot)

# Save results to file (uncomment to save)
# save(
#   hima_results,
#   hima_summary_long,
#   hima_selected_all,
#   file = "hima_by_fraction_results.rda"
# )
