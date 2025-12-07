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
#   Required packages: dplyr, stringr, tibble, purrr, CMAverse
#
# Goal:
#   Perform causal mediation analysis using CMAverse for three lipid exposures:
#   1) Triglyceride-rich lipoproteins (TRL)
#   2) Lipoprotein(a) (Lp(a))
#   
#
# Required input data:
#   - df_model_trl           # TRL analysis dataset
#   - df_model_lpa           # Lp(a) analysis dataset
#   - membership_trl         # Protein cluster membership for TRL
#   - membership_lpa         # Protein cluster membership for Lp(a)
#
# Main outputs:
#   - fit_trlp               # CMAverse fit object for TRL
#   - fit_lpa                # CMAverse fit object for Lp(a)
#   - trl_effects            # Effect estimates for TRL
#   - lpa_effects            # Effect estimates for Lp(a)
#   - med_data_trl           # Mediation dataset for TRL
#   - med_data_lpa           # Mediation dataset for Lp(a)
#
# ==============================================================================


# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(CMAverse)
})


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

#' Filter clusters with more than one protein
#'
#' @param membership_tbl Data frame with 'protein' and 'cluster' columns
#' @return Vector of cluster IDs that contain >1 protein
clusters_n_gt1 <- function(membership_tbl) {
  stopifnot(
    is.data.frame(membership_tbl),
    all(c("protein", "cluster") %in% names(membership_tbl))
  )
  
  membership_tbl %>%
    count(cluster, name = "n_prot") %>%
    filter(n_prot > 1) %>%
    pull(cluster)
}


#' Safely z-score a numeric vector
#'
#' @param x Numeric vector to standardize
#' @return Z-scored vector (NA if non-numeric or zero variance)
z_safe <- function(x) {
  if (!is.numeric(x)) {
    return(rep(NA_real_, length(x)))
  }
  
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  
  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  (x - m) / s
}


#' Create compact mediation dataset for CMAverse analysis
#'
#' Prepares data for causal mediation analysis by:
#' - Selecting exposure, mediators, outcome, and covariates
#' - Optionally filtering clusters to those with >1 protein
#' - Z-scoring exposure and mediators for standardized effects
#' - Removing missing data
#'
#' @param data Data frame containing all variables
#' @param exposure Name of exposure variable
#' @param outcome_time Name of survival time variable
#' @param event Name of event indicator variable
#' @param cluster_ids Vector of cluster IDs for mediators
#' @param covariates Vector of covariate names
#' @param membership Data frame with protein-cluster membership (optional)
#' @param enforce_n_gt1 Logical; enforce >1 protein per cluster?
#' @param make_z Logical; z-score mediators?
#' @param make_exposure_z Logical; z-score exposure?
#' @param col_pattern Pattern for mediator column names (use %s for cluster ID)
#' @param z_suffix Suffix for z-scored variables
#' @return List with mediation dataset and metadata
create_mediation_dataset <- function(data,
                                     exposure,
                                     outcome_time  = "cad2_ukb_survtime",
                                     event         = "cad2_ukb_ind",
                                     cluster_ids,
                                     covariates,
                                     membership    = NULL,
                                     enforce_n_gt1 = TRUE,
                                     make_z        = TRUE,
                                     make_exposure_z = TRUE,
                                     col_pattern   = "C%s_PC1",
                                     z_suffix      = "_z") {
  
  # Validate inputs
  stopifnot(
    is.data.frame(data),
    is.character(exposure),
    length(exposure) == 1
  )
  
  # Filter cluster IDs to those with >1 protein (if requested)
  if (enforce_n_gt1 && !is.null(membership)) {
    valid_ids <- clusters_n_gt1(membership)
    requested <- cluster_ids
    cluster_ids <- intersect(requested, valid_ids)
    
    if (!length(cluster_ids)) {
      stop(
        "No clusters left after enforcing n>1 proteins. Requested: ",
        paste(requested, collapse = ", ")
      )
    }
    
    dropped <- setdiff(requested, cluster_ids)
    if (length(dropped)) {
      message(
        "Dropped clusters (n<=1 proteins): ",
        paste(dropped, collapse = ", ")
      )
    }
  }
  
  # Identify mediator columns that exist in data
  mediator_cols <- sprintf(col_pattern, cluster_ids)
  mediator_cols <- intersect(mediator_cols, colnames(data))
  
  if (!length(mediator_cols)) {
    stop(
      "None of the requested mediator columns exist in 'data'. ",
      "Pattern: ", col_pattern, " IDs: ",
      paste(cluster_ids, collapse = ", ")
    )
  }
  
  # Check for missing required columns
  needed <- unique(c(exposure, outcome_time, event, covariates, mediator_cols))
  missing <- setdiff(needed, colnames(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Select required columns
  df <- data[, needed, drop = FALSE]
  
  # Z-score mediators
  mediator_z <- if (make_z) paste0(mediator_cols, z_suffix) else mediator_cols
  if (make_z) {
    for (i in seq_along(mediator_cols)) {
      mc <- mediator_cols[i]
      df[[mediator_z[i]]] <- z_safe(df[[mc]])
    }
  }
  
  # Z-score exposure
  exposure_use <- exposure
  if (make_exposure_z) {
    exposure_use <- paste0(exposure, z_suffix)
    df[[exposure_use]] <- z_safe(df[[exposure]])
  }
  
  # Identify variables used by cmest
  used_for_cmest <- c(outcome_time, event, exposure_use, covariates, mediator_z)
  
  # Remove missing data and tidy column order
  df <- df %>%
    dplyr::select(all_of(used_for_cmest)) %>%
    stats::na.omit() %>%
    dplyr::relocate(all_of(c(outcome_time, event, exposure_use))) %>%
    dplyr::relocate(all_of(covariates), .after = exposure_use)
  
  # Return results
  list(
    data             = df,
    exposure_var     = exposure_use,
    mediators_raw    = mediator_cols,
    mediators_z      = mediator_z,
    cluster_ids_kept = cluster_ids
  )
}


#' Run CMAverse mediation analysis
#'
#' Wrapper for CMAverse::cmest with standardized settings:
#' - astar = 0, a = 1 (comparing +1 SD vs. 0 SD on z-scored exposure)
#' - Rare binary outcome model with imputation-based inference
#' - Bootstrap confidence intervals
#'
#' @param data Data frame with complete cases
#' @param outcome Name of survival time variable
#' @param event Name of event indicator variable
#' @param exposure Name of (z-scored) exposure variable
#' @param mediators Vector of (z-scored) mediator variable names
#' @param covariates Vector of covariate names
#' @param nboot Number of bootstrap iterations
#' @return CMAverse fit object
run_cmest <- function(data,
                      outcome,
                      event,
                      exposure,
                      mediators,
                      covariates,
                      nboot = 1000L) {
  
  stopifnot(length(mediators) > 0)
  
  CMAverse::cmest(
    data       = data,
    model      = "rb",                                    # Rare binary outcome
    outcome    = outcome,
    event      = event,
    exposure   = exposure,
    mediator   = mediators,
    basec      = covariates,
    EMint      = FALSE,                                   # No exposure-mediator interaction
    mreg       = rep(list("linear"), length(mediators)),  # Linear models for mediators
    yreg       = "coxph",                                 # Cox proportional hazards
    astar      = 0,                                       # Reference level (0 SD)
    a          = 1,                                       # Comparison level (+1 SD)
    mval       = rep(list(0), length(mediators)),         # Mediator values at reference
    estimation = "imputation",
    inference  = "bootstrap",
    nboot      = nboot
  )
}


#' Extract effect estimates from CMAverse fit
#'
#' @param fit CMAverse fit object
#' @return Tibble with effect, estimate, and confidence intervals
extract_effects <- function(fit) {
  if (is.null(fit$effect.pe) || is.null(fit$effect.ci)) {
    return(NULL)
  }
  
  tibble(
    effect   = names(fit$effect.pe),
    estimate = unname(fit$effect.pe),
    lcl      = fit$effect.ci[, 1],
    ucl      = fit$effect.ci[, 2]
  )
}


# ==============================================================================
#  SHARED PARAMETERS
# ==============================================================================

# Covariates used in all models
covariates <- c(
  "age",
  "sex",
  "bmi_imp",
  "race_imp",
  "smoking_status_imp",
  "sbp_mean_imp",
  "conven_hba1c_imp",
  "med_lipid"
)

# Set random seed for reproducibility
set.seed(20251003)


# ==============================================================================
#  TRIGLYCERIDE-RICH LIPOPROTEINS (TRL) MEDIATION ANALYSIS
# ==============================================================================

message("\n=== TRL Mediation Analysis ===")

# Build mediation dataset
trl_build <- create_mediation_dataset(
  data            = df_model_trl,
  exposure        = "conven_trlp",
  outcome_time    = "cad2_ukb_survtime",
  event           = "cad2_ukb_ind",
  cluster_ids     = c(8, 13, 34),
  covariates      = covariates,
  membership      = membership_trl,
  enforce_n_gt1   = TRUE,
  make_z          = TRUE,
  make_exposure_z = TRUE
)

# Extract components
med_data_trl    <- trl_build$data
trl_exposure    <- trl_build$exposure_var
trl_mediators_z <- trl_build$mediators_z

message("TRL dataset: ", nrow(med_data_trl), " complete cases")
message("Mediators: ", paste(trl_mediators_z, collapse = ", "))

# Run mediation analysis
fit_trlp <- run_cmest(
  data       = med_data_trl,
  outcome    = "cad2_ukb_survtime",
  event      = "cad2_ukb_ind",
  exposure   = trl_exposure,
  mediators  = trl_mediators_z,
  covariates = covariates,
  nboot      = 1000
)

# Extract effects
trl_effects <- extract_effects(fit_trlp)


# ==============================================================================
#  LIPOPROTEIN(A) MEDIATION ANALYSIS
# ==============================================================================

message("\n=== Lp(a) Mediation Analysis ===")

# Build mediation dataset
lpa_build <- create_mediation_dataset(
  data            = df_model_lpa,
  exposure        = "conven_lpa_corr",
  outcome_time    = "cad2_ukb_survtime",
  event           = "cad2_ukb_ind",
  cluster_ids     = c(1, 3, 6),
  covariates      = covariates,
  membership      = membership_lpa,
  enforce_n_gt1   = TRUE,
  make_z          = TRUE,
  make_exposure_z = TRUE
)

# Extract components
med_data_lpa    <- lpa_build$data
lpa_exposure    <- lpa_build$exposure_var
lpa_mediators_z <- lpa_build$mediators_z

message("Lp(a) dataset: ", nrow(med_data_lpa), " complete cases")
message("Mediators: ", paste(lpa_mediators_z, collapse = ", "))

# Run mediation analysis
fit_lpa <- run_cmest(
  data       = med_data_lpa,
  outcome    = "cad2_ukb_survtime",
  event      = "cad2_ukb_ind",
  exposure   = lpa_exposure,
  mediators  = lpa_mediators_z,
  covariates = covariates,
  nboot      = 1000
)

# Extract effects
lpa_effects <- extract_effects(fit_lpa)


# ==============================================================================
#  DISPLAY RESULTS
# ==============================================================================

message("\n=== TRL Results ===")
print(fit_trlp)
if (!is.null(trl_effects)) {
  print(trl_effects)
}


message("\n=== Lp(a) Results ===")
print(fit_lpa)
if (!is.null(lpa_effects)) {
  print(lpa_effects)
}

message("\n=== Analysis Complete ===")
