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
#   Required packages: dplyr, stringr, purrr, tidyr, broom, survival, readr
#
# Goal:
#   Perform Cox proportional hazards regression analyses for ridge-derived
#   MPS - LDL, TRL, and Lp(a) - on incident CAD.
#   Seven progressive models are fitted:
#     1) Clinical covariates only
#     2) Clinical + conventional lipids
#     3) Clinical + other MPS
#     4) Clinical + other MPS + conventional lipids
#     5) Clinical + hsCRP 
#     6) Clinical + WBC 
#     7) Clinical + NLR 
#   
#   Outputs per-exposure hazard ratios and statistics in long format.
#
# Required input data:
#   - dataset.rda
#   - ridge_scores.rda
#
# Main outputs:
#   - results_mps_long.rda
#   - table-ukbb-coxph-mps.csv
#
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(broom)
  library(survival)
  library(readr)
})

# Merge datasets
dataset <- left_join(dataset, ridge_scores, by = "eid")

# ==============================================================================
# DEFINE COVARIATES AND COLUMN ASSUMPTIONS
# ==============================================================================

# Clinical covariates
num_covars <- c("age", "bmi_imp", "conven_hba1c_imp", "sbp_mean_imp")
cat_covars <- c("sex", "race_imp", "smoking_status_imp", "med_lipid")

present_cols <- names(dataset)

# Ridge MPS exposures
mps_candidates  <- c("LDL_ridge_full", "TRL_ridge_full", "TRL_rige_full", "LPA_ridge_full")
mps_present_raw <- intersect(mps_candidates, present_cols)
trl_mps         <- if ("TRL_ridge_full" %in% mps_present_raw) {
  "TRL_ridge_full"
} else {
  intersect("TRL_rige_full", mps_present_raw)
}

mps_present <- c(
  intersect("LDL_ridge_full", mps_present_raw),
  trl_mps,
  intersect("LPA_ridge_full", mps_present_raw)
) |>
  unname() |>
  unique()

if (!length(mps_present)) {
  stop("No ridge MPS exposures found in dataset.")
}

# Conventional lipid covariates (used as adjusters only)
conv_candidates  <- c("conven_ldlp", "conven_trlp", "conven_lpa_corr", "conven_lpa")
conv_present_raw <- intersect(conv_candidates, present_cols)
lpa_cv           <- if ("conven_lpa_corr" %in% conv_present_raw) {
  "conven_lpa_corr"
} else {
  intersect("conven_lpa", conv_present_raw)
}

conv_present <- c(
  intersect("conven_ldlp", conv_present_raw),
  intersect("conven_trlp", conv_present_raw),
  lpa_cv
) |>
  unname() |>
  unique()

# Inflammation / immune covariates
hscrp_present <- intersect("conven_hscrp", present_cols)
wbc_present   <- intersect("blood_wbc", present_cols)
neu_present   <- intersect("blood_neutro", present_cols)
lym_present   <- intersect("blood_lympho", present_cols)
have_nlr      <- length(neu_present) == 1 && length(lym_present) == 1

# ==============================================================================
# BUILD MODELING DATASET
# ==============================================================================

needed <- unique(c(
  "cad2_ukb_survtime", "cad2_ukb_ind",
  mps_present, num_covars, cat_covars,
  conv_present, hscrp_present, wbc_present, neu_present, lym_present
))

df_model <- dataset %>%
  select(any_of(needed)) %>%
  mutate(
    across(all_of(num_covars), as.numeric),
    across(all_of(mps_present), ~ suppressWarnings(as.numeric(.x))),
    across(
      any_of(c(conv_present, hscrp_present, wbc_present, neu_present, lym_present)),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# Derive neutrophil-to-lymphocyte ratio (NLR) if possible
if (have_nlr) {
  df_model <- df_model %>%
    mutate(
      nlr = ifelse(
        !is.na(blood_neutro) & !is.na(blood_lympho) & blood_lympho != 0,
        blood_neutro / blood_lympho,
        NA_real_
      )
    )
}

# ==============================================================================
# DEFINE MODEL SPECIFICATIONS
# ==============================================================================

build_specs_for_mps <- function(exposure_col) {
  other_mps <- setdiff(mps_present, exposure_col)
  
  specs <- list(
    list(name = "Model 1: Clinical covariates", add = character(0)),
    list(name = "Model 2: Clinical cov. + conventional", add = conv_present),
    list(name = "Model 3: Clinical cov. + other MPS", add = other_mps),
    list(name = "Model 4: Clinical cov. + other MPS + conventional", 
         add = c(other_mps, conv_present))
  )
  
  if (length(hscrp_present)) {
    specs <- append(specs, list(
      list(name = "Model 5: Clinical cov. + hsCRP", add = hscrp_present)
    ))
  }
  
  if (length(wbc_present)) {
    specs <- append(specs, list(
      list(name = "Model 6: Clinical cov. + WBC", add = wbc_present)
    ))
  }
  
  if ("nlr" %in% names(df_model)) {
    specs <- append(specs, list(
      list(name = "Model 7: Clinical cov. + NLR", add = "nlr")
    ))
  }
  
  return(specs)
}

# ==============================================================================
# FIT COX MODELS FOR EACH EXPOSURE
# ==============================================================================

fit_exposure_models <- function(exposure_col, df, num_covars, cat_covars) {
  specs <- build_specs_for_mps(exposure_col)
  rhs_base <- c(sprintf("scale(`%s`)", exposure_col), num_covars, cat_covars)
  
  map_dfr(specs, function(sp) {
    add_vars <- intersect(sp$add, names(df))
    rhs <- unique(c(rhs_base, add_vars))
    fml <- as.formula(
      paste0("Surv(cad2_ukb_survtime, cad2_ukb_ind) ~ ", paste(rhs, collapse = " + "))
    )
    
    mdl <- survival::coxph(fml, data = df, ties = "efron")
    
    term_regex <- paste0(
      "^scale\\((`)?",
      str_replace_all(exposure_col, "([\\W])", "\\\\\\1"),
      "(\\1)?\\)$"
    )
    
    td <- broom::tidy(mdl, exponentiate = FALSE, conf.int = TRUE) %>%
      filter(str_detect(term, term_regex)) %>%
      transmute(
        exposure   = exposure_col,
        model_name = sp$name,
        logHR      = estimate,
        SE         = std.error,
        z          = statistic,
        p          = p.value,
        logHR_LCL  = conf.low,
        logHR_UCL  = conf.high
      )
    
    if (!nrow(td)) {
      td <- tibble(
        exposure   = exposure_col,
        model_name = sp$name,
        logHR      = NA_real_,
        SE         = NA_real_,
        z          = NA_real_,
        p          = NA_real_,
        logHR_LCL  = NA_real_,
        logHR_UCL  = NA_real_
      )
    }
    
    mutate(td, HR = exp(logHR), HR_LCL = exp(logHR_LCL), HR_UCL = exp(logHR_UCL))
  })
}

# ==============================================================================
# RUN ANALYSES FOR ALL MPS EXPOSURES
# ==============================================================================

# Define exposure labels
exposure_labels <- c(
  LDL_ridge_full = "LDL MPS",
  TRL_ridge_full = "TRL MPS",
  LPA_ridge_full = "Lp(a) MPS"
)

# Fit models for all MPS exposures
results_mps_long <- map_dfr(
  mps_present,
  ~ fit_exposure_models(.x, df_model, num_covars, cat_covars)
) %>%
  mutate(
    Exposure      = .data$exposure,
    ExposureLabel = unname(exposure_labels[Exposure])
  ) %>%
  select(Exposure, ExposureLabel, model_name, HR, HR_LCL, HR_UCL, logHR, SE, z, p) %>%
  arrange(
    ExposureLabel,
    factor(
      model_name,
      levels = c(
        "Model 1: Clinical covariates",
        "Model 2: Clinical cov. + conventional",
        "Model 3: Clinical cov. + other MPS",
        "Model 4: Clinical cov. + other MPS + conventional",
        "Model 5: Clinical cov. + hsCRP",
        "Model 6: Clinical cov. + WBC",
        "Model 7: Clinical cov. + NLR"
      )
    )
  )

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save long-format results (RDS)
saveRDS(results_mps_long, "results_mps_long.rda")

# Prepare output table for CSV
results_mps_long_output <- results_mps_long %>%
  select(ExposureLabel, model_name, logHR, SE, p) %>%
  rename(
    Variable     = ExposureLabel,
    Model        = model_name,
    `logHR UKBB` = logHR,
    `SE UKBB`    = SE,
    `P UKBB`     = p
  )

# Save CSV output

write_csv(results_mps_long_output, file = "table-ukbb-coxph-mps.csv")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================