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
#   Perform Cox proportional hazards regression analysis to assess the 
#   association between conventional lipoproteins (LDL, TRL, Lp(a)) and 
#   incident CAD (cad2_ukb). Two modeling strategies:
#   - M1: Each lipoprotein adjusted for clinical covariates (one-at-a-time)
#   - M2: All lipoproteins mutually adjusted + clinical covariates
#
# Required input data:
#   - dataset.rda 
#
# Main outputs:
#   - table-ukbb-coxph-lipoproteins.csv 
#   - results_conv_output.rda
#
# ==============================================================================

# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(broom)
  library(survival)
  library(readr)
})

# Load  UK Biobank dataset
load("~PATH/dataset.rda")

# ==============================================================================
# DEFINE VARIABLES AND COLUMN ASSUMPTIONS
# ==============================================================================

# Time-to-event outcome: cad2_ukb_survtime (time), cad2_ukb_ind (event)

# Clinical covariates (adjustments)
num_covars <- c("age", "bmi_imp", "conven_hba1c_imp", "sbp_mean_imp")
cat_covars <- c("sex", "race_imp", "smoking_status_imp", "med_lipid")

# Conventional lipoprotein exposures (REQUIRED: conven_lpa_corr)
conv_candidates <- c("conven_ldlp", "conven_trlp", "conven_lpa_corr")

# Labels for display
conv_labels <- c(
  conven_ldlp     = "LDL",
  conven_trlp     = "TRL",
  conven_lpa_corr = "Lp(a)"
)

exposure_order <- c("LDL", "TRL", "Lp(a)")

# ==============================================================================
# VALIDATE AND SELECT EXPOSURE VARIABLES
# ==============================================================================

present_cols <- names(dataset)
conv_present_raw <- intersect(conv_candidates, present_cols)

# Ensure conven_lpa_corr is present
if (!("conven_lpa_corr" %in% conv_present_raw)) {
  stop("Required variable 'conven_lpa_corr' is missing from 'dataset'.")
}

# Build final exposure list (ensure Lp(a) is included)
conv_present <- c(
  intersect("conven_ldlp", conv_present_raw),
  intersect("conven_trlp", conv_present_raw),
  "conven_lpa_corr"
) |> 
  unname() |> 
  unique()

# ==============================================================================
# PREPARE MODELING DATASET
# ==============================================================================

needed <- unique(c(
  "cad2_ukb_survtime", 
  "cad2_ukb_ind", 
  num_covars, 
  cat_covars, 
  conv_present
))

df_model <- dataset %>%
  select(any_of(needed)) %>%
  mutate(
    across(all_of(num_covars), ~ suppressWarnings(as.numeric(.))),
    across(all_of(conv_present), ~ suppressWarnings(as.numeric(.))),
    across(all_of(cat_covars), ~ as.factor(.))
  )

# ==============================================================================
# HELPER FUNCTION: EXTRACT SCALED EXPOSURE TERMS
# ==============================================================================

tidy_keep_scaled <- function(mdl, keep_vars) {
  # Build regex to match scale(variable) terms
  parts <- paste0(
    "scale\\((`)?", 
    str_replace_all(keep_vars, "([\\W])", "\\\\\\1"), 
    "(\\1)?\\)"
  )
  rx <- paste0("^(", paste(parts, collapse = "|"), ")$")
  
  broom::tidy(mdl, exponentiate = FALSE, conf.int = TRUE) %>%
    filter(str_detect(term, rx))
}

# ==============================================================================
# MODEL 1: CLINICAL COVARIATES (ONE-AT-A-TIME)
# ==============================================================================

fit_M1 <- function(exposure_col, df) {
  # Build formula: scale(exposure) + clinical covariates
  rhs <- c(sprintf("scale(`%s`)", exposure_col), num_covars, cat_covars)
  fml <- as.formula(paste0(
    "Surv(cad2_ukb_survtime, cad2_ukb_ind) ~ ", 
    paste(rhs, collapse = " + ")
  ))
  
  # Fit Cox model
  mdl <- survival::coxph(fml, data = df, ties = "efron")
  
  # Extract and format results
  tidy_keep_scaled(mdl, exposure_col) %>%
    transmute(
      exposure   = exposure_col,
      model_name = "M1: Clinical",
      logHR      = estimate,
      SE         = std.error,
      z          = statistic,
      p          = p.value,
      logHR_LCL  = conf.low,
      logHR_UCL  = conf.high
    ) %>%
    mutate(
      HR     = exp(logHR), 
      HR_LCL = exp(logHR_LCL), 
      HR_UCL = exp(logHR_UCL)
    )
}

# Fit M1 for each lipoprotein
results_M1 <- purrr::map_dfr(conv_present, fit_M1, df = df_model)

# ==============================================================================
# MODEL 2: CLINICAL + MUTUAL ADJUSTMENT FOR LIPOPROTEINS
# ==============================================================================

fit_M2 <- function(exposure_cols, df) {
  # Build formula: scale(all exposures) + clinical covariates
  rhs <- c(sprintf("scale(`%s`)", exposure_cols), num_covars, cat_covars)
  fml <- as.formula(paste0(
    "Surv(cad2_ukb_survtime, cad2_ukb_ind) ~ ", 
    paste(rhs, collapse = " + ")
  ))
  
  # Fit Cox model
  mdl <- survival::coxph(fml, data = df, ties = "efron")
  
  # Extract and format results
  tidy_keep_scaled(mdl, exposure_cols) %>%
    transmute(
      exposure   = str_replace_all(term, "^scale\\((`)?|(`)?\\)$", ""),
      model_name = "M2: Clinical + mutual adjustment for lipoproteins",
      logHR      = estimate,
      SE         = std.error,
      z          = statistic,
      p          = p.value,
      logHR_LCL  = conf.low,
      logHR_UCL  = conf.high
    ) %>%
    mutate(
      HR     = exp(logHR), 
      HR_LCL = exp(logHR_LCL), 
      HR_UCL = exp(logHR_UCL)
    )
}

# Fit M2 with all lipoproteins
results_M2 <- fit_M2(conv_present, df_model)

# ==============================================================================
# COMBINE AND FORMAT RESULTS
# ==============================================================================

results_conv_long <- bind_rows(results_M1, results_M2) %>%
  mutate(
    Exposure      = exposure,
    ExposureLabel = unname(conv_labels[Exposure]),
    model_name    = factor(
      model_name,
      levels = c(
        "M1: Clinical",
        "M2: Clinical + mutual adjustment for lipoproteins"
      )
    ),
    ExposureLabel = factor(ExposureLabel, levels = exposure_order)
  ) %>%
  arrange(ExposureLabel, model_name) %>%
  select(Exposure, ExposureLabel, model_name, HR, HR_LCL, HR_UCL, 
         logHR, SE, z, p)

# Save raw results
save(results_conv_output, file = "results_conv_output.rda")

# ==============================================================================
# CREATE COMPACT OUTPUT TABLE
# ==============================================================================

results_conv_output <- results_conv_long %>%
  select(ExposureLabel, model_name, logHR, SE, p) %>%
  rename(
    Variable     = ExposureLabel,
    Model        = model_name,
    `logHR UKBB` = logHR,
    `SE UKBB`    = SE,
    `P UKBB`     = p
  )

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write_csv(results_conv_output, file = "PATH/table-ukbb-coxph-lipoproteins.csv")