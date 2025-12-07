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
#   Required packages: tidyverse, gtsummary
#
# Goal:
#   Generate baseline characteristics table for UK Biobank cohort including:
#   1) Converting HbA1c from IFCC (mmol/mol) to NGSP (%) units
#   2) Creating a summary table of demographic, clinical, and laboratory variables
#   3) Formatting variables with appropriate labels and precision
#
# Required input data:
#   - dataset.rda (final merged analysis dataset from main pipeline)
#
# Outputs:
#   - Formatted baseline characteristics table (gtsummary object)
#
# Notes:
#   - HbA1c conversion formula: NGSP% = (IFCC mmol/mol / 10.929) + 2.15
#   - Table includes demographics, anthropometrics, metabolic markers,
#     lipid panel, and medication use
#
# ==============================================================================

# Load required libraries
library(tidyverse)
library(gtsummary)

# Load dataset
load("dataset.rda")

# ==============================================================================
#  DATA PREPARATION
# ==============================================================================

# Convert HbA1c from IFCC (mmol/mol) to NGSP (%)
dataset <- dataset %>%
  mutate(conven_hba1c_percent = conven_hba1c / 10.929 + 2.15)

# ==============================================================================
#  BASELINE CHARACTERISTICS TABLE
# ==============================================================================

# Generate summary table
baseline_table <- dataset %>%
  select(
    age,
    sex,
    bmi,
    race,
    smoking_status,
    conven_hba1c_percent,
    sbp_mean,
    med_statin,
    med_lipid,
    conven_apob,
    conven_ldl,
    conven_rchol,
    conven_tg,
    conven_lpa_corr,
    conven_hdl
  ) %>%
  tbl_summary(
    label = list(
      age                  ~ "Age, years",
      sex                  ~ "Female sex",
      race                 ~ "Race/ethnicity",
      bmi                  ~ "Body mass index, kg/m²",
      smoking_status       ~ "Current smokers",
      conven_hba1c_percent ~ "Glycated hemoglobin, %",
      sbp_mean             ~ "Systolic blood pressure, mmHg",
      med_statin           ~ "Use of statins",
      med_lipid            ~ "Use of any lipid-lowering medication",
      conven_apob          ~ "Apolipoprotein B, g/L",
      conven_ldl           ~ "LDL cholesterol, mmol/L",
      conven_rchol         ~ "Remnant cholesterol, mmol/L",
      conven_tg            ~ "Triglycerides, mmol/L",
      conven_lpa_corr      ~ "Lipoprotein(a), nmol/L",
      conven_hdl           ~ "HDL cholesterol, mmol/L"
    ),
    missing = "no",
    digits = list(
      conven_hba1c_percent ~ 1,
      conven_rchol         ~ 2,
      conven_ldl           ~ 1,
      conven_hdl           ~ 1
    ),
    value = list(
      sex            ~ "Female",
      smoking_status ~ "Current"
    )
  )

# Display table
baseline_table
