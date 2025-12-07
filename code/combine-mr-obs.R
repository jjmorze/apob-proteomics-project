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
#   Required packages: readxl, dplyr, tidyr, stringr, readr
#
# Goal:
#   Process and validate Mendelian randomization (MR) results by:
#   1) Reading and cleaning MR results from Excel
#   2) Rescaling estimates to 1-SD of lipid values in UKBB (z-score scale)
#   3) Merging with observational UKBB protein associations
#   4) Assessing agreement between observational and genetic estimates
#   5) Reshaping data to long format for downstream analysis
#   6) Exporting validated results
#
# Required input data:
#   - mr_results.xlsx
#   - model_1.rda
#
# Main outputs:
#   - ukbb_mr_valid_output.csv 
#   - combined.rda
#   - combined_long.rda
#
# ==============================================================================

# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ==============================================================================
# READ MR RESULTS
# ==============================================================================

mr_results <- read_excel(
  "~/mr_results.xlsx"
)

# ==============================================================================
# CLEAN AND STANDARDIZE MR RESULTS
# ==============================================================================

# Rescale estimates to 1-SD of lipid values in UKBB dataset (N = 35,475)
# Original MR estimates:
#   - LDL: per 1 mmol/L
#   - TRL: per 1 mmol/L
#   - Lp(a): per 1 nmol/L
# UKBB standard deviations:
#   - TRL: 75.33189 mg/dL
#   - LDL: 432.5011 mg/dL
#   - Lp(a): 74.66814 nmol/L

mr_results_clean <- mr_results %>%
  rename(
    # Protein identifiers
    protein_name = protein_name,
    uniprot      = uniprot,
    gene         = gene,
    
    # TRL-related variables
    n.snp.trl       = NSNPs_TRL,
    estimate.mr.trl = MR_beta_TRL,
    ci.lb.mr.trl    = MR_CILower_TRL,
    ci.ub.mr.trl    = MR_CIUpper_TRL,
    p.value.mr.trl  = MR_P_Value_TRL,
    
    # LDL-related variables
    n.snp.ldl       = NSNPs_LDL,
    estimate.mr.ldl = MR_beta_LDL,
    ci.lb.mr.ldl    = MR_CILower_LDL,
    ci.ub.mr.ldl    = MR_CIUpper_LDL,
    p.value.mr.ldl  = MR_P_Value_LDL,
    
    # Lp(a)-related variables
    n.snp.lpa       = NSNPs_LPA,
    estimate.mr.lpa = MR_beta_LPA,
    ci.lb.mr.lpa    = MR_CILower_LPA,
    ci.ub.mr.lpa    = MR_CIUpper_LPA,
    p.value.mr.lpa  = MR_P_Value_LPA
  ) %>%
  mutate(
    # Calculate standard errors from 95% confidence intervals
    std.error.mr.trl = (ci.ub.mr.trl - ci.lb.mr.trl) / (2 * 1.96),
    std.error.mr.ldl = (ci.ub.mr.ldl - ci.lb.mr.ldl) / (2 * 1.96),
    std.error.mr.lpa = (ci.ub.mr.lpa - ci.lb.mr.lpa) / (2 * 1.96)
  ) %>%
  # Rescale to 1-SD units based on UKBB lipid distributions
  mutate(
    # TRL: multiply by SD = 75.33189 mg/dL (0.07533189 mmol/L)
    estimate.mr.trl  = estimate.mr.trl * 0.07533189,
    std.error.mr.trl = std.error.mr.trl * 0.07533189,
    
    # LDL: multiply by SD = 432.5011 mg/dL (0.4331006 mmol/L)
    estimate.mr.ldl  = estimate.mr.ldl * 0.4331006,
    std.error.mr.ldl = std.error.mr.ldl * 0.4331006,
    
    # Lp(a): multiply by SD = 74.66814 nmol/L
    estimate.mr.lpa  = estimate.mr.lpa * 74.66814,
    std.error.mr.lpa = std.error.mr.lpa * 74.66814
  ) %>%
  select(
    protein_name,
    uniprot,
    gene,
    # LDL block
    n.snp.ldl,
    estimate.mr.ldl,
    std.error.mr.ldl,
    p.value.mr.ldl,
    # TRL block
    n.snp.trl,
    estimate.mr.trl,
    std.error.mr.trl,
    p.value.mr.trl,
    # Lp(a) block
    n.snp.lpa,
    estimate.mr.lpa,
    std.error.mr.lpa,
    p.value.mr.lpa
  )

# ==============================================================================
# MERGE WITH OBSERVATIONAL UKBB RESULTS
# ==============================================================================

# Select proteins with significant Bonferroni-corrected associations
# from the observational analysis (model_1)
combined <- left_join(
  model_1 %>%
    filter(
      pbonf.ukbb.ldl < 0.05 | 
        pbonf.ukbb.trl < 0.05 | 
        pbonf.ukbb.lpa < 0.05
    ) %>%
    select(-outcome),
  mr_results_clean,
  by = c("protein_name", "uniprot", "gene")
)

# ==============================================================================
# ASSESS AGREEMENT BETWEEN OBSERVATIONAL AND GENETIC ESTIMATES
# ==============================================================================

# Function to classify agreement between UKBB and MR estimates
# Agreement levels:
#   "++"  : Same direction, Bonferroni-significant (p < 0.05/1417)
#   "+"   : Same direction, nominally significant (p < 0.05)
#   "0"   : Not significant or disagreement (p >= 0.05)
#   "-"   : Opposite direction, nominally significant (p < 0.05)
#   "--"  : Opposite direction, Bonferroni-significant (p < 0.05/1417)

agree_label <- function(beta_ukbb, beta_mr, p_mr, n_tests = 1417) {
  p_bonf <- p_mr * n_tests
  case_when(
    is.na(beta_ukbb) | is.na(beta_mr) | is.na(p_mr) ~ NA_character_,
    p_mr > 0.05                                      ~ "0",
    sign(beta_ukbb) == sign(beta_mr) & p_bonf < 0.05 ~ "++",
    sign(beta_ukbb) == sign(beta_mr) & p_mr   < 0.05 ~ "+",
    sign(beta_ukbb) != sign(beta_mr) & p_bonf < 0.05 ~ "--",
    sign(beta_ukbb) != sign(beta_mr) & p_mr   < 0.05 ~ "-",
    TRUE ~ "0"
  )
}

# Apply agreement classification
combined <- combined %>%
  mutate(
    agree.ldl = agree_label(estimate.ukbb.ldl, estimate.mr.ldl, p.value.mr.ldl),
    agree.trl = agree_label(estimate.ukbb.trl, estimate.mr.trl, p.value.mr.trl),
    agree.lpa = agree_label(estimate.ukbb.lpa, estimate.mr.lpa, p.value.mr.lpa),
    # Set agreement to NA if UKBB association is not Bonferroni-significant
    agree.ldl = ifelse(p.value.ukbb.ldl > 0.05 / 1417, NA_character_, agree.ldl),
    agree.trl = ifelse(p.value.ukbb.trl > 0.05 / 1417, NA_character_, agree.trl),
    agree.lpa = ifelse(p.value.ukbb.lpa > 0.05 / 1417, NA_character_, agree.lpa)
  )

# Save results
save(combined, file = "combined.rda")

# ==============================================================================
# RESHAPE TO LONG FORMAT
# ==============================================================================

combined_long <- combined %>%
  # Pivot all lipid fraction columns to long format
  pivot_longer(
    cols = matches("^(estimate\\.ukbb|std\\.error\\.ukbb|p\\.value\\.ukbb|pbonf\\.ukbb|estimate\\.mr|std\\.error\\.mr|p\\.value\\.mr|agree)\\.(ldl|trl|lpa)$"),
    names_to = c(".value", "fraction"),
    names_pattern = "^(estimate\\.ukbb|std\\.error\\.ukbb|p\\.value\\.ukbb|pbonf\\.ukbb|estimate\\.mr|std\\.error\\.mr|p\\.value\\.mr|agree)\\.(ldl|trl|lpa)$"
  ) %>%
  mutate(
    # Standardize fraction labels
    fraction = recode(
      str_to_upper(fraction),
      "LDL" = "LDL",
      "TRL" = "TRL",
      "LPA" = "Lp(a)"
    ),
    # Determine direction of UKBB association
    direction.ukbb = case_when(
      estimate.ukbb > 0 ~ "positive",
      estimate.ukbb < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  # Retain only Bonferroni-significant UKBB associations
  filter(pbonf.ukbb < 0.05) %>%
  # Clean and factor the agreement variable
  mutate(
    agree = recode(
      str_replace_all(as.character(agree), "\\s+", ""),
      "++" = "++",
      "+"  = "+",
      "0"  = "0",
      "-"  = "-",
      "--" = "--",
      .default = NA_character_
    ),
    agree = factor(agree, levels = c("--", "-", "0", "+", "++"), ordered = TRUE)
  ) %>%
  select(
    protein_name,
    uniprot,
    gene,
    fraction,
    direction.ukbb,
    estimate.mr,
    std.error.mr,
    p.value.mr,
    agree
  ) %>%
  arrange(fraction, protein_name)

# Save results
save(combined_long, file = "combined_long.rda")

# ==============================================================================
# PREPARE OUTPUT TABLE
# ==============================================================================

ukbb_mr_valid_output <- combined_long %>%
  rename(
    `Protein - long name`   = protein_name,
    `Protein - short name`  = gene,
    `UniProt ID`            = uniprot,
    `Fraction`              = fraction,
    `UKBB association sign` = direction.ukbb,
    `Estimate MR`           = estimate.mr,
    `SE MR`                 = std.error.mr,
    `P MR`                  = p.value.mr,
    `Agreement`             = agree
  )

# ==============================================================================
# EXPORT RESULTS
# ==============================================================================

write_csv(ukbb_mr_valid_output, file = "ukbb_mr_valid_output.csv")