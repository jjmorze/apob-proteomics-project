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
#   Required packages: readxl, dplyr, readr
#
# Goal:
#   Integrate MESA validation results with UK Biobank Cox proportional hazards
#   analysis for incident CHD. Merge proteomic associations from both cohorts,
#   harmonize variable names, and export a combined results table for 
#   cross-cohort validation of protein biomarkers.
#
# Required input data:
#   - MESA Validation 2025-8-17.xlsx (sheet: "Incident CHD") 
#   - cox_results_long 
#
# Main outputs:
#   - table-mesa-coxph-proteins.csv 
#   - combined_ukbb_mesa_coxph.rda
#
# ==============================================================================

# ==============================================================================
#  LOAD PACKAGES
# ==============================================================================

library(tidyverse)
library(readxl)

# ==============================================================================
#  LOAD AND PREPARE MESA RESULTS
# ==============================================================================

# Load MESA validation results for incident CHD
mesa_coxph <- read_excel(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Desktop – Jakub's Mac mini/Projects/UK Biobank/Olink Proteomics/natmed-apob-proteomics/data/mesa-validation/MESA Validation 2025-8-17.xlsx",
  sheet = "Incident CHD"
)

# Rename columns for clarity and consistency
mesa_coxph <- mesa_coxph %>%
  dplyr::rename(
    outcome.type     = outcome,
    estimate.mesa    = estimate,
    std.error.mesa   = std.error,
    p.value.mesa     = p.value,
    n.mesa           = n,
    n.event.mesa     = nevent,
    gene             = Assay,
    uniprot          = UniProt
  )


# ==============================================================================
#  MERGE UK BIOBANK AND MESA RESULTS
# ==============================================================================

# Combine UKBB Cox results with MESA Hard CHD outcomes
combined_ukbb_mesa_coxph <- left_join(
  cox_results_long,
  mesa_coxph %>% filter(outcome.type == "Hard CHD"),
  by = c("gene", "uniprot")
)

# Save raw table

save(combined_ukbb_mesa_coxph, file = "combined_ukbb_mesa_coxph.rda")

# ==============================================================================
#  FORMAT OUTPUT TABLE
# ==============================================================================

# Select and rename columns for final output table
combined_ukbb_mesa_coxph_output <- combined_ukbb_mesa_coxph %>%
  dplyr::select(
    protein_name,
    uniprot,
    gene,
    fraction,
    direction.ukbb,
    n.mesa,
    n.event.mesa,
    estimate.mesa,
    std.error.mesa,
    p.value.mesa
  ) %>%
  dplyr::rename(
    `Protein - long name`   = protein_name,
    `Protein - short name`  = gene,
    `UniProt ID`            = uniprot,
    `Fraction`              = fraction,
    `UKBB association sign` = direction.ukbb,
    `n MESA`                = n.mesa,
    `Events MESA`           = n.event.mesa,
    `log(HR) MESA`          = estimate.mesa,
    `SE MESA`               = std.error.mesa,
    `P MESA`                = p.value.mesa
  )


# ==============================================================================
#  SAVE RESULTS
# ==============================================================================

# Export combined UKBB-MESA proteomic results
write_csv(
  combined_ukbb_mesa_coxph_output,
  file = "PATH/table-mesa-coxph-proteins.csv")
)

# ==============================================================================
#  END OF SCRIPT
# ==============================================================================