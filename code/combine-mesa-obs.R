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
#   Validate UK Biobank proteomic associations against MESA cohort data by:
#   1) Loading MESA validation results
#   2) Harmonizing variable names and fraction labels
#   3) Merging with UK Biobank significant associations
#   4) Assessing directional agreement between cohorts
#   5) Exporting consolidated validation table
#
# Required input data:
#   - MESA Validation 2025-8-17.xlsx # RESTRICTED
#   - combined_long.rda
#
# Main outputs:
#   - mesa-validation-with-ukbb-output.csv
#   - mesa_validation_with_ukbb.rda
#   - combined_long_signif.rda
#
# ==============================================================================

# ==============================================================================
#  PACKAGE LOADING
# ==============================================================================

library(readxl)
library(dplyr)
library(readr)

# ==============================================================================
# 1. LOAD MESA VALIDATION DATA
# ==============================================================================

mesa_validation <- read_excel(
  "PATH/MESA Validation 2025-8-17.xlsx"
)

# ==============================================================================
# 2. RENAME AND HARMONIZE VARIABLES
# ==============================================================================

mesa_validation <- mesa_validation %>% 
  dplyr::rename(
    fraction      = term,
    estimate.mesa = estimate,
    std.error.mesa = std.error,
    statistic.mesa = statistic,
    p.value.mesa  = p.value,
    n.mesa        = n,
    gene          = Assay,
    uniprot       = UniProt
  ) %>% 
  mutate(
    fraction = case_when(
      fraction == "ldl_z" ~ "LDL",
      fraction == "trl_z" ~ "TRL",
      fraction == "lpa_z" ~ "Lp(a)",
      TRUE ~ fraction
    )
  )

# ==============================================================================
# 3. FILTER UKBB SIGNIFICANT ASSOCIATIONS
# ==============================================================================

combined_long_signif <- combined_long %>% 
  filter(agree == "++") %>% 
  dplyr::select(protein_name, uniprot, gene, fraction, direction.ukbb, agree)

# Save the file with raw formatting
save(
  combined_long_signif, 
  file = "PATH/combined_long_signif.rda"
)

# ==============================================================================
# 4. MERGE UKBB AND MESA DATA
# ==============================================================================

mesa_validation_with_ukbb <- left_join(
  combined_long_signif, 
  mesa_validation,
  by = c("gene", "fraction")
)

# ==============================================================================
# 5. ASSESS DIRECTIONAL AGREEMENT WITH MESA
# ==============================================================================

# Function to determine agreement between UKBB and MESA results
agree_mesa <- function(direction_ukbb, beta_mesa, p_mesa, alpha = 0.05) {
  # Convert UKBB direction to numeric sign
  s_ukbb <- ifelse(direction_ukbb == "positive", 1L,
                   ifelse(direction_ukbb == "negative", -1L, NA_integer_))
  
  # Determine agreement category
  dplyr::case_when(
    is.na(s_ukbb) | is.na(beta_mesa) | is.na(p_mesa) ~ NA_character_,
    p_mesa > alpha                                    ~ "0",
    sign(beta_mesa) == s_ukbb                         ~ "+",
    sign(beta_mesa) != s_ukbb                         ~ "-",
    TRUE                                              ~ NA_character_
  )
}

# Apply agreement function
mesa_validation_with_ukbb <- mesa_validation_with_ukbb %>%
  mutate(
    agree.mesa = agree_mesa(direction.ukbb, estimate.mesa, p.value.mesa)
  )

# Save the file with raw formatting
save(
  mesa_validation_with_ukbb, 
  file = "PATH/mesa_validation_with_ukbb.rda"
)

# ==============================================================================
# 6. PREPARE AND EXPORT FINAL OUTPUT TABLE
# ==============================================================================

mesa_validation_with_ukbb_output <- mesa_validation_with_ukbb %>%
  dplyr::select(
    protein_name, 
    gene, 
    uniprot.x, 
    fraction, 
    direction.ukbb, 
    estimate.mesa, 
    std.error.mesa, 
    p.value.mesa, 
    n.mesa, 
    agree.mesa
  ) %>% 
  rename(
    `Protein - long name`   = protein_name,
    `Protein - short name`  = gene,
    `UniProt ID`            = uniprot.x,
    `Fraction`              = fraction,
    `UKBB association sign` = direction.ukbb,
    `Estimate MESA`         = estimate.mesa,
    `SE MESA`               = std.error.mesa,
    `P MESA`                = p.value.mesa,
    `N of participants`     = n.mesa,
    `Agreement`             = agree.mesa
  )

# Save the validation table
write_csv(
  mesa_validation_with_ukbb_output, 
  file = "PATH/mesa-validation-with-ukbb-output.csv"
)