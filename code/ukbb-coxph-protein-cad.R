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
#   Required packages: dplyr, purrr, tidyr, stringr, survival, broom, tibble
#
# Goal:
#   Perform Cox proportional hazards regression for fraction-specific proteins
#   (LDL, TRL, Lp(a)) to assess associations with CAD outcomes in UK Biobank.
#   Models adjust for demographic and clinical covariates but exclude 
#   fraction-specific lipid measurements.
#
# Required input data:
#   - dataset.rda               # Main analysis dataset with proteins & covariates
#   - combined_long_signif.rda  # Protein metadata with fraction assignments
#
# Main outputs:
#   - cox_results_long.rda              # Raw results with metadata
#   - table-ukbb-coxph-proteins.csv            # Formatted results table
#
# ==============================================================================


# ==============================================================================
#  LOAD PACKAGES
# ==============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(survival)
library(broom)
library(tibble)


# ==============================================================================
#  DEFINE PROTEIN SIGNATURES BY FRACTION
# ==============================================================================
# Note: Assumes 'combined_long_signif' is already loaded in workspace
# Extract unique gene names for each lipoprotein fraction

protein_ldl <- combined_long_signif %>% 
  filter(fraction == "LDL") %>% 
  pull(gene) %>% 
  unique()

protein_trl <- combined_long_signif %>% 
  filter(fraction == "TRL") %>% 
  pull(gene) %>% 
  unique()

protein_lpa <- combined_long_signif %>% 
  filter(fraction == "Lp(a)") %>% 
  pull(gene) %>% 
  unique()


# ==============================================================================
#  DEFINE COVARIATES
# ==============================================================================
# Continuous covariates (will be standardized)
covars_num_cox <- c(
  "age",
  "bmi_imp",
  "conven_hba1c_imp",
  "sbp_mean_imp"
)

# Categorical covariates
covars_categorical <- c(
  "sex",
  "race_imp",
  "smoking_status_imp",
  "med_lipid"
)

# Note: TRL/LDL/Lp(a) fraction measurements intentionally excluded from models


# ==============================================================================
#  PREPARE MODELING DATASET
# ==============================================================================
# Extract protein column names (assumed to be columns 48:2965)
outcome_cols_all <- names(dataset)[48:2965]

# Define all required columns
base_needed <- c(
  "cad2_ukb_survtime",
  "cad2_ukb_ind",
  covars_num_cox,
  covars_categorical
)

# Prepare dataset: standardize continuous variables and ensure numeric types
df_model_base <- dataset %>%
  select(all_of(unique(c(base_needed, outcome_cols_all)))) %>%
  mutate(
    across(all_of(covars_num_cox), ~ as.numeric(scale(.))),
    across(all_of(outcome_cols_all), ~ suppressWarnings(as.numeric(.)))
  )


# ==============================================================================
#  ORGANIZE PROTEIN SIGNATURES
# ==============================================================================
# Create named list of gene signatures by fraction
genes_by_fraction <- list(
  "LDL"   = protein_ldl,
  "TRL"   = protein_trl,
  "Lp(a)" = protein_lpa
)

# Extract protein metadata for downstream joining
prot_meta <- combined_long_signif %>%
  distinct(protein_name, uniprot, gene, fraction, direction.ukbb)

# Filter to only genes present in the dataset
genes_present_by_fraction <- lapply(
  genes_by_fraction,
  function(gv) intersect(gv, colnames(df_model_base))
)


# ==============================================================================
#  DEFINE COX REGRESSION FUNCTION
# ==============================================================================
# Fit Cox PH model for a single protein gene
#
# Args:
#   gene_col: Character string of gene/protein column name
#
# Returns:
#   Tibble with gene name, log(HR), SE, and p-value

fit_cox_one <- function(gene_col) {
  # Construct formula: Surv ~ scaled_protein + covariates
  rhs <- paste(
    c(
      sprintf("scale(`%s`)", gene_col),
      covars_num_cox,
      covars_categorical
    ),
    collapse = " + "
  )
  
  fml <- as.formula(paste0("Surv(cad2_ukb_survtime, cad2_ukb_ind) ~ ", rhs))
  
  # Fit Cox proportional hazards model
  mdl <- coxph(fml, data = df_model_base)
  
  # Extract protein coefficient using regex to handle special characters
  term_regex <- paste0(
    "^scale\\((`)?",
    str_replace_all(gene_col, "([\\W])", "\\\\\\1"),
    "(\\1)?\\)$"
  )
  
  td <- tidy(mdl) %>%
    filter(str_detect(term, term_regex)) %>%
    transmute(
      gene                  = gene_col,
      estimate.ukbb.protein = estimate,      # log(HR) per 1 SD increase
      std.ukbb.protein      = std.error,
      p.value.protein       = p.value
    )
  
  # Handle case where protein term was dropped from model
  if (nrow(td) == 0) {
    td <- tibble(
      gene                  = gene_col,
      estimate.ukbb.protein = NA_real_,
      std.ukbb.protein      = NA_real_,
      p.value.protein       = NA_real_
    )
  }
  
  return(td)
}


# ==============================================================================
#  RUN COX MODELS ACROSS ALL FRACTIONS
# ==============================================================================
# Apply Cox regression to each protein within each fraction

cox_results_list <- imap(genes_present_by_fraction, function(genes_vec, frac_label) {
  # Skip empty fractions
  if (length(genes_vec) == 0) return(tibble())
  
  # Fit models for all genes in this fraction
  map_dfr(genes_vec, function(g) {
    tryCatch(
      fit_cox_one(g),
      error = function(e) {
        # Return NA row if model fails
        tibble(
          gene                  = g,
          estimate.ukbb.protein = NA_real_,
          std.ukbb.protein      = NA_real_,
          p.value.protein       = NA_real_
        )
      }
    )
  }) %>%
    mutate(fraction = frac_label, .before = 1)
})

# Combine results from all fractions
cox_results_raw <- bind_rows(cox_results_list)


# ==============================================================================
#  JOIN WITH PROTEIN METADATA
# ==============================================================================
# Add protein names, UniProt IDs, and association directions

cox_results_long <- cox_results_raw %>%
  left_join(prot_meta, by = c("gene", "fraction")) %>%
  select(
    protein_name,
    uniprot,
    gene,
    fraction,
    direction.ukbb,
    estimate.ukbb.protein,
    std.ukbb.protein,
    p.value.protein
  ) %>%
  arrange(fraction, gene)

# Save intermediate results
save(cox_results_long, file = "cox_results_long.rda")


# ==============================================================================
#  FORMAT OUTPUT TABLE
# ==============================================================================
# Rename columns for publication-ready output

cox_results_long_output <- cox_results_long %>%
  rename(
    `Protein - long name`   = protein_name,
    `Protein - short name`  = gene,
    `UniProt ID`            = uniprot,
    `Fraction`              = fraction,
    `UKBB association sign` = direction.ukbb,
    `log(HR) UKBB`          = estimate.ukbb.protein,
    `SE UKBB`               = std.ukbb.protein,
    `P UKBB`                = p.value.protein
  )


# ==============================================================================
#  SAVE RESULTS
# ==============================================================================
# Save formatted results


write_csv(cox_results_long_output, file = "PATH/table-ukbb-coxph-proteins.csv")
