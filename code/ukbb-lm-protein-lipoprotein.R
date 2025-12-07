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
#   Required packages: dplyr, purrr, broom, tidyr, readr
#
# Goal:
#   Perform linear regression analysis of UK Biobank Olink proteomic data
#   against lipoprotein fractions (Lp(a), LDL, TRL), adjusting for key
#   covariates. The workflow:
#   1) Loads pre-built analysis dataset with proteomics and covariates
#   2) Performs PCA to determine optimal Bonferroni correction factor
#   3) Fits linear models for each protein outcome against lipoprotein fractions
#   4) Applies Bonferroni correction based on effective dimensionality
#   5) Exports formatted results table
#
# Required input data:
#   - dataset.rda # RESTRICTED
#   - name_uniprot_gene_pathway.xlsx
#
# Outputs:
#   - model_1.rda (raw model results)
#   - table-ukbb-lm-protein-lipoprotein.csv (formatted results table)
#
# ==============================================================================

# Load required libraries
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(readr)

# Load  UK Biobank dataset
load("~PATH/dataset.rda")


# ==============================================================================
# DETERMINE EFFECTIVE DIMENSIONALITY VIA PCA
# ==============================================================================

# Run PCA on protein variables (columns 48:2965)
pca_res <- prcomp(
  dataset[, 48:2965],
  center = TRUE,
  scale. = TRUE
)

# Calculate proportion of variance explained by each PC
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Calculate cumulative variance
cum_var <- cumsum(var_explained)

# Determine number of PCs needed to explain ≥90% of variance
n_pcs_90 <- which(cum_var >= 0.90)[1]

# Report result
cat("Number of PCs explaining at least 90% of variance:", n_pcs_90, "\n")
# Result: 1417 PCs (used for Bonferroni correction)


# ==============================================================================
# PREPARE MODELING DATASET
# ==============================================================================

# Define outcome variables (protein columns)
outcome_cols <- names(dataset)[48:2965]

# Define continuous covariates (to be scaled)
covars_num_to_scale <- c(
  "conven_lpa_corr",      # Lp(a) concentration (corrected)
  "conven_ldlp",          # LDL particle concentration
  "conven_trlp",          # TRL particle concentration
  "age",                  # Age at assessment
  "bmi_imp",              # Body mass index (imputed)
  "conven_hba1c_imp",     # HbA1c (imputed)
  "sbp_mean_imp"          # Systolic blood pressure (imputed)
)

# Define categorical covariates
covars_categorical <- c(
  "sex",                  # Biological sex
  "race_imp",             # Race/ethnicity (imputed)
  "smoking_status_imp",   # Smoking status (imputed)
  "med_lipid"             # Lipid-lowering medication use
)

# Combine all required variables
vars_needed <- c(covars_num_to_scale, covars_categorical, outcome_cols)

# Create modeling dataset with standardized covariates
df_model <- dataset %>%
  select(all_of(vars_needed)) %>%
  mutate(
    # Standardize continuous covariates (mean=0, SD=1)
    across(all_of(covars_num_to_scale), ~ as.numeric(scale(.))),
    # Ensure outcomes are numeric (will be scaled within model formula)
    across(all_of(outcome_cols), as.numeric)
  )


# ==============================================================================
# DEFINE MODELING FRAMEWORK
# ==============================================================================

# Create mapping for clean output variable names
term_map <- c(
  "conven_lpa_corr" = "lpa",
  "conven_ldlp"     = "ldl",
  "conven_trlp"     = "trl"
)

# Define function to fit linear model for one protein outcome
fit_one <- function(outcome) {
  
  # Construct formula: scale(protein) ~ lipoprotein fractions + covariates
  fml <- as.formula(
    paste0(
      "scale(`", outcome, "`) ~ ",
      paste(
        c(
          names(term_map),
          covars_num_to_scale[!covars_num_to_scale %in% names(term_map)],
          covars_categorical
        ),
        collapse = " + "
      )
    )
  )
  
  # Fit linear model
  mdl <- lm(fml, data = df_model)
  
  # Extract and format results for lipoprotein fractions only
  td <- tidy(mdl) %>%
    filter(term %in% names(term_map)) %>%
    mutate(key = term_map[term]) %>%
    transmute(key, estimate, std.error, p.value)
  
  # Pivot to wide format with informative column names
  wide <- td %>%
    pivot_wider(
      names_from = key,
      values_from = c(estimate, std.error, p.value),
      names_glue = "{.value}.ukbb.{key}"
    )
  
  # Apply Bonferroni correction (capped at 1.0)
  for (k in unname(term_map)) {
    pcol <- paste0("p.value.ukbb.", k)
    wide[[paste0("pbonf.ukbb.", k)]] <- pmin(wide[[pcol]] * 1417, 1.0)
  }
  
  # Add outcome identifier
  wide %>%
    mutate(outcome = outcome, .before = 1)
}


# ==============================================================================
# RUN MODELS ACROSS ALL PROTEINS
# ==============================================================================

# Fit models for all protein outcomes
model_1 <- map_dfr(outcome_cols, fit_one) %>%
  mutate(name_uniprot_gene_pathway = name_uniprot_gene_pathway, .before = 1) %>%
  unnest_wider(name_uniprot_gene)  # Expands to: protein_name, uniprot, gene

# Save raw model results
save(model_1, file = "model_1.rda")


# ==============================================================================
# SUMMARIZE SIGNIFICANT ASSOCIATIONS
# ==============================================================================

# Count Bonferroni-significant proteins for each lipoprotein fraction
significance_summary <- model_1 %>%
  summarise(
    lpa_sig = sum(pbonf.ukbb.lpa < 0.05, na.rm = TRUE),
    ldl_sig = sum(pbonf.ukbb.ldl < 0.05, na.rm = TRUE),
    trl_sig = sum(pbonf.ukbb.trl < 0.05, na.rm = TRUE)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "fraction",
    values_to = "n_proteins"
  )

print(significance_summary)
# Expected output:
# fraction    n_proteins
# lpa_sig            221
# ldl_sig           1058
# trl_sig           2020


# ==============================================================================
# FORMAT OUTPUT TABLE
# ==============================================================================

# Create publication-ready results table
model_1_output <- model_1 %>%
  select(
    protein_name,
    uniprot,
    gene,
    # LDL associations
    estimate.ukbb.ldl,
    std.error.ukbb.ldl,
    p.value.ukbb.ldl,
    pbonf.ukbb.ldl,
    # TRL associations
    estimate.ukbb.trl,
    std.error.ukbb.trl,
    p.value.ukbb.trl,
    pbonf.ukbb.trl,
    # Lp(a) associations
    estimate.ukbb.lpa,
    std.error.ukbb.lpa,
    p.value.ukbb.lpa,
    pbonf.ukbb.lpa
  ) %>%
  rename(
    `Protein - long name`  = protein_name,
    `Protein - short name` = gene,
    `UniProt ID`           = uniprot,
    `Estimate LDL`         = estimate.ukbb.ldl,
    `SE LDL`               = std.error.ukbb.ldl,
    `P LDL`                = p.value.ukbb.ldl,
    `Padj LDL`             = pbonf.ukbb.ldl,
    `Estimate TRL`         = estimate.ukbb.trl,
    `SE TRL`               = std.error.ukbb.trl,
    `P TRL`                = p.value.ukbb.trl,
    `Padj TRL`             = pbonf.ukbb.trl,
    `Estimate Lp(a)`       = estimate.ukbb.lpa,
    `SE Lp(a)`             = std.error.ukbb.lpa,
    `P Lp(a)`              = p.value.ukbb.lpa,
    `Padj Lp(a)`           = pbonf.ukbb.lpa
  ) %>%
  select(
    `Protein - long name`,
    `Protein - short name`,
    `UniProt ID`,
    everything()
  ) %>%
  ungroup() %>%
  tibble::as_tibble()

# Export formatted results table
write_csv(model_1_output, file = "table-ukbb-lm-protein-lipoprotein.csv")

