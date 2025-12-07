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
#   Required packages: dplyr, purrr, broom, tidyr, readr, tibble
#
# Goal:
#   Perform sensitivity analyses examining associations between lipoprotein
#   fractions and plasma proteins using:
#   1) NMR-derived lipoprotein particle concentrations (VLDL-P, IDL+LDL-P)
#   2) Olink-measured LPA protein concentration
#   Both analyses adjust for conventional lipoprotein measures and covariates,
#   with Bonferroni correction for multiple testing.
#
# Required input data:
#   - dataset.rda                   
#   - dataset_side1.rda
#   - dataset_side2.rda  
#
# Main outputs:
#   - model_2.rda                                  # NMR sensitivity results
#   - model_3.rda                                  # Olink LPA sensitivity results
#   - table-ukbb-lm-protein-sensitvity-nmr.csv           # NMR results table
#   - table-ukbb-lm-protein-sensitvity-olinklpa.csv      # Olink LPA results table
#
# ==============================================================================

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(readr)

# ==============================================================================
# SENSITIVITY ANALYSIS 1: NMR-DERIVED LIPOPROTEIN PARTICLES
# ==============================================================================

# ------------------------------------------------------------------------------
# Prepare NMR particle data
# ------------------------------------------------------------------------------

load("~/dataset_side2.rda")


# Merge NMR data with main dataset
dataset <- left_join(dataset, dataset_side2)

# Filter to participants with complete NMR data
dataset_sub1 <- dataset %>% 
  filter(!is.na(VLDL_P) & !is.na(IDL_LDL_P))

# ------------------------------------------------------------------------------
# Define model components for NMR analysis
# ------------------------------------------------------------------------------

outcome_cols_s1 <- names(dataset_sub1)[48:2965]

covars_num_to_scale_s1 <- c(
  "conven_lpa_corr",
  "IDL_LDL_P",
  "VLDL_P",
  "age",
  "bmi_imp",
  "conven_hba1c_imp",
  "sbp_mean_imp"
)

covars_categorical_s1 <- c(
  "sex",
  "race_imp",
  "smoking_status_imp",
  "med_lipid"
)

vars_needed_s1 <- c(
  covars_num_to_scale_s1,
  covars_categorical_s1,
  outcome_cols_s1
)

# ------------------------------------------------------------------------------
# Prepare modeling dataset
# ------------------------------------------------------------------------------

df_model_s1 <- dataset_sub1 %>%
  select(all_of(vars_needed_s1)) %>%
  mutate(
    across(all_of(covars_num_to_scale_s1), ~ as.numeric(scale(.))),
    across(all_of(outcome_cols_s1), as.numeric)
  )

# ------------------------------------------------------------------------------
# Define term mapping for clean output names
# ------------------------------------------------------------------------------

term_map_s1 <- c(
  "conven_lpa_corr" = "lpa",
  "IDL_LDL_P"       = "ldl",
  "VLDL_P"          = "trl"
)

# ------------------------------------------------------------------------------
# Model fitting function
# ------------------------------------------------------------------------------

fit_one_s1 <- function(outcome) {
  
  # Construct formula
  fml <- as.formula(
    paste0(
      "scale(`", outcome, "`) ~ ",
      paste(
        c(
          names(term_map_s1),
          covars_num_to_scale_s1[!covars_num_to_scale_s1 %in% names(term_map_s1)],
          covars_categorical_s1
        ),
        collapse = " + "
      )
    )
  )
  
  # Fit model
  mdl <- lm(fml, data = df_model_s1)
  
  # Extract and reshape results
  td <- broom::tidy(mdl) %>%
    dplyr::filter(term %in% names(term_map_s1)) %>%
    dplyr::mutate(key = term_map_s1[term]) %>%
    dplyr::transmute(key, estimate, std.error, p.value)
  
  wide <- tidyr::pivot_wider(
    td,
    names_from = key,
    values_from = c(estimate, std.error, p.value),
    names_glue = "{.value}.ukbb.nmr.{key}"
  )
  
  # Apply Bonferroni correction (capped at 1.0)
  for (k in unname(term_map_s1)) {
    pcol <- paste0("p.value.ukbb.nmr.", k)
    wide[[paste0("pbonf.ukbb.nmr.", k)]] <- pmin(wide[[pcol]] * 1417, 1.0)
  }
  
  wide %>%
    dplyr::mutate(outcome = outcome, .before = 1)
}

# ------------------------------------------------------------------------------
# Run models across all protein outcomes
# ------------------------------------------------------------------------------

model_2 <- purrr::map_dfr(outcome_cols_s1, fit_one_s1) %>%
  dplyr::mutate(name_uniprot_gene = name_uniprot_gene, .before = 1) %>%
  tidyr::unnest_wider(name_uniprot_gene)

save(model_2, file = "model_2.rda")

# ------------------------------------------------------------------------------
# Count Bonferroni-significant associations
# ------------------------------------------------------------------------------

model_2 %>%
  dplyr::summarise(
    lpa_sig = sum(pbonf.ukbb.nmr.lpa < 0.05, na.rm = TRUE),
    ldl_sig = sum(pbonf.ukbb.nmr.ldl < 0.05, na.rm = TRUE),
    trl_sig = sum(pbonf.ukbb.nmr.trl < 0.05, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to = "fraction",
    values_to = "n_proteins"
  )

# ------------------------------------------------------------------------------
# Format output table
# ------------------------------------------------------------------------------

model_2_output <- model_2 %>%
  select(
    protein_name,
    uniprot,
    gene,
    estimate.ukbb.nmr.ldl,
    std.error.ukbb.nmr.ldl,
    p.value.ukbb.nmr.ldl,
    pbonf.ukbb.nmr.ldl,
    estimate.ukbb.nmr.trl,
    std.error.ukbb.nmr.trl,
    p.value.ukbb.nmr.trl,
    pbonf.ukbb.nmr.trl,
    estimate.ukbb.nmr.lpa,
    std.error.ukbb.nmr.lpa,
    p.value.ukbb.nmr.lpa,
    pbonf.ukbb.nmr.lpa
  ) %>%
  rename(
    `Protein - long name`  = protein_name,
    `Protein - short name` = gene,
    `UniProt ID`           = uniprot,
    `Estimate LDL`         = estimate.ukbb.nmr.ldl,
    `SE LDL`               = std.error.ukbb.nmr.ldl,
    `P LDL`                = p.value.ukbb.nmr.ldl,
    `Padj LDL`             = pbonf.ukbb.nmr.ldl,
    `Estimate TRL`         = estimate.ukbb.nmr.trl,
    `SE TRL`               = std.error.ukbb.nmr.trl,
    `P TRL`                = p.value.ukbb.nmr.trl,
    `Padj TRL`             = pbonf.ukbb.nmr.trl,
    `Estimate Lp(a)`       = estimate.ukbb.nmr.lpa,
    `SE Lp(a)`             = std.error.ukbb.nmr.lpa,
    `P Lp(a)`              = p.value.ukbb.nmr.lpa,
    `Padj Lp(a)`           = pbonf.ukbb.nmr.lpa
  ) %>%
  select(
    `Protein - long name`,
    `Protein - short name`,
    `UniProt ID`,
    everything()
  ) %>%
  ungroup() %>%
  tibble::as_tibble()

write_csv(model_2_output, file = "table-ukbb-lm-protein-sensitvity-nmr.csv")

# ==============================================================================
# SENSITIVITY ANALYSIS 2: OLINK LPA PROTEIN
# ==============================================================================

# ------------------------------------------------------------------------------
# Prepare Olink LPA data
# ------------------------------------------------------------------------------

load("~/dataset_side1_092625.rda")

dataset_sub2 <- left_join(dataset, dataset_side1)

# ------------------------------------------------------------------------------
# Define model components for Olink LPA analysis
# ------------------------------------------------------------------------------

outcome_cols_s2 <- names(dataset_sub2)[48:2965]

covars_num_to_scale_s2 <- c(
  "LPA",
  "conven_ldlp",
  "conven_trlp",
  "age",
  "bmi_imp",
  "conven_hba1c_imp",
  "sbp_mean_imp"
)

covars_categorical_s2 <- c(
  "sex",
  "race_imp",
  "smoking_status_imp",
  "med_lipid"
)

vars_needed_s2 <- c(
  covars_num_to_scale_s2,
  covars_categorical_s2,
  outcome_cols_s2
)

# ------------------------------------------------------------------------------
# Prepare modeling dataset
# ------------------------------------------------------------------------------

df_model_s2 <- dataset_sub2 %>%
  select(all_of(vars_needed_s2)) %>%
  mutate(
    across(all_of(covars_num_to_scale_s2), ~ as.numeric(scale(.))),
    across(all_of(outcome_cols_s2), as.numeric)
  )

# ------------------------------------------------------------------------------
# Define term mapping for clean output names
# ------------------------------------------------------------------------------

term_map_s2 <- c(
  "LPA"         = "lpa",
  "conven_ldlp" = "ldl",
  "conven_trlp" = "trl"
)

# ------------------------------------------------------------------------------
# Model fitting function
# ------------------------------------------------------------------------------

fit_one_s2 <- function(outcome) {
  
  # Construct formula
  fml <- as.formula(
    paste0(
      "scale(`", outcome, "`) ~ ",
      paste(
        c(
          names(term_map_s2),
          covars_num_to_scale_s2[!covars_num_to_scale_s2 %in% names(term_map_s2)],
          covars_categorical_s2
        ),
        collapse = " + "
      )
    )
  )
  
  # Fit model
  mdl <- lm(fml, data = df_model_s2)
  
  # Extract and reshape results
  td <- tidy(mdl) %>%
    filter(term %in% names(term_map_s2)) %>%
    mutate(key = term_map_s2[term]) %>%
    transmute(key, estimate, std.error, p.value)
  
  wide <- td %>%
    pivot_wider(
      names_from = key,
      values_from = c(estimate, std.error, p.value),
      names_glue = "{.value}.ukbb.olinklpa.{key}"
    )
  
  # Apply Bonferroni correction (capped at 1.0)
  for (k in unname(term_map_s2)) {
    pcol <- paste0("p.value.ukbb.olinklpa.", k)
    wide[[paste0("pbonf.ukbb.olinklpa.", k)]] <- pmin(wide[[pcol]] * 1417, 1.0)
  }
  
  wide %>%
    mutate(outcome = outcome, .before = 1)
}

# ------------------------------------------------------------------------------
# Run models across all protein outcomes
# ------------------------------------------------------------------------------

model_3 <- purrr::map_dfr(outcome_cols_s2, fit_one_s2) %>%
  dplyr::mutate(name_uniprot_gene = name_uniprot_gene, .before = 1) %>%
  tidyr::unnest_wider(name_uniprot_gene)

save(model_3, file = "model_3.rda")

# ------------------------------------------------------------------------------
# Count Bonferroni-significant associations
# ------------------------------------------------------------------------------

model_3 %>%
  dplyr::summarise(
    lpa_sig = sum(pbonf.ukbb.olinklpa.lpa < 0.05, na.rm = TRUE),
    ldl_sig = sum(pbonf.ukbb.olinklpa.ldl < 0.05, na.rm = TRUE),
    trl_sig = sum(pbonf.ukbb.olinklpa.trl < 0.05, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to = "fraction",
    values_to = "n_proteins"
  )

# ------------------------------------------------------------------------------
# Format output table
# ------------------------------------------------------------------------------

model_3_output <- model_3 %>%
  select(
    protein_name,
    uniprot,
    gene,
    estimate.ukbb.olinklpa.ldl,
    std.error.ukbb.olinklpa.ldl,
    p.value.ukbb.olinklpa.ldl,
    pbonf.ukbb.olinklpa.ldl,
    estimate.ukbb.olinklpa.trl,
    std.error.ukbb.olinklpa.trl,
    p.value.ukbb.olinklpa.trl,
    pbonf.ukbb.olinklpa.trl,
    estimate.ukbb.olinklpa.lpa,
    std.error.ukbb.olinklpa.lpa,
    p.value.ukbb.olinklpa.lpa,
    pbonf.ukbb.olinklpa.lpa
  ) %>%
  rename(
    `Protein - long name`  = protein_name,
    `Protein - short name` = gene,
    `UniProt ID`           = uniprot,
    `Estimate LDL`         = estimate.ukbb.olinklpa.ldl,
    `SE LDL`               = std.error.ukbb.olinklpa.ldl,
    `P LDL`                = p.value.ukbb.olinklpa.ldl,
    `Padj LDL`             = pbonf.ukbb.olinklpa.ldl,
    `Estimate TRL`         = estimate.ukbb.olinklpa.trl,
    `SE TRL`               = std.error.ukbb.olinklpa.trl,
    `P TRL`                = p.value.ukbb.olinklpa.trl,
    `Padj TRL`             = pbonf.ukbb.olinklpa.trl,
    `Estimate Lp(a)`       = estimate.ukbb.olinklpa.lpa,
    `SE Lp(a)`             = std.error.ukbb.olinklpa.lpa,
    `P Lp(a)`              = p.value.ukbb.olinklpa.lpa,
    `Padj Lp(a)`           = pbonf.ukbb.olinklpa.lpa
  ) %>%
  select(
    `Protein - long name`,
    `Protein - short name`,
    `UniProt ID`,
    everything()
  ) %>%
  ungroup() %>%
  tibble::as_tibble()

write_csv(model_3_output, file = "table-ukbb-lm-protein-sensitvity-olinklpa.csv")
