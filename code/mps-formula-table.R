# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-26
# Last modified: 2025-11-28
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, stringr, tibble, purrr
#
# Goal:
#   Build formula strings for multi-protein scores (MPS) from ridge regression
#   models. Converts ridge_boot_* objects (LDL, TRL, LPA) into human-readable
#   mathematical formulas showing protein coefficients.
#
# Required input data:
#   - ridge_boot_LDL_full
#   - ridge_boot_TRL_full
#   - ridge_boot_LPA_full
#
# Main outputs:
#   - mps-formulas.csv 
#   - mps-formulas.rda 
#
# ==============================================================================


# ==============================================================================
# PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})


# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================

#' Convert ridge regression object to formula string
#'
#' @param ridge_obj Ridge regression bootstrap object with summary_df component
#' @param signature_label Character label for the signature (e.g., "LDL")
#' @param coef_col Character name of coefficient column (default: "coef_full_lambda")
#'
#' @return Tibble with columns: Signature, n_terms, Formula
ridge_to_formula <- function(ridge_obj, 
                             signature_label, 
                             coef_col = "coef_full_lambda") {
  
  # Validate input object
  if (is.null(ridge_obj) || is.null(ridge_obj$summary_df)) {
    return(tibble(
      Signature = signature_label,
      n_terms = NA_integer_,
      Formula = NA_character_
    ))
  }
  
  df <- ridge_obj$summary_df
  
  # Check required columns
  if (!all(c("protein", coef_col) %in% names(df))) {
    return(tibble(
      Signature = signature_label,
      n_terms = NA_integer_,
      Formula = NA_character_
    ))
  }
  
  # Filter non-zero coefficients
  df_terms <- df %>%
    transmute(
      protein_clean = str_remove(protein, "_rz$"),
      coef = .data[[coef_col]]
    ) %>%
    filter(!is.na(coef), coef != 0)
  
  # Handle empty result
  if (nrow(df_terms) == 0) {
    return(tibble(
      Signature = signature_label,
      n_terms = 0L,
      Formula = "0"
    ))
  }
  
  # Helper: format absolute coefficient value
  abs_fmt <- function(x) {
    formatC(abs(x), digits = 4, format = "f")
  }
  
  # Build first term (without leading sign for positive)
  first <- df_terms[1, ]
  first_prefix <- if (first$coef < 0) "- " else ""
  first_term <- paste0(
    first_prefix,
    abs_fmt(first$coef),
    " * ",
    first$protein_clean
  )
  
  # Build remaining terms with explicit +/- signs
  rest <- NULL
  if (nrow(df_terms) > 1) {
    rest <- df_terms[-1, ] %>%
      mutate(
        sign = ifelse(coef < 0, "- ", "+ "),
        term = paste0(sign, abs_fmt(coef), " * ", protein_clean)
      ) %>%
      pull(term)
  }
  
  # Combine all terms
  formula_str <- paste(c(first_term, rest), collapse = " ")
  
  tibble(
    Signature = signature_label,
    n_terms = nrow(df_terms),
    Formula = formula_str
  )
}


# ==============================================================================
# BUILD MPS FORMULAS
# ==============================================================================

# Define mapping of signature labels to R object names
objects_map <- c(
  "LDL" = "ridge_boot_LDL_full",
  "TRL" = "ridge_boot_TRL_full",
  "Lp(a)" = "ridge_boot_LPA_full"
)

# Process all signatures
mps_formulas <- imap_dfr(objects_map, function(obj_name, sig_label) {
  ridge_obj <- get0(obj_name, ifnotfound = NULL)
  ridge_to_formula(ridge_obj, sig_label)
})

# Save raw table 
save(mps_formulas, file = "mps_formulas.rda")

# Rename columns for output
mps_formulas_output <- mps_formulas %>%
  rename(
    Fraction = Signature,
    `N of proteins` = n_terms,
    Formula = Formula
  )


# ==============================================================================
# SAVE RESULTS
# ==============================================================================

write.csv(mps_formulas_output, "mps-formulas.csv", row.names = FALSE)
