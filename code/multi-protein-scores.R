# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2024-09-26
# Last modified: 2025-11-28
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: dplyr, glmnet
#
# Goal:
#   Perform ridge regression with bootstrap validation for lipoprotein-specific
#   protein panels:
#   1) Extract fraction-specific proteins (LDL, TRL, Lp(a))
#   2) Z-score transform protein measurements
#   3) Generate residualized lipoprotein outcomes (adjusted for cross-fraction effects)
#   4) Run ridge regression with nested cross-validation and bootstrap resampling
#   5) Compute weighted protein scores for each fraction using full-data and
#      bootstrap-median coefficients
#   6) Merge all scores into a unified dataset
#
# Required input data:
#   - combined_long_signif.rda
#   - dataset.rda
#
# Main outputs:
#   - coef_summary_LDL.csv
#   - coef_summary_TRL.csv
#   - coef_summary_LPA.csv
#   - ridge_boot_LDL_full.rds
#   - ridge_boot_TRL_full.rds
#   - ridge_boot_LPA_full.rds
#   - ridge_scores.rda
#   - ldl_protein_data.rda
#   - trl_protein_data.rda
#   - lpa_protein_data.rda
#
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
})


# ==============================================================================
# PROTEIN SELECTION AND Z-SCORE TRANSFORMATION
# ==============================================================================

# Select fraction-specific proteins
protein_ldl <- combined_long_signif %>%
  filter(fraction == "LDL") %>%
  pull(gene)

protein_trl <- combined_long_signif %>%
  filter(fraction == "TRL") %>%
  pull(gene)

protein_lpa <- combined_long_signif %>%
  filter(fraction == "Lp(a)") %>%
  pull(gene)

# Extract protein measurements
protein_data <- dataset %>%
  dplyr::select(eid, A1BG:ZPR1)

# Z-score transformation
zscore_protein <- protein_data %>%
  mutate(
    across(
      -eid,
      ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
      .names = "{.col}_rz"
    )
  ) %>%
  dplyr::select(eid, ends_with("_rz"))

# Save z-scored proteins
save(zscore_protein, file = "zscore_protein.rda")

# Create z-scored protein name vectors
protein_ldl_rz <- paste0(protein_ldl, "_rz")
protein_trl_rz <- paste0(protein_trl, "_rz")
protein_lpa_rz <- paste0(protein_lpa, "_rz")

# Extract fraction-specific z-scored protein data
protein_ldl_data <- zscore_protein %>%
  dplyr::select(all_of(protein_ldl_rz))

protein_trl_data <- zscore_protein %>%
  dplyr::select(all_of(protein_trl_rz))

protein_lpa_data <- zscore_protein %>%
  dplyr::select(all_of(protein_lpa_rz))

# Save selected z-score protein datasets
save(protein_ldl_data, file = "ldl_protein_data.rda")

save(protein_trl_data, file = "trl_protein_data.rda")

save(protein_lpa_data, file = "lpa_protein_data.rda")

# ==============================================================================
# RESIDUALIZED LIPOPROTEIN OUTCOMES
# ==============================================================================

# Regress each lipoprotein fraction on the other two to obtain residuals
model_trlp <- lm(conven_trlp ~ conven_ldlp + conven_lpa_corr, data = dataset)
model_ldlp <- lm(conven_ldlp ~ conven_lpa_corr + conven_trlp, data = dataset)
model_lpa_corr <- lm(conven_lpa_corr ~ conven_ldlp + conven_trlp, data = dataset)

# Add residuals to dataset
dataset <- dataset %>%
  mutate(
    resid_trlp = model_trlp$residuals,
    resid_ldlp = model_ldlp$residuals,
    resid_lpa_corr = model_lpa_corr$residuals
  )

# Extract residualized outcomes for modeling
resid_ldl_data <- dataset %>% select(resid_ldlp)
resid_trl_data <- dataset %>% select(resid_trlp)
resid_lpa_data <- dataset %>% select(resid_lpa_corr)


# ==============================================================================
# RIDGE REGRESSION WITH BOOTSTRAP
# ==============================================================================

# ---------- Utility Functions ----------

# Create deterministic cross-validation fold assignments
.make_foldid <- function(n, nfolds, seed) {
  set.seed(seed)
  foldid <- rep(seq_len(nfolds), length.out = n)
  sample(foldid, n)  # shuffle assignment
}

# ---------- Core Ridge Regression Runner ----------
#
# Arguments:
#   X_df           - Data frame of proteins (already Z-scored)
#   y_df           - Data frame with one column (covariate-adjusted outcome)
#   frac_tag       - Fraction label ("LDL" / "TRL" / "Lp(a)") for file naming
#   B              - Number of bootstrap resamples
#   nfolds         - Inner cross-validation folds
#   lambda_rule    - "1se" (stable) or "min" (best-CV)
#   base_seed      - Random seed for reproducibility
#   save_boot_path - Path to save per-bootstrap coefficients (optional)
#   save_every     - Checkpointing frequency
#
# Returns:
#   List containing:
#     - summary_df: Data frame with coefficient summaries
#     - full_lambda: Lambda selected from full-data CV
#     - full_intercept: Intercept from full-data model
#     - full_coef_all: All coefficients (including intercept)
#     - boot_coef_matrix: Bootstrap coefficient matrix (proteins x B)
#     - boot_lambdas: Lambda values from each bootstrap sample
#
run_ridge_with_boot <- function(
    X_df,
    y_df,
    frac_tag,
    B = 500,
    nfolds = 10,
    lambda_rule = c("1se", "min"),
    base_seed = 20250924L,
    save_boot_path = NULL,
    save_every = 50
) {
  lambda_rule <- match.arg(lambda_rule)
  
  # Build matrices
  X <- makeX(X_df)
  y_mat <- makeX(y_df)
  
  if (is.data.frame(y_mat)) {
    y_mat <- as.matrix(y_mat)
  }
  
  if (NCOL(y_mat) != 1L) {
    stop("y must have exactly one column after makeX().")
  }
  
  y <- as.numeric(y_mat[, 1])
  
  if (NROW(X) != length(y)) {
    stop("X and y have different numbers of rows.")
  }
  
  p <- ncol(X)
  protein_names <- colnames(X)
  
  if (is.null(protein_names)) {
    protein_names <- paste0("feat_", seq_len(p))
  }
  
  # Full-data inner CV to select lambda
  set.seed(base_seed)
  foldid_full <- .make_foldid(n = nrow(X), nfolds = nfolds, seed = base_seed)
  
  cv_full <- cv.glmnet(
    x = X,
    y = y,
    alpha = 0,
    family = "gaussian",
    nfolds = nfolds,
    foldid = foldid_full,
    intercept = TRUE,
    standardize = FALSE,
    parallel = FALSE
  )
  
  lambda_star_full <- if (lambda_rule == "1se") {
    cv_full$lambda.1se
  } else {
    cv_full$lambda.min
  }
  
  # Fit ridge at selected lambda on full data
  fit_full <- glmnet(
    x = X,
    y = y,
    alpha = 0,
    family = "gaussian",
    lambda = lambda_star_full,
    intercept = TRUE,
    standardize = FALSE
  )
  
  coef_full <- as.numeric(coef(fit_full))
  names(coef_full) <- c("(Intercept)", protein_names)
  coef_full_noint <- coef_full[-1]
  
  # Bootstrap with inner CV in each resample
  boot_coef_mat <- matrix(
    NA_real_,
    nrow = p,
    ncol = B,
    dimnames = list(protein_names, paste0("b", seq_len(B)))
  )
  lambdas_boot <- numeric(B)
  
  # Initialize checkpointing file if requested
  if (!is.null(save_boot_path)) {
    saveRDS(
      list(
        partial = TRUE,
        coef = boot_coef_mat[, 0, drop = FALSE],
        lambdas = lambdas_boot[0]
      ),
      file = save_boot_path
    )
  }
  
  # Bootstrap loop
  for (b in seq_len(B)) {
    # Resample with replacement
    set.seed(base_seed + b)
    idx <- sample.int(n = nrow(X), size = nrow(X), replace = TRUE)
    
    # Inner CV on bootstrap sample
    foldid_b <- .make_foldid(
      n = length(idx),
      nfolds = nfolds,
      seed = base_seed + 10L * b
    )
    
    cv_b <- cv.glmnet(
      x = X[idx, , drop = FALSE],
      y = y[idx],
      alpha = 0,
      family = "gaussian",
      nfolds = nfolds,
      foldid = foldid_b,
      intercept = TRUE,
      standardize = FALSE,
      parallel = FALSE
    )
    
    lambda_star_b <- if (lambda_rule == "1se") {
      cv_b$lambda.1se
    } else {
      cv_b$lambda.min
    }
    lambdas_boot[b] <- lambda_star_b
    
    # Fit at lambda* on bootstrap sample
    fit_b <- glmnet(
      x = X[idx, , drop = FALSE],
      y = y[idx],
      alpha = 0,
      family = "gaussian",
      lambda = lambda_star_b,
      intercept = TRUE,
      standardize = FALSE
    )
    
    cf <- as.numeric(coef(fit_b))[-1]
    boot_coef_mat[, b] <- cf
    
    # Memory management
    rm(idx, foldid_b, cv_b, fit_b, cf)
    gc(verbose = FALSE)
    
    # Periodic checkpointing
    if (!is.null(save_boot_path) && (b %% save_every == 0 || b == B)) {
      saveRDS(
        list(
          partial = (b < B),
          coef = boot_coef_mat[, seq_len(b), drop = FALSE],
          lambdas = lambdas_boot[seq_len(b)]
        ),
        file = save_boot_path
      )
    }
  }
  
  # Summarize bootstrap distributions
  boot_median <- apply(boot_coef_mat, 1, stats::median, na.rm = TRUE)
  boot_var <- apply(boot_coef_mat, 1, stats::var, na.rm = TRUE)
  
  summary_df <- data.frame(
    protein = protein_names,
    coef_full_lambda = as.numeric(coef_full_noint),
    boot_median = as.numeric(boot_median),
    boot_variance = as.numeric(boot_var),
    stringsAsFactors = FALSE,
    row.names = protein_names
  )
  
  list(
    summary_df = summary_df,
    full_lambda = lambda_star_full,
    full_intercept = unname(coef_full[1]),
    full_coef_all = coef_full,
    boot_coef_matrix = boot_coef_mat,
    boot_lambdas = lambdas_boot
  )
}


# ==============================================================================
# RUN RIDGE REGRESSION FOR ALL THREE FRACTIONS
# ==============================================================================

# Set parameters
B <- 1             # Number of bootstrap resamples (increase for production)
nfolds <- 10       # Cross-validation folds
lambda_rule <- "1se"  # Use lambda.1se for stability (or "min" for lambda.min)

# LDL fraction
res_ldl <- run_ridge_with_boot(
  X_df = protein_ldl_data,
  y_df = resid_ldl_data,
  frac_tag = "LDL",
  B = B,
  nfolds = nfolds,
  lambda_rule = lambda_rule,
  base_seed = 123L,
  save_boot_path = "boot_coefs_LDL.rds",
  save_every = 50
)

# TRL fraction
res_trl <- run_ridge_with_boot(
  X_df = protein_trl_data,
  y_df = resid_trl_data,
  frac_tag = "TRL",
  B = B,
  nfolds = nfolds,
  lambda_rule = lambda_rule,
  base_seed = 123L,
  save_boot_path = "boot_coefs_TRL.rds",
  save_every = 50
)

# Lp(a) fraction
res_lpa <- run_ridge_with_boot(
  X_df = protein_lpa_data,
  y_df = resid_lpa_data,
  frac_tag = "LPA",
  B = B,
  nfolds = nfolds,
  lambda_rule = lambda_rule,
  base_seed = 123L,
  save_boot_path = "boot_coefs_LPA.rds",
  save_every = 50
)


# ==============================================================================
# EXTRACT AND SAVE COEFFICIENT TABLES
# ==============================================================================

# Extract coefficient summary tables
coef_table_LDL <- res_ldl$summary_df
coef_table_TRL <- res_trl$summary_df
coef_table_LPA <- res_lpa$summary_df

# Save coefficient tables as CSV
write.csv(coef_table_LDL, file = "coef_summary_LDL.csv", row.names = FALSE)
write.csv(coef_table_TRL, file = "coef_summary_TRL.csv", row.names = FALSE)
write.csv(coef_table_LPA, file = "coef_summary_LPA.csv", row.names = FALSE)

# Save full result objects (contains all bootstrapped coefficients)
saveRDS(res_ldl, file = "ridge_boot_LDL_full.rds")
saveRDS(res_trl, file = "ridge_boot_TRL_full.rds")
saveRDS(res_lpa, file = "ridge_boot_LPA_full.rds")

# Quick preview
head(coef_table_LDL)
head(coef_table_TRL)
head(coef_table_LPA)


# ==============================================================================
# SCORE CALCULATION HELPER FUNCTIONS
# ==============================================================================

# Check for duplicate values
stop_if_dups <- function(x, what) {
  d <- duplicated(x)
  if (any(d, na.rm = TRUE)) {
    dup_vals <- unique(x[d])
    stop(sprintf(
      "Duplicate %s detected: %s",
      what,
      paste(head(dup_vals, 10), collapse = ", ")
    ))
  }
}

# Extract participant IDs from protein data frame
extract_ids <- function(X_df, dataset, df_name) {
  # Priority: explicit eid column > match with dataset$eid > rownames
  if ("eid" %in% names(X_df)) {
    id <- X_df$eid
  } else if ("eid" %in% names(dataset) && nrow(X_df) == nrow(dataset)) {
    id <- dataset$eid
  } else if (!is.null(rownames(X_df))) {
    id <- rownames(X_df)
  } else {
    stop(sprintf(
      "Cannot determine participant IDs for %s. Add an 'eid' column or ensure row counts match 'dataset'.",
      df_name
    ))
  }
  
  if (any(is.na(id))) {
    stop(sprintf("Missing eid values detected for %s.", df_name))
  }
  
  stop_if_dups(id, sprintf("ids in %s", df_name))
  
  return(id)
}

# Compute weighted protein scores
align_and_score <- function(X_df, coef_table, prefix, df_name) {
  # Validate coefficient table structure
  required_cols <- c("protein", "coef_full_lambda", "boot_median")
  if (!all(required_cols %in% names(coef_table))) {
    stop("coef_table must have columns: protein, coef_full_lambda, boot_median")
  }
  
  # Check for duplicates
  stop_if_dups(colnames(X_df), sprintf("column names in %s", df_name))
  stop_if_dups(coef_table$protein, sprintf("proteins in coef_table for %s", prefix))
  
  # Find common proteins between data and coefficient table
  common <- intersect(colnames(X_df), coef_table$protein)
  
  if (length(common) == 0) {
    stop(sprintf(
      "[%s] No overlapping protein names between data and coefficient table.",
      prefix
    ))
  }
  
  # Build aligned matrices and coefficient vectors
  X <- as.matrix(X_df[, common, drop = FALSE])
  storage.mode(X) <- "double"
  
  w_full <- setNames(coef_table$coef_full_lambda, coef_table$protein)[common]
  w_boot <- setNames(coef_table$boot_median, coef_table$protein)[common]
  
  # Calculate scores as weighted sums (no intercept)
  score_full <- as.vector(X %*% w_full)
  score_boot <- as.vector(X %*% w_boot)
  
  # Return as data frame with descriptive column names
  out <- data.frame(
    setNames(
      list(score_full, score_boot),
      c(
        paste0(prefix, "_ridge_full"),
        paste0(prefix, "_ridge_boot")
      )
    ),
    check.names = FALSE
  )
  
  return(out)
}

# Merge multiple score data frames by eid
merge_by_eid <- function(dataset, ...) {
  # Start with dataset$eid
  out <- data.frame(eid = dataset$eid, check.names = FALSE)
  
  # Merge each additional data frame
  dots <- list(...)
  for (df in dots) {
    out <- merge(out, df, by = "eid", all = TRUE, sort = FALSE)
  }
  
  return(out)
}


# ==============================================================================
# COMPUTE PROTEIN SCORES
# ==============================================================================

# Extract participant IDs for each fraction
eid_ldl <- extract_ids(protein_ldl_data, dataset, "protein_ldl_data")
eid_trl <- extract_ids(protein_trl_data, dataset, "protein_trl_data")
eid_lpa <- extract_ids(protein_lpa_data, dataset, "protein_lpa_data")

# Compute scores for each fraction (full-data and bootstrap-median weights)
scores_ldl_only <- align_and_score(
  protein_ldl_data,
  coef_table_LDL,
  "LDL",
  "protein_ldl_data"
)

scores_trl_only <- align_and_score(
  protein_trl_data,
  coef_table_TRL,
  "TRL",
  "protein_trl_data"
)

scores_lpa_only <- align_and_score(
  protein_lpa_data,
  coef_table_LPA,
  "LPA",
  "protein_lpa_data"
)

# Add eid column to each score data frame
scores_ldl <- data.frame(
  eid = eid_ldl,
  scores_ldl_only,
  check.names = FALSE
)

scores_trl <- data.frame(
  eid = eid_trl,
  scores_trl_only,
  check.names = FALSE
)

scores_lpa <- data.frame(
  eid = eid_lpa,
  scores_lpa_only,
  check.names = FALSE
)


# ==============================================================================
# MERGE SCORES INTO FINAL DATASET
# ==============================================================================

# Merge all scores with dataset by eid (full outer join)
final_scores <- merge_by_eid(dataset, scores_ldl, scores_trl, scores_lpa)

# Save final scores
save(final_scores, file = "ridge_scores_092625.rda")

# Quick summary
str(final_scores)
summary(final_scores[, -1])  # Summaries of the six score columns
