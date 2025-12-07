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
#   Required packages: dplyr, purrr, stringr, survival, broom, tibble
#
# Goal:
#   Cox proportional hazards models for lipoprotein-associated protein cluster
#   PC1 scores. The script:
#   1) Filters clusters with >1 protein
#   2) Builds modeling datasets with cluster scores
#   3) Computes correlations between cluster scores and conventional fractions
#   4) Fits Cox models (CAD2 outcome) adjusting for demographics and covariates
#   5) Generates summary tables with HR estimates and correlations
#
# Required input data:
#   - dataset              : main phenotypic + proteomic dataset
#   - rev_scores_trl       : TRL cluster PC1 scores
#   - rev_scores_lpa       : Lp(a) cluster PC1 scores
#   - membership_trl       : TRL cluster membership (protein, cluster)
#   - membership_lpa       : Lp(a) cluster membership (protein, cluster)
#
# Intermediate outputs:
#   - df_model_trl         : TRL modeling dataset
#   - df_model_lpa         : Lp(a) modeling dataset
#
# Main outputs:
#   - cox_clusters_long                        : combined results (all fractions)
#   - table-ukbb-coxph-network-cluster.csv     : formatted table
#   - df_trl_hima.rda                          : TRL dataset
#   - df_lpa_hima.rda                          : Lp(a) dataset
#
# ==============================================================================


# ==============================================================================
#  LOAD PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(survival)
  library(broom)
  library(tibble)
})


# ==============================================================================
#  CONFIGURATION
# ==============================================================================

# --- Outcome and covariates ---
covars_num_cox       <- c("age", "bmi_imp", "conven_hba1c_imp", "sbp_mean_imp")
covars_categorical   <- c("sex", "race_imp", "smoking_status_imp", "med_lipid")

# --- Conventional fraction variables (for correlations only) ---
conven_cols <- c("conven_trlp", "conven_ldlp", "conven_lpa_corr")


# ==============================================================================
#  PREPARE BASE MODELING DATASET
# ==============================================================================

# Select required columns and scale numeric covariates
base_needed <- c(
  "cad2_ukb_survtime", 
  "cad2_ukb_ind",
  covars_num_cox, 
  covars_categorical, 
  conven_cols
)

df_model_base <- dataset %>%
  dplyr::select(all_of(intersect(base_needed, colnames(dataset)))) %>%
  dplyr::mutate(
    dplyr::across(
      all_of(intersect(covars_num_cox, colnames(.))), 
      ~ as.numeric(scale(.))
    )
  )


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# Filter clusters with >1 protein and prepare metadata
# ------------------------------------------------------------------------------
prep_cluster_meta <- function(membership_tbl) {
  # Args:
  #   membership_tbl: tibble with columns (protein, cluster)
  # Returns:
  #   tibble with cluster, n_prot, proteins (comma-separated)
  
  membership_tbl %>%
    mutate(protein_clean = sub("_rz$", "", protein)) %>%
    group_by(cluster) %>%
    summarise(
      n_prot   = n(),
      proteins = paste(sort(protein_clean), collapse = ", "),
      .groups  = "drop"
    ) %>%
    filter(n_prot > 1) %>%
    arrange(cluster)
}


# ------------------------------------------------------------------------------
# Build per-fraction modeling dataset
# ------------------------------------------------------------------------------
build_fraction_df <- function(scores_df, cluster_meta) {
  # Args:
  #   scores_df: data frame with cluster PC1 scores (columns: C1_PC1, C2_PC1, ...)
  #   cluster_meta: tibble from prep_cluster_meta()
  # Returns:
  #   list with df (merged dataset) and keep_cols (cluster score column names)
  
  # Identify cluster score columns for retained clusters
  keep_clusters <- cluster_meta$cluster
  keep_cols <- paste0("C", keep_clusters, "_PC1")
  keep_cols <- intersect(keep_cols, colnames(scores_df))
  
  if (length(keep_cols) == 0) {
    warning("No cluster score columns found after filtering.")
    return(list(df = df_model_base, keep_cols = character(0)))
  }
  
  # Merge base dataset with cluster scores
  df_mod <- dplyr::bind_cols(
    df_model_base,
    scores_df %>% dplyr::select(all_of(keep_cols))
  )
  
  list(df = df_mod, keep_cols = keep_cols)
}


# ------------------------------------------------------------------------------
# Compute correlations with conventional fraction
# ------------------------------------------------------------------------------
corr_vs_conven <- function(df_mod, cluster_cols, conven_col) {
  # Args:
  #   df_mod: modeling data frame
  #   cluster_cols: character vector of cluster score column names
  #   conven_col: conventional fraction variable name
  # Returns:
  #   tibble with ClusterScore and cor_with_conven
  
  if (!(conven_col %in% colnames(df_mod))) {
    warning(sprintf("Conventional column '%s' not found; correlations set to NA.", 
                    conven_col))
    return(tibble(ClusterScore = cluster_cols, cor_with_conven = NA_real_))
  }
  
  purrr::map_dfr(cluster_cols, function(cc) {
    if (!(cc %in% colnames(df_mod))) {
      tibble(ClusterScore = cc, cor_with_conven = NA_real_)
    } else {
      v1 <- df_mod[[cc]]
      v2 <- df_mod[[conven_col]]
      cor_val <- suppressWarnings(cor(v1, v2, use = "pairwise.complete.obs"))
      tibble(ClusterScore = cc, cor_with_conven = unname(cor_val))
    }
  })
}


# ------------------------------------------------------------------------------
# Fit Cox model for single cluster score
# ------------------------------------------------------------------------------
fit_cox_cluster <- function(df_mod, cluster_col) {
  # Args:
  #   df_mod: modeling data frame
  #   cluster_col: cluster score column name
  # Returns:
  #   tibble with ClusterScore, estimate (log HR), std.error, p.value
  
  # Build formula
  rhs <- paste(
    c(
      sprintf("scale(`%s`)", cluster_col),
      intersect(covars_num_cox, colnames(df_mod)),
      intersect(covars_categorical, colnames(df_mod))
    ),
    collapse = " + "
  )
  fml <- as.formula(paste0("Surv(cad2_ukb_survtime, cad2_ukb_ind) ~ ", rhs))
  
  # Fit model
  mdl <- survival::coxph(fml, data = df_mod)
  
  # Extract results for cluster score term
  term_regex <- paste0(
    "^scale\\((`)?",
    stringr::str_replace_all(cluster_col, "([\\W])", "\\\\\\1"),
    "(\\1)?\\)$"
  )
  
  td <- broom::tidy(mdl) %>%
    dplyr::filter(stringr::str_detect(term, term_regex)) %>%
    dplyr::transmute(
      ClusterScore = cluster_col,
      estimate     = estimate,   # log(HR) per 1 SD of cluster score
      std.error    = std.error,
      p.value      = p.value
    )
  
  # Handle missing results
  if (nrow(td) == 0) {
    td <- tibble(
      ClusterScore = cluster_col,
      estimate     = NA_real_,
      std.error    = NA_real_,
      p.value      = NA_real_
    )
  }
  
  td
}


# ------------------------------------------------------------------------------
# Complete pipeline for one lipoprotein fraction
# ------------------------------------------------------------------------------
run_fraction <- function(frac_label, scores_df, membership_tbl, conven_var) {
  # Args:
  #   frac_label: character label for fraction (e.g., "TRL", "Lp(a)")
  #   scores_df: cluster PC1 scores
  #   membership_tbl: cluster membership table
  #   conven_var: conventional fraction variable name
  # Returns:
  #   list with df_model and results tibble
  
  # Step 1: Prepare cluster metadata (filter to >1 protein)
  meta <- prep_cluster_meta(membership_tbl)
  
  # Step 2: Build modeling dataset
  built <- build_fraction_df(scores_df, meta)
  df_mod <- built$df
  present_cols <- built$keep_cols
  
  if (length(present_cols) == 0) {
    return(list(
      df_model = df_mod,
      results  = tibble(
        Fraction            = character(0),
        Cluster             = integer(0),
        proteins            = character(0),
        logHR               = numeric(0),
        SE                  = numeric(0),
        P                   = numeric(0),
        cor_with_conven     = numeric(0)
      )
    ))
  }
  
  # Step 3: Compute correlations with conventional fraction
  cors <- corr_vs_conven(df_mod, present_cols, conven_var)
  
  # Step 4: Fit Cox models for each cluster
  res <- purrr::map_dfr(present_cols, ~ {
    tryCatch(
      fit_cox_cluster(df_mod, .x),
      error = function(e) {
        tibble(
          ClusterScore = .x,
          estimate     = NA_real_,
          std.error    = NA_real_,
          p.value      = NA_real_
        )
      }
    )
  })
  
  # Step 5: Assemble final results
  res <- res %>%
    left_join(cors, by = "ClusterScore") %>%
    mutate(
      Cluster = as.integer(stringr::str_match(ClusterScore, "^C(\\d+)_PC1$")[, 2])
    ) %>%
    select(Cluster, logHR = estimate, SE = std.error, P = p.value, cor_with_conven) %>%
    left_join(meta %>% select(Cluster = cluster, proteins), by = "Cluster") %>%
    mutate(Fraction = frac_label, .before = 1) %>%
    arrange(Cluster)
  
  list(
    df_model = df_mod,
    results  = res
  )
}


# ==============================================================================
#  RUN ANALYSES
# ==============================================================================

# --- TRL fraction ---
out_trl <- run_fraction("TRL", rev_scores_trl, membership_trl, "conven_trlp")

# --- Lp(a) fraction ---
out_lpa <- run_fraction("Lp(a)", rev_scores_lpa, membership_lpa, "conven_lpa_corr")


# ==============================================================================
#  EXTRACT MODELING DATASETS
# ==============================================================================

df_model_trl <- out_trl$df_model
df_model_lpa <- out_lpa$df_model


# ==============================================================================
#  COMBINE AND FORMAT RESULTS
# ==============================================================================

# --- Combine all fraction results ---
cox_clusters_long <- bind_rows(
  out_trl$results, 
  out_lpa$results
) %>%
  select(Fraction, Cluster, proteins, logHR, SE, P, cor_with_conven)

# --- Create formatted output table ---
cox_clusters_long_output <- cox_clusters_long %>%
  rename(
    Proteins                                = proteins,
    `Correlation with lipoprotein fraction` = cor_with_conven
  ) %>%
  select(
    Fraction,
    Cluster,
    Proteins,
    `Correlation with lipoprotein fraction`,
    logHR,
    SE,
    P
  )


# ==============================================================================
#  SAVE OUTPUTS
# ==============================================================================

# --- Save formatted table ---
write_csv(
  cox_clusters_long_output, 
  file = "table-ukbb-coxph-network-cluster.csv"
)

# --- Save modeling datasets ---
save(df_model_trl, file = "df_trl_hima.rda")
save(df_model_lpa, file = "df_lpa_hima.rda")

# --- Display results preview ---
print(dplyr::glimpse(cox_clusters_long))
