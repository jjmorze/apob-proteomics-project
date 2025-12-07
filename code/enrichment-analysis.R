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
#   Required packages: enrichR, GO.db, AnnotationDbi, dplyr
#
# Goal:
#   Perform GO enrichment analysis (Enrichr) with GO leaf-term filtering
#   for protein biomarkers across three lipoprotein fractions (LDL, TRL, Lp(a)).
#   The analysis:
#   1) Queries Enrichr database for GO Biological Process enrichment
#   2) Filters results to leaf terms (no enriched descendants)
#   3) Generates fraction-specific and combined enrichment tables
#
# Required input data:
#   - combined_long_signif.rda
#
# Outputs:
#   - ea_ldl.rda
#   - ea_trl.rda
#   - ea_lpa.rda
#   - merged_full.rda
#   - table-pathway-enrichment.csv: export-ready enrichment table
#
# Notes:
#   - Requires internet access for Enrichr API queries
#   - GO leaf filtering removes parent terms with enriched children
#   - Database: GO_Biological_Process_2023
#
# ==============================================================================


# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(enrichR)
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
})


# ==============================================================================
# FUNCTION
# ==============================================================================

# Run Enrichr with GO Leaf Term Filtering
#
# Performs GO enrichment via Enrichr and filters to leaf terms only
# (terms with no enriched descendants in the ontology).
#
# @param gene_symbols Character vector of gene symbols
# @param database Character, Enrichr database name
# @return List with $full (all enriched terms) and $leaf (filtered terms)
get_enrichr_with_leaf_filtering <- function(gene_symbols,
                                            database = "GO_Biological_Process_2023") {
  
  # Input validation
  gene_symbols <- unique(na.omit(gene_symbols))
  if (length(gene_symbols) == 0L) {
    warning("No valid gene symbols supplied.")
    return(list(full = NULL, leaf = NULL))
  }
  
  # Step 1: Run Enrichr
  enriched_full <- tryCatch(
    {
      res <- enrichr(gene_symbols, database)
      res[[database]]
    },
    error = function(e) {
      warning("enrichr() failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(enriched_full) || nrow(enriched_full) == 0) {
    warning("No enrichment results found.")
    return(list(full = NULL, leaf = NULL))
  }
  
  # Step 2: Extract GO IDs from term names (format: "... (GO:#######)")
  enriched_full <- enriched_full %>%
    mutate(
      GOID = sub(".*\\((GO:\\d{7})\\)", "\\1", Term, perl = TRUE)
    ) %>%
    filter(grepl("^GO:\\d{7}$", GOID))
  
  if (nrow(enriched_full) == 0) {
    warning("No GO IDs could be parsed from Term.")
    return(list(full = NULL, leaf = NULL))
  }
  
  # Step 3: Load GO BP offspring DAG
  go_offspring_env <- get("GOBPOFFSPRING", envir = asNamespace("GO.db"))
  
  # Step 4: Get descendants of each enriched term
  enriched_goids <- enriched_full$GOID
  go_descendants <- AnnotationDbi::mget(
    enriched_goids,
    envir = go_offspring_env,
    ifnotfound = NA
  )
  
  # Step 5: Identify terms with enriched descendants
  has_enriched_descendants <- vapply(go_descendants, function(desc) {
    if (all(is.na(desc))) return(FALSE)
    any(desc %in% enriched_goids)
  }, logical(1))
  
  # Step 6: Keep only leaf GO terms (no enriched descendants)
  leaf_go_ids <- enriched_goids[!has_enriched_descendants]
  enriched_leaf <- enriched_full %>% 
    filter(GOID %in% leaf_go_ids)
  
  # Return both full and leaf results
  list(full = enriched_full, leaf = enriched_leaf)
}


# Drop Legacy P-value Columns
#
# @param df Data frame
# @return Data frame with Old.P.value and Old.Adjusted.P.value removed
drop_old_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    dplyr::select(-dplyr::any_of(c("Old.P.value", "Old.Adjusted.P.value")))
}


#' Add Fraction Label to Enrichment Results
#'
# @param df Data frame of enrichment results
# @param frac_label Character, fraction name
# @return Data frame with added 'fraction' column
bind_with_fraction <- function(df, frac_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% mutate(fraction = frac_label)
}


# ==============================================================================
# PREPARE GENE LISTS BY FRACTION
# ==============================================================================

protein_ldl <- combined_long_signif %>%
  filter(fraction == "LDL") %>%
  pull(gene)

protein_trl <- combined_long_signif %>%
  filter(fraction == "TRL") %>%
  pull(gene)

protein_lpa <- combined_long_signif %>%
  filter(fraction == "Lp(a)") %>%
  pull(gene)


# ==============================================================================
# RUN ENRICHMENT ANALYSIS PER FRACTION
# ==============================================================================

ea_ldl <- get_enrichr_with_leaf_filtering(
  protein_ldl, 
  database = "GO_Biological_Process_2023"
)

ea_trl <- get_enrichr_with_leaf_filtering(
  protein_trl, 
  database = "GO_Biological_Process_2023"
)

ea_lpa <- get_enrichr_with_leaf_filtering(
  protein_lpa, 
  database = "GO_Biological_Process_2023"
)

# Save fraction-specific results in a raw form
save(
  ea_ldl, 
  file = "PATH/ea_ldl.rda"
)

save(
  ea_trl, 
  file = "PATH/ea_trl.rda"
)

save(
  ea_lpa, 
  file = "PATH/ea_lpa.rda"
)

# ==============================================================================
# MERGE RESULTS ACROSS FRACTIONS
# ==============================================================================

# Combine full results
merged_full <- bind_rows(
  bind_with_fraction(drop_old_cols(ea_ldl$full), "LDL"),
  bind_with_fraction(drop_old_cols(ea_trl$full), "TRL"),
  bind_with_fraction(drop_old_cols(ea_lpa$full), "Lp(a)")
)

# Combine leaf-only results
merged_leaf <- bind_rows(
  bind_with_fraction(drop_old_cols(ea_ldl$leaf), "LDL"),
  bind_with_fraction(drop_old_cols(ea_trl$leaf), "TRL"),
  bind_with_fraction(drop_old_cols(ea_lpa$leaf), "Lp(a)")
)

# Sort by fraction and adjusted p-value
if (!is.null(merged_full) && nrow(merged_full) > 0 && 
    "Adjusted.P.value" %in% names(merged_full)) {
  merged_full <- merged_full %>% 
    arrange(fraction, Adjusted.P.value)
}

if (!is.null(merged_leaf) && nrow(merged_leaf) > 0 && 
    "Adjusted.P.value" %in% names(merged_leaf)) {
  merged_leaf <- merged_leaf %>% 
    arrange(fraction, Adjusted.P.value)
}

# Save raw merged table with enrichment analysis
save(
  merged_full, 
  file = "PATH/merged_full.rda"
)

# ==============================================================================
# EXPORT TABLE
# ==============================================================================

# Create export-ready table with ordered columns and clean names
# Column order: Fraction, Term, Overlap, Combined Score, P, P adjusted, Genes
merged_full_output <- merged_full %>%
  # Enforce fraction factor ordering
  mutate(fraction = factor(fraction, levels = c("LDL", "TRL", "Lp(a)"))) %>%
  arrange(fraction, Adjusted.P.value) %>%
  # Select and order columns (drop Odds.Ratio and GOID)
  dplyr::select(
    fraction, 
    Term, 
    Overlap, 
    Combined.Score, 
    P.value, 
    Adjusted.P.value, 
    Genes
  ) %>%
  # Rename for clean export
  dplyr::rename(
    Fraction = fraction,
    `Combined Score` = Combined.Score,
    P = P.value,
    `P adjusted` = Adjusted.P.value
  )

# Write to CSV
write_csv(merged_full_output, file = "table-pathway-enrichment.csv")