# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Jakub Morze, MD, PhD
# Created: 2025-09-26
# Last modified: 2025-11-26
# Contact: jakub.morze@chalmers.se / jjmorze
# Version: 1.0.1
# License: MIT
#
# Dependencies:
#   R >= 4.0.0
#   Required packages: readr, readxl, tidyr, purrr, lubridate, tidyverse, missMDA, VIM, stringr
#
# Goal:
#   Build an analysis-ready UK Biobank dataset by:
#   1) deriving core tabularized variables
#   2) cleaning and imputing Olink proteomic data,
#   3) constructing CVD/metabolic endpoints from hospital, procedure,
#      mortality and self-report data,
#   4) deriving medication class indicators, and
#   5) merging all components, applying exclusions and covariate imputation,
#      and harmonising protein names to gene-based labels.
#
# Required input data:
#   - basic_covariates.csv # RESTRICTED
#   - olink_instance_0.csv # RESTRICTED
#   - icd10_registry.csv # RESTRICTED
#   - icd9_registry.csv # RESTRICTED
#   - opcs4_registry.csv # RESTRICTED
#   - mortality_registry.csv # RESTRICTED
#   - srd_registry.csv # RESTRICTED
#   - sro_registry.csv # RESTRICTED
#   - process_data.csv # RESTRICTED
#   - medications.csv # RESTRICTED
#   - ukbb-drug-coding.xlsx
#   - lpa_values_outofrange[rebadged_570811].csv # RESTRICTED
#   - name_uniprot_gene_pathway.xlsx
#
# Intermediate outputs:
#   - basic_covariates_output.rda
#   - olink_filtered.rda
#   - olink_filtered_imputed.rda
#   - cad_output.rda, 
#   - cad2_output.rda
#   - stroke_output.rda, 
#   - ischemic_stroke_output.rda
#   - pad_output.rda
#   - dm2_output.rda
#   - medications_output.rda
#
# Main outputs:
#   - dataset_side1.rda   (Olink APOB/LPA saved separately)
#   - dataset_side2.rda   (NMR data)
#   - dataset.rda         (final merged analysis dataset)
#
# ==============================================================================

# list of required packages
packages <- c(
  "readr", "readxl", "tidyr", "purrr", "lubridate",
  "tidyverse", "missMDA", "VIM", "stringr"
)

# install any missing packages, then load them
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    require(pkg, character.only = TRUE)
  }
}

# Libraries
library(readr)
library(readxl)
library(tidyr)
library(purrr)
library(lubridate)
library(tidyverse)
library(missMDA)
library(VIM)
library(stringr)


##### PREPARATION OF BASIC COVARIATES #####

basic_covariates <- read_csv("basic_covariates.csv")

basic_covariates <- basic_covariates %>% 
  rename(eid = `Participant ID`,
         age = `Age at recruitment`, 
         sex = Sex,
         ethnic_background = `Ethnic background | Instance 0`, 
         bmi = `Body mass index (BMI) | Instance 0`, 
         smoking_status = `Smoking status | Instance 0`, 
         conven_tchol = `Cholesterol | Instance 0`, 
         conven_ldl = `LDL direct | Instance 0`, 
         conven_tg = `Triglycerides | Instance 0`, 
         conven_hdl = `HDL cholesterol | Instance 0`, 
         conven_lpa = `Lipoprotein A | Instance 0`, 
         conven_apob = `Apolipoprotein B | Instance 0`,
         conven_hba1c = `Glycated haemoglobin (HbA1c) | Instance 0`, 
         conven_hscrp = `C-reactive protein | Instance 0`,
         conven_creat = `Creatinine | Instance 0`,
         blood_wbc = `White blood cell (leukocyte) count | Instance 0`, 
         blood_rdw = `Red blood cell (erythrocyte) distribution width | Instance 0`, 
         blood_lympho = `Lymphocyte count | Instance 0`, 
         blood_neutro = `Neutrophill count | Instance 0`,
         sbp_ar_a0 = `Systolic blood pressure, automated reading | Instance 0 | Array 0`,
         sbp_ar_a1 = `Systolic blood pressure, automated reading | Instance 0 | Array 1`, 
         sbp_mr_a0 = `Systolic blood pressure, manual reading | Instance 0 | Array 0`,
         sbp_mr_a1 = `Systolic blood pressure, manual reading | Instance 0 | Array 1`,
         dbp_ar_a0 = `Diastolic blood pressure, automated reading | Instance 0 | Array 0`, 
         dbp_ar_a1 = `Diastolic blood pressure, automated reading | Instance 0 | Array 1`,  
         dbp_mr_a0 = `Diastolic blood pressure, manual reading | Instance 0 | Array 0`, 
         dbp_mr_a1 = `Diastolic blood pressure, manual reading | Instance 0 | Array 1`
  )

# Create "race" variable
basic_covariates <- basic_covariates %>%
  mutate(
    race = case_when(
      ethnic_background %in% c("White", "British", "Irish", "Any other white background") ~ "White",
      ethnic_background %in% c("Black or Black British", "Caribbean", "African", "Any other Black background") ~ "Black",
      ethnic_background %in% c("Asian or Asian British", "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background") ~ "Asian",
      ethnic_background %in% c("Mixed", "White and Black Caribbean", "White and Black African", "White and Asian", 
                               "Any other mixed background", "Other ethnic group", "Do not know") ~ "Mixed/Other",
      ethnic_background %in% c("", "Prefer not to answer") ~ NA_character_,
      TRUE ~ "X"
    ),
    race = factor(race, levels = c("White", "Black", "Asian", "Mixed/Other"))
  )

# Recode "smoking_status" variable
basic_covariates <- basic_covariates %>%
  mutate(
    smoking_status = case_when(
      smoking_status %in% c("", "Prefer not to answer") ~ NA_character_,
      TRUE ~ smoking_status
    ),
    smoking_status = factor(smoking_status, 
                            levels = c("Never", "Previous", "Current"))
  )

# Create "sbp_mean" and "dbp_mean" variables
basic_covariates <- basic_covariates %>%
  mutate(
    sbp_auto_mean   = rowMeans(across(c(sbp_ar_a0, sbp_ar_a1)), na.rm = TRUE),
    sbp_manual_mean = rowMeans(across(c(sbp_mr_a0, sbp_mr_a1)), na.rm = TRUE),
    sbp_mean        = rowMeans(across(c(sbp_auto_mean, sbp_manual_mean)), na.rm = TRUE),
    
    dbp_auto_mean   = rowMeans(across(c(dbp_ar_a0, dbp_ar_a1)), na.rm = TRUE),
    dbp_manual_mean = rowMeans(across(c(dbp_mr_a0, dbp_mr_a1)), na.rm = TRUE),
    dbp_mean        = rowMeans(across(c(dbp_auto_mean, dbp_manual_mean)), na.rm = TRUE)
  )

# Create "egfr_mdrd" variable - requires geenrated "race". Creatinine is in umol/L
calculate_gfr_mdrd <- function(creat_umol, age, sex, race = NA_character_) {
  # Convert to mg/dL
  scr_mg_dl <- creat_umol / 88.4
  
  base <- 175 * (scr_mg_dl ^ -1.154) * (age ^ -0.203)
  
  female_factor <- if_else(sex == "Female", 0.742, 1, missing = 1)
  race_factor   <- if_else(race == "Black", 1.212, 1, missing = 1)
  
  out <- base * female_factor * race_factor
  
  if_else(is.na(creat_umol) | is.na(age) | is.na(sex), NA_real_, out)
}

basic_covariates <- basic_covariates %>%
  mutate(gfr_mdrd = calculate_gfr_mdrd(conven_creat, age, sex, race))

# Select final set of basic covariates
basic_covariates_output <- basic_covariates %>%
  select(
    eid,
    age,
    sex,
    bmi,
    race,
    smoking_status,
    conven_tchol,
    conven_ldl,
    conven_tg,
    conven_lpa,
    conven_hdl,
    conven_apob,
    conven_hba1c,
    conven_hscrp,
    conven_creat,
    gfr_mdrd,
    blood_wbc,
    blood_rdw,
    blood_lympho,
    blood_neutro,
    sbp_mean,
    dbp_mean
  )

# Remove temporary file
rm(basic_covariates)

# Save basic covariates file
save(basic_covariates_output, file = "basic_covariates_output.rda")

##### PROTEOMIC DATA PREPARATION #####

olink_instance_0 <- read_csv("olink_instance_0.csv")

# Count NAs per row #
na_count_per_row <- rowSums(is.na(olink_instance_0))

# Count NAs per column
na_count_per_column <- colSums(is.na(olink_instance_0))

# Define thresholds
threshold_column = 0.30 * nrow(olink_instance_0)  # 30% NAs for columns
threshold_row = 0.30 * ncol(olink_instance_0)     # 30% NAs for rows

# Select columns with < 30% NAs
valid_columns <- na_count_per_column < threshold_column
filtered_data_columns <- olink_instance_0[, valid_columns]

# Select rows with < 30% NAs
valid_rows <- na_count_per_row < threshold_row
filtered_data <- filtered_data_columns[valid_rows, ]

# Packages
library(tidyverse)
library(missMDA)

# 1) Drop ID and prep numeric data
data_to_impute <- filtered_data %>%
  select(-`Participant ID`) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

# 3) Impute with PCA
set.seed(123)
imp <- imputePCA(data_to_impute, ncp = 2, scale = TRUE)

# 4) Recombine with ID (ID is not used in imputation)
imputed_data <- filtered_data %>%
  select(`Participant ID`) %>%
  bind_cols(as_tibble(imp$completeObs))

# Rename file and ID variable
olink_filtered <- filtered_data %>% rename(eid = `Participant ID`)
olink_filtered_imputed <- imputed_data %>% rename(eid = `Participant ID`)

# Save imputed files 
save(olink_filtered, file = "olink_filtered.rda")
save(olink_filtered_imputed, file = "olink_filtered_imputed.rda")

##### HEALTH EVENT DATA PREPARATION #####

## Read raw extracted data ##

# Read ICD-10 hospitalisation data
read_icd10_registry <- function(file_path) {
  # Read the full dataset with extended guess_max
  df <- read_csv(file_path, guess_max = 502128, show_col_types = FALSE)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify all ICD10 date columns
  date_cols <- grep("Date of first in-patient diagnosis - ICD10", names(df), value = TRUE)
  
  # Convert date columns to Date class
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.Date(.x, format = "%Y-%m-%d")))
  
  df <- df %>% 
    mutate(ListD = str_split(`Diagnoses - ICD10`, "\\|"))
  
  return(df)
}

icd10_registry <- read_icd10_registry("icd10_registry.csv")

# Read ICD-9 hospitalisation data
read_icd9_registry <- function(file_path) {
  # Read the full dataset with extended guess_max
  df <- read_csv(file_path, guess_max = 502128, show_col_types = FALSE)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify all ICD9 date columns
  date_cols <- grep("Date of first in-patient diagnosis - ICD9", names(df), value = TRUE)
  
  # Convert date columns to Date class
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.Date(.x, format = "%Y-%m-%d")))
  
  df <- df %>% 
    mutate(ListD = str_split(`Diagnoses - ICD9`, "\\|"))
  
  return(df)
}

icd9_registry <- read_icd9_registry("icd9_registry.csv")

# Read OPCS-4 oerative procedures data
read_opcs4_registry <- function(file_path) {
  # Read the full dataset with extended guess_max
  df <- read_csv(file_path, guess_max = 502128, show_col_types = FALSE)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify all opcs4 date columns
  date_cols <- grep("Date of first operative procedure - OPCS4", names(df), value = TRUE)
  
  # Convert date columns to Date class
  df <- df %>%
    mutate(across(all_of(date_cols), ~ as.Date(.x, format = "%Y-%m-%d")))
  
  df <- df %>% 
    mutate(ListD = str_split(`Operative procedures - OPCS4`, "\\|"))
  
  return(df)
}

opcs4_registry <- read_opcs4_registry("opcs4_registry.csv")

# Read mortality data
read_mortality_registry <- function(file_path) {
  # Read full dataset with comma delimiter and extended guess_max
  df <- read_csv(file_path, guess_max = 502128, show_col_types = FALSE)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify the date column
  date_col <- "Date of death | Instance 0"
  
  # Identify all columns except eid and the date column
  cols_to_char <- setdiff(names(df), c("eid", date_col))
  
  # Coerce all non-id, non-date columns to character
  df <- df %>%
    mutate(
      across(all_of(cols_to_char), as.character),
      `{date_col}` := as.Date(.data[[date_col]], format = "%Y-%m-%d")
    )
  
  return(df)
}

mortality_registry <- read_mortality_registry("mortality_registry.csv")

# Read self-reported diseases
read_srd_registry <- function(file_path) {
  # Read column names first to customize col_types
  col_names <- names(read_csv(file_path, n_max = 0, show_col_types = FALSE))
  
  # Define types: default is character, except Participant ID which is numeric
  col_types <- rep("c", length(col_names))
  pid_index <- which(col_names == "Participant ID")
  if (length(pid_index) == 1) {
    col_types[pid_index] <- "d"  # numeric
  }
  
  # Read full dataset
  df <- read_csv(file_path, col_types = paste0(col_types, collapse = ""), show_col_types = FALSE, guess_max = 502128)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify all columns corresponding to year/age first occurred
  year_cols <- grep("year/age first occurred", names(df), value = TRUE)
  
  # Replace special responses with "1"
  df <- df %>%
    mutate(across(all_of(year_cols), ~ifelse(.x %in% c("Time uncertain/unknown", "Preferred not to answer"), "1", .x)))
  
  return(df)
}

srd_registry <- read_srd_registry("srd_registry.csv")

# Read self-reported operative procedures
read_sro_registry <- function(file_path) {
  # Read column names first to set types
  col_names <- names(read_csv(file_path, n_max = 0, show_col_types = FALSE))
  
  # Define column types: all character, except Participant ID
  col_types <- rep("c", length(col_names))
  pid_index <- which(col_names == "Participant ID")
  if (length(pid_index) == 1) {
    col_types[pid_index] <- "d"
  }
  
  # Read full dataset with original column names
  df <- read_csv(file_path, col_types = paste0(col_types, collapse = ""), show_col_types = FALSE, guess_max = 502128)
  
  # Rename Participant ID to eid
  df <- df %>% rename(eid = `Participant ID`)
  
  # Identify "Operation year/age first occurred" columns
  year_cols <- grep("Operation year/age first occurred", names(df), value = TRUE)
  
  # Replace "Time uncertain/unknown" and "Preferred not to answer" with "1"
  df <- df %>%
    mutate(across(all_of(year_cols), ~ifelse(.x %in% c("Time uncertain/unknown", "Preferred not to answer"), "1", .x)))
  
  return(df)
}

sro_registry <- read_sro_registry("sro_registry.csv")

# Read procedural data (study checkpoint dates)
read_process_data <- function(file_path) {
  # Read the full dataset
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Rename variables as specified
  df <- df %>%
    rename(
      eid         = `Participant ID`,
      age         = `Age at recruitment`,
      sex         = `Sex`,
      date_attend = `Date of attending assessment centre | Instance 0`,
      center      = `UK Biobank assessment centre | Instance 0`,
      date_lost   = `Date lost to follow-up`,
      reason_lost = `Reason lost to follow-up`,
      date_death  = `Date of death | Instance 0`
    )
  
  return(df)
}

process_data <- read_process_data("process_data.csv")

#### Load extraction functions ####

# Function - Extract ICD-10 dates from hospitalization registry
extract_icd10 <- function(digit4, data) {
  # Automatically extract 3-digit codes from 4-digit input
  digit3 <- sub("\\..*$", "", digit4)
  
  # Base data frame with all eids
  data_base <- data.frame(eid = data$eid)
  
  # List to store date columns
  date_columns_all <- list()
  
  for (i in seq_along(digit4)) {
    d3 <- digit3[i]
    d4 <- digit4[i]
    
    # Filter rows containing the current 3-digit code
    data1 <- data %>% filter(grepl(d3, `Diagnoses - ICD10`))
    data2 <- data1 %>% select(eid)
    
    # Select date columns and transpose into list of rows
    date_columns <- data1 %>% select(-eid, -`Diagnoses - ICD10`, -ListD)
    date_list <- transpose(as.list(date_columns))
    
    # Match current 4-digit code
    m_idx <- map(data1$ListD, ~pmatch(d4, .x))
    
    # Extract corresponding dates
    matched_dates <- map2(date_list, m_idx, function(x, idx) {
      if (is.na(idx) || is.null(idx) || length(idx) == 0) {
        return(NA_real_)
      } else {
        return(as.numeric(x[[idx]]))
      }
    }) %>%
      unlist() %>%
      as.Date(origin = "1970-01-01")
    
    # Merge matched dates back with full eid list
    col_full <- data.frame(eid = data$eid) %>%
      left_join(bind_cols(data2, tmp = matched_dates), by = "eid") %>%
      select(tmp)
    
    # Assign formatted column name
    colname <- paste0(d4, "_ukb_icd10_date")
    colnames(col_full) <- colname
    
    # Store in output list
    date_columns_all[[colname]] <- col_full
  }
  
  # Bind all into one data frame
  final_df <- bind_cols(data_base, date_columns_all)
  return(final_df)
}

# Function - Extract ICD-9 dates from hospitalization registry
extract_icd9 <- function(digit4, data) {
  # Extract 3-digit version (first 3 chars, no dots in ICD9)
  digit3 <- substr(digit4, 1, 3)
  
  # Base dataframe for joining later
  data_base <- data.frame(eid = data$eid)
  date_columns_all <- list()
  
  for (i in seq_along(digit4)) {
    d4 <- digit4[i]
    d3 <- digit3[i]
    
    # Filter participants whose ListD (ICD9 codes) contain the 3-digit code
    data_match <- data %>%
      filter(str_detect(`Diagnoses - ICD9`, fixed(d3))) %>%
      select(eid, `Diagnoses - ICD9`, ListD, starts_with("Date"))  # include only date arrays
    
    if (nrow(data_match) == 0) {
      # If no matches, fill with NAs
      date_columns_all[[paste0(d4, "_ukb_ICD9_date")]] <- rep(NA_Date_, nrow(data))
      next
    }
    
    # Transpose date columns to get row-wise lists
    date_lists <- transpose(as.list(select(data_match, starts_with("Date"))))
    
    # Match 4-digit code in each ListD and extract date
    matched_dates <- map2(data_match$ListD, date_lists, function(icd_list, date_row) {
      idx <- pmatch(d4, icd_list)
      if (is.na(idx) || idx < 1 || idx > length(date_row)) return(NA_Date_)
      as.Date(date_row[[idx]])
    })
    
    # Create output column
    tmp_df <- data.frame(eid = data_match$eid, tmp = unlist(matched_dates))
    col_full <- left_join(data_base, tmp_df, by = "eid") %>% pull(tmp)
    col_name <- paste0(d4, "_ukb_ICD9_date")
    date_columns_all[[col_name]] <- col_full
  }
  
  # Return merged data frame
  result <- bind_cols(data_base, as.data.frame(date_columns_all))
  return(result)
}

# Function - Extract OPCS-4 dates from operative procedure registry
extract_opcs4 <- function(digit4, data) {
  # Automatically extract 3-digit codes from 4-digit input
  digit3 <- sub("\\..*$", "", digit4)
  
  # Base data frame with all eids
  data_base <- data.frame(eid = data$eid)
  
  # List to store date columns
  date_columns_all <- list()
  
  for (i in seq_along(digit4)) {
    d3 <- digit3[i]
    d4 <- digit4[i]
    
    # Filter rows containing the current 3-digit code
    data1 <- data %>% filter(grepl(d3, `Operative procedures - OPCS4`))
    data2 <- data1 %>% select(eid)
    
    # Select date columns and transpose into list of rows
    date_columns <- data1 %>% select(-eid, -`Operative procedures - OPCS4`, -ListD)
    date_list <- transpose(as.list(date_columns))
    
    # Match current 4-digit code
    m_idx <- map(data1$ListD, ~pmatch(d4, .x))
    
    # Extract corresponding dates
    matched_dates <- map2(date_list, m_idx, function(x, idx) {
      if (is.na(idx) || is.null(idx) || length(idx) == 0) {
        return(NA_real_)
      } else {
        return(as.numeric(x[[idx]]))
      }
    }) %>%
      unlist() %>%
      as.Date(origin = "1970-01-01")
    
    # Merge matched dates back with full eid list
    col_full <- data.frame(eid = data$eid) %>%
      left_join(bind_cols(data2, tmp = matched_dates), by = "eid") %>%
      select(tmp)
    
    # Assign formatted column name
    colname <- paste0(d4, "_ukb_opcs4_date")
    colnames(col_full) <- colname
    
    # Store in output list
    date_columns_all[[colname]] <- col_full
  }
  
  # Bind all into one data frame
  final_df <- bind_cols(data_base, date_columns_all)
  return(final_df)
}

# Function - Extract ICD-10 from mortality registry
extract_death <- function(digit3, data) {
  # Collapse multiple ICD codes into a regex pattern
  pattern <- paste(digit3, collapse = "|")
  
  # Filter rows with a match in any of the ICD10 fields
  data1 <- data %>%
    filter(if_any(`Underlying (primary) cause of death: ICD10 | Instance 0`:`Contributory (secondary) causes of death: ICD10 | Instance 0 | Array 14`, 
                  ~grepl(pattern, .)))
  
  # Create indicator for primary cause
  data1 <- data1 %>%
    mutate(
      prim = if_else(if_any(contains("primary"), ~grepl(pattern, .)), 1, 0),
      cause = 1
    )
  
  # Select and rename
  data2 <- data1 %>%
    select(eid, `Date of death | Instance 0`, cause, prim) %>%
    rename(date_death = `Date of death | Instance 0`) %>%
    rename_with(~c("date_death", "retidm_ukb_death_any", "retidm_ukb_death_primary"), -eid)
  
  # Join back to original registry to ensure full sample
  data2 <- left_join(mortality_registry %>% select(eid), data2, by = "eid")
  
  return(data2)
}

# Function - Extract SR codes from self-reported diseases
extract_srd <- function(df, ages_df, disease_names) {
  name_cols <- names(df)[2:35]
  date_cols <- names(df)[36:69]
  
  # Reshape to long format
  long_df <- df %>%
    select(eid, all_of(name_cols), all_of(date_cols)) %>%
    mutate(row = row_number()) %>%
    pivot_longer(cols = -c(eid, row),
                 names_to = "col_name",
                 values_to = "value") %>%
    mutate(type = if_else(col_name %in% name_cols, "name", "date"),
           col_index = str_extract(col_name, "\\d+$") %>% as.integer())
  
  # Separate into names and dates and join by row and col_index
  long_names <- long_df %>%
    filter(type == "name") %>%
    rename(disease_name = value) %>%
    select(eid, row, col_index, disease_name)
  
  long_dates <- long_df %>%
    filter(type == "date") %>%
    rename(date_info = value) %>%
    select(eid, row, col_index, date_info)
  
  merged <- long_names %>%
    left_join(long_dates, by = c("eid", "row", "col_index")) %>%
    filter(!is.na(disease_name),
           str_to_lower(disease_name) %in% str_to_lower(disease_names))
  
  # Merge with age and date_attend
  merged <- merged %>%
    left_join(ages_df, by = "eid")  # expects columns eid, age, date_attend
  
  # Compute diagnosis date with automatic age/year detection
  merged <- merged %>%
    mutate(
      date_info_num = suppressWarnings(as.numeric(date_info)),
      diagnosis_date = case_when(
        !is.na(date_info_num) & date_info_num < 100 & !is.na(age) & !is.na(date_attend) ~ {
          year_est <- year(date_attend - as.difftime((age - date_info_num) * 365.25, units = "days"))
          as.Date(paste0(year_est, "-12-31"))
        },
        !is.na(date_info_num) & date_info_num >= 100 ~ as.Date(paste0(date_info_num, "-12-31")),
        TRUE ~ as.Date(NA)
      )
    )
  
  # Pivot wider to get one column per disease
  output <- merged %>%
    group_by(eid, disease_name) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(colname = paste0(make.names(disease_name), "_ukb_srd_date")) %>%
    select(eid, colname, diagnosis_date) %>%
    pivot_wider(names_from = colname, values_from = diagnosis_date)
  
  # Return all eids
  final <- df %>%
    select(eid) %>%
    left_join(output, by = "eid")
  
  return(final)
}

# Function - Extract SR codes from self-reported operative procedures
extract_sro <- function(df, ages_df, operation_names) {
  name_cols <- names(df)[2:33]
  date_cols <- names(df)[34:65]
  
  # Reshape to long format
  long_df <- df %>%
    select(eid, all_of(name_cols), all_of(date_cols)) %>%
    mutate(row = row_number()) %>%
    pivot_longer(cols = -c(eid, row),
                 names_to = "col_name",
                 values_to = "value") %>%
    mutate(type = if_else(col_name %in% name_cols, "name", "date"),
           col_index = str_extract(col_name, "\\d+$") %>% as.integer())
  
  # Separate into names and dates and join by row and col_index
  long_names <- long_df %>%
    filter(type == "name") %>%
    rename(operation_name = value) %>%
    select(eid, row, col_index, operation_name)
  
  long_dates <- long_df %>%
    filter(type == "date") %>%
    rename(date_info = value) %>%
    select(eid, row, col_index, date_info)
  
  merged <- long_names %>%
    left_join(long_dates, by = c("eid", "row", "col_index")) %>%
    filter(!is.na(operation_name),
           str_to_lower(operation_name) %in% str_to_lower(operation_names))
  
  # Merge with age and date_attend
  merged <- merged %>%
    left_join(ages_df, by = "eid")  # expects columns eid, age, date_attend
  
  # Compute operation date with automatic age/year detection
  merged <- merged %>%
    mutate(
      date_info_num = suppressWarnings(as.numeric(date_info)),
      operation_date = case_when(
        !is.na(date_info_num) & date_info_num < 100 & !is.na(age) & !is.na(date_attend) ~ {
          year_est <- year(date_attend - as.difftime((age - date_info_num) * 365.25, units = "days"))
          as.Date(paste0(year_est, "-12-31"))
        },
        !is.na(date_info_num) & date_info_num >= 100 ~ as.Date(paste0(date_info_num, "-12-31")),
        TRUE ~ as.Date(NA)
      )
    )
  
  # Pivot wider to get one column per operation
  output <- merged %>%
    group_by(eid, operation_name) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(colname = paste0(make.names(operation_name), "_ukb_sro_date")) %>%
    select(eid, colname, operation_date) %>%
    pivot_wider(names_from = colname, values_from = operation_date)
  
  # Return result for all eids
  final <- df %>%
    select(eid) %>%
    left_join(output, by = "eid")
  
  return(final)
}

# Function - Combine multi-source information to joint outcome variable set
derive_combined_outcome <- function(prefix,
                                    srd_data, sro_data = NULL,
                                    icd10_data, icd9_data = NULL, opcs4_data = NULL,
                                    mortality_data,
                                    meta_data,
                                    censor_date = as.Date("2022-12-31")) {
  
  extract_earliest_date <- function(df, prefix, suffix = "date") {
    stopifnot("eid" %in% colnames(df))
    date_cols <- setdiff(colnames(df), "eid")
    df[date_cols] <- lapply(df[date_cols], function(x) {
      x <- suppressWarnings(as.numeric(x))
      as.Date(x, origin = "1970-01-01")
    })
    df[[paste0(prefix, "_", suffix)]] <- do.call(pmin, c(df[date_cols], list(na.rm = TRUE)))
    df[, c("eid", paste0(prefix, "_", suffix))]
  }
  
  df <- meta_data
  df$date_attend <- as.Date(df$date_attend)
  df$date_lost   <- as.Date(df$date_lost)
  df$date_death  <- as.Date(df$date_death)
  
  srd <- extract_earliest_date(srd_data, prefix = paste0(prefix, "_srd"))
  df <- dplyr::left_join(df, srd, by = "eid")
  
  if (!is.null(sro_data)) {
    sro <- extract_earliest_date(sro_data, prefix = paste0(prefix, "_sro"))
    df <- dplyr::left_join(df, sro, by = "eid")
  } else {
    df[[paste0(prefix, "_sro_date")]] <- NA
  }
  
  icd10 <- extract_earliest_date(icd10_data, prefix = paste0(prefix, "_icd10"))
  df <- dplyr::left_join(df, icd10, by = "eid")
  
  if (!is.null(icd9_data)) {
    icd9 <- extract_earliest_date(icd9_data, prefix = paste0(prefix, "_icd9"))
    df <- dplyr::left_join(df, icd9, by = "eid")
  } else {
    df[[paste0(prefix, "_icd9_date")]] <- NA
  }
  
  if (!is.null(opcs4_data)) {
    opcs4 <- extract_earliest_date(opcs4_data, prefix = paste0(prefix, "_opcs4"))
    df <- dplyr::left_join(df, opcs4, by = "eid")
  } else {
    df[[paste0(prefix, "_opcs4_date")]] <- NA
  }
  
  # Death handling
  death_any_cols <- grep("_ukb_death_any$", names(mortality_data), value = TRUE)
  mortality_data$date_death <- as.Date(mortality_data$date_death, origin = "1970-01-01")
  death_flag <- apply(mortality_data[death_any_cols], 1, function(x) any(x == 1, na.rm = TRUE))
  death_date <- rep(NA_Date_, nrow(mortality_data))
  death_date[death_flag] <- mortality_data$date_death[death_flag]
  death_df <- data.frame(eid = mortality_data$eid, death_date)
  df <- dplyr::left_join(df, death_df, by = "eid")
  names(df)[names(df) == "death_date"] <- paste0(prefix, "_ukb_death_date")
  
  attend <- df$date_attend
  lost   <- df$date_lost
  death  <- df$date_death
  
  df$endpoint <- pmin(lost, death, censor_date, na.rm = TRUE)
  
  comphosp_date <- pmin(
    df[[paste0(prefix, "_icd10_date")]],
    df[[paste0(prefix, "_icd9_date")]],
    df[[paste0(prefix, "_opcs4_date")]],
    df[[paste0(prefix, "_ukb_death_date")]],
    na.rm = TRUE
  )
  comphosp_date[is.infinite(comphosp_date)] <- NA
  df[[paste0(prefix, "_ukb_comphosp_date")]] <- comphosp_date
  
  sr_date <- pmin(df[[paste0(prefix, "_srd_date")]],
                  df[[paste0(prefix, "_sro_date")]],
                  na.rm = TRUE)
  sr_date[is.infinite(sr_date)] <- NA
  
  # *** MODIFIED LOGIC ***
  # If self-report date is prevalent (<= baseline), it overrides comphosp_date.
  # Otherwise, fall back to comphosp_date.
  pdate <- ifelse(
    !is.na(sr_date) & sr_date <= attend,
    sr_date,
    comphosp_date
  )
  # ifelse() strips class; fix to Date
  pdate <- as.Date(pdate, origin = "1970-01-01")
  df[[paste0(prefix, "_ukb_pdate")]] <- pdate
  
  fulldate <- pmin(pdate, df$endpoint, na.rm = TRUE)
  fulldate[is.infinite(fulldate)] <- NA
  df[[paste0(prefix, "_ukb_fulldate")]] <- fulldate
  
  df[[paste0(prefix, "_ukb_prev")]] <- as.integer(!is.na(fulldate) & fulldate <= attend)
  df[[paste0(prefix, "_ukb_ind")]]  <- as.integer(!is.na(pdate))
  df[[paste0(prefix, "_ukb_survtime")]] <- round(as.numeric(fulldate - attend) / 365.25, 3)
  
  return(df)
}



#### Extract coronary artery disease (broad defintion) ####

cad_srd_names <- c("heart attack/myocardial infarction")

cad_sro_names <- c("coronary angioplasty", "coronary artery bypass grafts", "triple heart bypass")

cad_codes_icd10 <- c(
  # I21: Acute myocardial infarction
  "I21.0", "I21.1", "I21.2", "I21.3", "I21.4", "I21.9",
  # I22: Subsequent myocardial infarction
  "I22.0", "I22.1", "I22.8", "I22.9",
  # I23: Certain current complications following acute MI
  "I23.0", "I23.1", "I23.2", "I23.3", "I23.4", "I23.5", "I23.6", "I23.8",
  # I24: Other acute ischemic heart diseases
  "I24.0", "I24.1", "I24.8", "I24.9",
  # I25: Chronic ischemic heart disease
  "I25.0", "I25.1", "I25.2", "I25.3", "I25.4", "I25.5", "I25.6", "I25.8", "I25.9"
)

cad_codes_icd9 <- c(
  # 410: Acute myocardial infarction
  "4109",
  # 411: Other acute and subacute forms of ischemic heart disease
  "4119",
  # 412: Old myocardial infarction
  "4129",
  # 413: Angina pectoris
  "4139",
  # 414: Other forms of chronic ischemic heart disease
  "4140", "4141", "4148", "4149"
)

cad_codes_opcs4 <- c(
  # K40 – Saphenous vein graft replacement of coronary artery
  "K40.1", # Replacement of one coronary artery
  "K40.2", # Replacement of two coronary arteries
  "K40.3", # Replacement of three coronary arteries
  "K40.4", # Replacement of four or more coronary arteries
  "K40.8", # Other specified replacement
  "K40.9", # Unspecified replacement
  
  # K41 – Other autograft replacement of coronary artery
  "K41.1", # Replacement of one coronary artery
  "K41.2", # Replacement of two coronary arteries
  "K41.3", # Replacement of three coronary arteries
  "K41.4", # Replacement of four or more coronary arteries
  "K41.8", # Other specified replacement
  "K41.9", # Unspecified replacement
  
  # K45 – Other bypass of coronary artery
  "K45.1", # Single internal mammary–coronary artery bypass
  "K45.2", # Double internal mammary–coronary artery bypass
  "K45.3", # Other thoracic artery–coronary artery bypass
  "K45.4", # Gastroepiploic artery–coronary artery bypass
  "K45.5", # Other specified artery–coronary artery bypass
  "K45.6", # Coronary–coronary artery bypass
  "K45.8", # Other specified bypass
  "K45.9", # Unspecified bypass
  
  # K49 – Other graft replacement of coronary artery
  "K49.1", # Replacement of one coronary artery
  "K49.2", # Replacement of two coronary arteries
  "K49.3", # Replacement of three coronary arteries
  "K49.4", # Replacement of four or more coronary arteries
  "K49.8", # Other specified replacement
  "K49.9", # Unspecified replacement
  
  # K50 – Anastomosis of coronary artery
  "K50.2", # Other specified anastomosis of coronary artery
  
  # K75 – Other therapeutic transluminal operations on coronary artery
  "K75.1", # Percutaneous transluminal balloon angioplasty (PTCA)
  "K75.2", # PTCA + stent insertion
  "K75.3", # Primary transluminal coronary stent insertion
  "K75.4", # Other specified transluminal coronary stent insertion
  "K75.8", # Other specified transluminal operation
  "K75.9"  # Unspecified transluminal operation
)


cad_srd <- extract_srd(df = srd_registry, 
                       ages_df = process_data, 
                       disease_names = c("heart attack/myocardial infarction")
)

cad_sro <- extract_sro(df = sro_registry, 
                       ages_df = process_data, 
                       operation_names = c("coronary angioplasty (ptca) +/- stent", "coronary artery bypass grafts (cabg)", "triple heart bypass")
)

cad_icd10 <- extract_icd10(cad_codes_icd10, icd10_registry)

cad_icd9 <- extract_icd9(cad_codes_icd9, icd9_registry)

cad_opcs4 <- extract_opcs4(cad_codes_opcs4, opcs4_registry)

cad_death <- extract_death(cad_codes_icd10, mortality_registry)

cad_output <- derive_combined_outcome(prefix = "cad",
                                      srd_data = cad_srd, 
                                      sro_data = cad_sro,
                                      icd10_data = cad_icd10,
                                      icd9_data = cad_icd9,
                                      opcs4_data = cad_opcs4,
                                      mortality_data = cad_death,
                                      meta_data = process_data, # must have eid, date_attend, date_lost, date_death
                                      censor_date = as.Date("2022-12-31")
)

cad_output <- cad_output %>%
  select(eid, cad_ukb_fulldate:cad_ukb_survtime)

save(cad_output, file = "cad_output.rda")

#### Extract coronary artery disease (narrow definition) ####

cad2_srd_names <- c("heart attack/myocardial infarction")

cad2_sro_names <- c("coronary angioplasty", "coronary artery bypass grafts", "triple heart bypass")

cad2_codes_icd10 <- c(
  # I21 – Acute myocardial infarction (MI)
  "I21.0", # Acute transmural MI of anterior wall
  "I21.1", # Acute transmural MI of inferior wall
  "I21.2", # Acute transmural MI of other sites
  "I21.3", # Acute transmural MI of unspecified site
  "I21.4", # Acute subendocardial MI
  "I21.9", # Acute MI, unspecified
  
  # I22 – Subsequent myocardial infarction
  "I22.0", # Subsequent MI of anterior wall
  "I22.1", # Subsequent MI of inferior wall
  "I22.8", # Subsequent MI of other sites
  "I22.9", # Subsequent MI of unspecified site
  
  # I23 – Certain current complications following acute MI
  "I23.0", # Hemopericardium as current complication
  "I23.1", # Atrial septal defect as current complication
  "I23.2", # Ventricular septal defect as current complication
  "I23.3", # Rupture of cardiac wall without hemopericardium
  "I23.5", # Rupture of chordae tendineae
  "I23.6", # Rupture of papillary muscle
  "I23.8", # Other current complications following acute MI
  
  # I24 – Other acute ischemic heart diseases
  "I24.1", # Dressler’s syndrome (post-MI pericarditis)
  
  # I25 – Chronic ischemic heart disease
  "I25.2"  # Old myocardial infarction
)

cad2_codes_icd9 <- c(
  # 410: Acute myocardial infarction
  "4109",
  # 411: Other acute and subacute forms of ischemic heart disease
  "4119",
  # 412: Old myocardial infarction
  "4129"
)

cad2_codes_opcs4 <- c(
  # K40 – Saphenous vein graft replacement of coronary artery
  "K40.1", # Replacement of one coronary artery
  "K40.2", # Replacement of two coronary arteries
  "K40.3", # Replacement of three coronary arteries
  "K40.4", # Replacement of four or more coronary arteries
  "K40.8", # Other specified replacement
  "K40.9", # Unspecified replacement
  
  # K41 – Other autograft replacement of coronary artery
  "K41.1", # Replacement of one coronary artery
  "K41.2", # Replacement of two coronary arteries
  "K41.3", # Replacement of three coronary arteries
  "K41.4", # Replacement of four or more coronary arteries
  "K41.8", # Other specified replacement
  "K41.9", # Unspecified replacement
  
  # K45 – Other bypass of coronary artery
  "K45.1", # Single internal mammary–coronary artery bypass
  "K45.2", # Double internal mammary–coronary artery bypass
  "K45.3", # Other thoracic artery–coronary artery bypass
  "K45.4", # Gastroepiploic artery–coronary artery bypass
  "K45.5", # Other specified artery–coronary artery bypass
  "K45.6", # Coronary–coronary artery bypass
  "K45.8", # Other specified bypass
  "K45.9", # Unspecified bypass
  
  # K49 – Other graft replacement of coronary artery
  "K49.1", # Replacement of one coronary artery
  "K49.2", # Replacement of two coronary arteries
  "K49.3", # Replacement of three coronary arteries
  "K49.4", # Replacement of four or more coronary arteries
  "K49.8", # Other specified replacement
  "K49.9", # Unspecified replacement
  
  # K50 – Anastomosis of coronary artery
  "K50.2", # Other specified anastomosis of coronary artery
  
  # K75 – Other therapeutic transluminal operations on coronary artery
  "K75.1", # Percutaneous transluminal balloon angioplasty (PTCA)
  "K75.2", # PTCA + stent insertion
  "K75.3", # Primary transluminal coronary stent insertion
  "K75.4", # Other specified transluminal coronary stent insertion
  "K75.8", # Other specified transluminal operation
  "K75.9"  # Unspecified transluminal operation
)


cad2_srd <- extract_srd(df = srd_registry, 
                        ages_df = process_data, 
                        disease_names = c("heart attack/myocardial infarction")
)

cad2_sro <- extract_sro(df = sro_registry, 
                        ages_df = process_data, 
                        operation_names = c("coronary angioplasty (ptca) +/- stent", "coronary artery bypass grafts (cabg)", "triple heart bypass")
)

cad2_icd10 <- extract_icd10(cad2_codes_icd10, icd10_registry)

cad2_icd9 <- extract_icd9(cad2_codes_icd9, icd9_registry)

cad2_opcs4 <- extract_opcs4(cad2_codes_opcs4, opcs4_registry)

cad2_death <- extract_death(cad2_codes_icd10, mortality_registry)

cad2_output <- derive_combined_outcome(prefix = "cad2",
                                       srd_data = cad2_srd, 
                                       sro_data = cad2_sro,
                                       icd10_data = cad2_icd10,
                                       icd9_data = cad2_icd9,
                                       opcs4_data = cad2_opcs4,
                                       mortality_data = cad2_death,
                                       meta_data = process_data, # must have eid, date_attend, date_lost, date_death
                                       censor_date = as.Date("2022-12-31")
)

cad2_output <- cad2_output %>%
  select(eid, cad2_ukb_fulldate:cad2_ukb_survtime)

save(cad2_output, file = "cad2_output.rda")

#### Extract total stroke ####

stroke_srd_names <- c(
  "stroke",                    # Field 20002 Code 1081
  "subarachnoid haemorrhage",  # Field 20002 Code 1086
  "brain haemorrhage",         # Field 20002 Code 1491
  "ischaemic stroke"           # Field 20002 Code 1583
)

stroke_codes_icd10 <- c(
  # Subarachnoid haemorrhage
  "I60.0", "I60.1", "I60.2", "I60.3", "I60.4", "I60.5",
  "I60.6", "I60.7", "I60.8", "I60.9",
  
  # Intracerebral haemorrhage
  "I61.0", "I61.1", "I61.2", "I61.3", "I61.4",
  "I61.5", "I61.6", "I61.8", "I61.9",
  
  # Cerebral infarction
  "I63.0", "I63.1", "I63.2", "I63.3", "I63.4",
  "I63.5", "I63.6", "I63.8", "I63.9",
  
  # Stroke NOS
  "I64"    # Stroke, not specified as haemorrhage or infarction
)

# Stroke: ICD-9
stroke_codes_icd9 <- c(
  # Subarachnoid haemorrhage
  "4309",
  # Intracerebral haemorrhage
  "4319",
  # Occlusion of cerebral arteries
  "4349",
  # Acute but ill-defined cerebrovascular disease
  "4369"
)

# Stroke extraction
stroke_srd <- extract_srd(df = srd_registry, 
                          ages_df = process_data, 
                          disease_names = stroke_srd_names
)

stroke_icd10 <- extract_icd10(stroke_codes_icd10, icd10_registry)

stroke_icd9 <- extract_icd9(stroke_codes_icd9, icd9_registry)

stroke_death <- extract_death(stroke_codes_icd10, mortality_registry)

stroke_output <- derive_combined_outcome(prefix = "stroke",
                                         srd_data = stroke_srd, 
                                         icd10_data = stroke_icd10,
                                         icd9_data = stroke_icd9,
                                         mortality_data = stroke_death,
                                         meta_data = process_data,
                                         censor_date = as.Date("2022-12-31")
)

stroke_output <- stroke_output %>%
  select(eid, stroke_ukb_fulldate:stroke_ukb_survtime)

save(stroke_output, file = "stroke_output.rda")

#### Extract ischaemic stroke ####

# Ischaemic Stroke: UK Biobank Self-report
ischemic_stroke_srd_names <- c(
  "ischaemic stroke"   # Field 20002 Code 1583
)

# Ischaemic Stroke: ICD-10
ischemic_stroke_codes_icd10 <- c(
  # Cerebral infarction
  "I63.0", "I63.1", "I63.2", "I63.3", "I63.4",
  "I63.5", "I63.6", "I63.8", "I63.9",
  
  # Stroke NOS (classified as ischaemic in your table)
  "I64"
)

# Ischaemic Stroke: ICD-9
ischemic_stroke_codes_icd9 <- c(
  # Occlusion of cerebral arteries
  "4349", 
  # Acute but ill-defined cerebrovascular disease (treated as ischaemic)
  "4369"
)

# Ischaemic Stroke extraction
ischemic_stroke_srd <- extract_srd(df = srd_registry, 
                                   ages_df = process_data, 
                                   disease_names = ischemic_stroke_srd_names
)

ischemic_stroke_icd10 <- extract_icd10(ischemic_stroke_codes_icd10, icd10_registry)

ischemic_stroke_icd9 <- extract_icd9(ischemic_stroke_codes_icd9, icd9_registry)

ischemic_stroke_death <- extract_death(ischemic_stroke_codes_icd10, mortality_registry)

ischemic_stroke_output <- derive_combined_outcome(prefix = "ischemic_stroke",
                                                  srd_data = ischemic_stroke_srd, 
                                                  icd10_data = ischemic_stroke_icd10,
                                                  icd9_data = ischemic_stroke_icd9,
                                                  mortality_data = ischemic_stroke_death,
                                                  meta_data = process_data,
                                                  censor_date = as.Date("2022-12-31")
)

ischemic_stroke_output <- ischemic_stroke_output %>%
  select(eid, ischemic_stroke_ukb_fulldate:ischemic_stroke_ukb_survtime)

save(ischemic_stroke_output, file = "ischemic_stroke_output.rda")

#### Extract peripheral artery disease ####

# PAD: UK Biobank Self-report (names, not codes)
pad_srd_names <- c(
  "peripheral vascular disease",
  "leg claudication/ intermittent claudication",
  "arterial embolism"
)

pad_sro_names <- c(
  "fem-pop bypass/leg artery bypass",
  "leg artery angioplasty",
  "amputation of leg"
)


# PAD: ICD-10 codes (4–5 digit)
pad_codes_icd10 <- c(
  # I70
  "I70.0", "I70.00", "I70.01",
  "I70.2", "I70.20", "I70.21",
  "I70.8", "I70.80",
  "I70.9", "I70.90",
  # I73
  "I73.8", "I73.9"
)

# PAD: ICD-9 codes 
pad_codes_icd9 <- c(
  # 440
  "4400", "4402",
  # 443
  "4438", "4439"
)

# PAD: OPCS-4 procedure codes
pad_codes_opcs4 <- c(
  # L21
  "L21.6",
  # L51
  "L51.3", "L51.6", "L51.8",
  # L52
  "L52.1", "L52.2",
  # L54
  "L54.1", "L54.4", "L54.8",
  # L59
  "L59.1", "L59.2", "L59.3", "L59.4", "L59.5", "L59.6", "L59.7", "L59.8",
  # L60
  "L60.1", "L60.2",
  # L63
  "L63.1", "L63.5", "L63.9",
  # L66
  "L66.7",
  # X09
  "X09.3", "X09.4", "X09.5"
)


# PAD extraction
pad_srd <- extract_srd(df = srd_registry, 
                       ages_df = process_data, 
                       disease_names = pad_srd_names
)

pad_sro <- extract_sro(df = sro_registry, 
                       ages_df = process_data, 
                       operation_names = pad_sro_names
)

pad_icd10 <- extract_icd10(pad_codes_icd10, icd10_registry)

pad_icd9 <- extract_icd9(pad_codes_icd9, icd9_registry)

pad_opcs4 <- extract_opcs4(pad_codes_opcs4, opcs4_registry)

pad_death <- extract_death(pad_codes_icd10, mortality_registry)

pad_output <- derive_combined_outcome(prefix = "pad",
                                      srd_data = pad_srd, 
                                      sro_data = pad_sro,
                                      icd10_data = pad_icd10,
                                      icd9_data = pad_icd9,
                                      opcs4_data = pad_opcs4,
                                      mortality_data = pad_death,
                                      meta_data = process_data,
                                      censor_date = as.Date("2022-12-31")
)

pad_output <- pad_output %>%
  select(eid, pad_ukb_fulldate:pad_ukb_survtime)

save(pad_output, file = "pad_output.rda")



#### Extract type 2 diabetes ####

dm2_codes_icd10 <- c(
  "E11.0",    # Type 2 diabetes mellitus with hyperosmolarity
  "E11.1",    # Type 2 diabetes mellitus with ketoacidosis
  "E11.2",    # Type 2 diabetes mellitus with kidney complications
  "E11.3",    # Type 2 diabetes mellitus with ophthalmic complications
  "E11.4",    # Type 2 diabetes mellitus with neurological complications
  "E11.5",    # Type 2 diabetes mellitus with peripheral circulatory complications
  "E11.6",    # Type 2 diabetes mellitus with other specified complications
  "E11.7",    # Type 2 diabetes mellitus with multiple complications
  "E11.8",    # Type 2 diabetes mellitus with unspecified complications
  "E11.9"     # Type 2 diabetes mellitus without complications
)


dm2_codes_icd9 <- c("25000", "25010")


dm2_srd <- extract_srd(df = srd_registry, 
                       ages_df = process_data, 
                       disease_names = c("type 2 diabetes")
)


dm2_icd10 <- extract_icd10(dm2_codes_icd10, icd10_registry)

dm2_icd9 <- extract_icd9(dm2_codes_icd9, icd9_registry)

dm2_death <- extract_death(dm2_codes_icd10, mortality_registry)


dm2_output <- derive_combined_outcome(prefix = "dm2",
                                      srd_data = dm2_srd, 
                                      icd10_data = dm2_icd10,
                                      icd9_data = dm2_icd9,
                                      mortality_data = dm2_death,
                                      meta_data = process_data,
                                      censor_date = as.Date("2022-12-31")
)


dm2_output <- dm2_output %>%
  select(eid, dm2_ukb_fulldate:dm2_ukb_survtime)

save(dm2_output, file = "dm2_output.rda")


##### MEDICATION DATA PREPARATION #####

# Load individual level medication use
medications <- read_csv("medications.csv", guess_max = 502128)

# Load class coding file
ukbb_drug_coding <- read_excel("ukbb-drug-coding.xlsx")

## Load extraction functions

## 1) Build the exact label set for a group from the coding table
get_group_labels_exact <- function(coding,
                                   group_col,          # e.g. "STATIN"
                                   label_col = "Category") {
  if (!group_col %in% names(coding)) stop(sprintf("'%s' not in coding.", group_col), call. = FALSE)
  if (!label_col %in% names(coding)) stop(sprintf("'%s' not in coding.", label_col), call. = FALSE)
  
  keep <- !is.na(coding[[group_col]]) & coding[[group_col]] == group_col
  labs <- coding[[label_col]][keep]
  labs <- as.character(labs)
  unique(labs[!is.na(labs) & nzchar(labs)])
}

## 2) Compute 0/1 indicator: any of the selected columns equals any label (EXACT match)
med_indicator_exact <- function(prescriptions,
                                labels,     # character vector of exact labels
                                drug_cols,  # character vector: your 48 column names
                                use_fastmatch = TRUE) {
  if (length(labels) == 0L) return(rep.int(0L, nrow(prescriptions)))
  if (!all(drug_cols %in% names(prescriptions))) {
    miss <- setdiff(drug_cols, names(prescriptions))
    stop("Missing drug_cols in prescriptions: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  M <- as.matrix(prescriptions[, drug_cols, drop = FALSE])
  storage.mode(M) <- "character"
  M[is.na(M)] <- ""
  
  if (use_fastmatch && requireNamespace("fastmatch", quietly = TRUE)) {
    L <- vapply(seq_len(ncol(M)),
                function(j) fastmatch::fmatch(M[, j], labels, nomatch = 0L) > 0L,
                logical(nrow(M)))
  } else {
    L <- vapply(seq_len(ncol(M)),
                function(j) M[, j] %in% labels,
                logical(nrow(M)))
  }
  as.integer(rowSums(L) > 0L)
}

## 3) Convenience wrapper: add the indicator as a new column
add_med_indicator_exact <- function(prescriptions,
                                    coding,
                                    group_col,              # e.g. "STATIN"
                                    out_col   = group_col,  # output column name
                                    drug_cols,              # the 48 prescription columns
                                    label_col = "Category",
                                    use_fastmatch = TRUE) {
  labels <- get_group_labels_exact(coding, group_col, label_col)
  prescriptions[[out_col]] <- med_indicator_exact(prescriptions, labels, drug_cols, use_fastmatch)
  prescriptions
}


# 1) Select 48 prescription columns by name pattern (base R)
drug_cols <- grep("^Treatment/medication code \\| Instance 0 \\| Array [0-9]+$",
                  names(medications), value = TRUE)

## Extract lipid medications
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "LIPID",
  out_col       = "med_lipid",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract statin medications
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "STATIN",
  out_col       = "med_statin",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract statin medications
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "DIAB",
  out_col       = "med_diab",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract ACE inhibitor use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "ACE",
  out_col       = "med_ace",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract CRB use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "CRB",
  out_col       = "med_crb",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract BB use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "BB",
  out_col       = "med_bb",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract ARB use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "ARB",
  out_col       = "med_arb",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract THI use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "THI",
  out_col       = "med_thi",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract LOP use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "LOP",
  out_col       = "med_lop",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)

## Extract POT use
medications <- add_med_indicator_exact(
  prescriptions = medications,
  coding        = ukbb_drug_coding,  # your coding table
  group_col     = "POT",
  out_col       = "med_pot",      # name the output column as you like
  drug_cols     = drug_cols,
  label_col     = "Category"
)


medications_output <- medications %>%
  mutate(
    # Hypertension meds: ACEi, ARB, BB, CCB, Thiazide, Potassium-sparing, Loop
    med_hyperten = case_when(
      med_ace == 1 ~ 1L,
      med_arb == 1 ~ 1L,
      med_bb  == 1 ~ 1L,
      med_crb == 1 ~ 1L,  # CCB (rename if your column is med_ccb)
      med_thi == 1 ~ 1L,  # Thiazide
      med_pot == 1 ~ 1L,  # Potassium-sparing
      med_lop == 1 ~ 1L,  # Loop
      TRUE ~ 0L
    ),
    # Manual factor conversion for each variable
    med_lipid   = factor(med_lipid,   levels = c(0, 1), labels = c("No", "Yes")),
    med_statin  = factor(med_statin,  levels = c(0, 1), labels = c("No", "Yes")),
    med_diab    = factor(med_diab,    levels = c(0, 1), labels = c("No", "Yes")),
    med_hyperten= factor(med_hyperten,levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  select(`Participant ID`, med_statin, med_lipid, med_diab, med_hyperten) %>% 
  rename(eid = `Participant ID`)

# Save the file
save(medications_output, file="medications_output.rda")

##### BUILD WORK DATASET #####

# Load previosuly extracted data
load("basic_covariates_output.rda")
load("olink_filtered_imputed.rda")
load("cad_output.rda")
load("cad2_output.rda")
load("stroke_output.rda")
load("ischemic_stroke_output.rda")
load("pad_output.rda")
load("dm2_output.rda")
load("medications_output.rda")

# Add lipoprotein(a) values from Return 2321
lpa_values_outofrange_rebadged_570811_ <- read_csv("lpa_values_outofrange[rebadged_570811].csv")


lpa_values_outofrange_rebadged_570811_ <- lpa_values_outofrange_rebadged_570811_ %>% rename(eid = eid_570811,
                                                                                            conven_lpa_corr = LPA_oval_b)
# Combine with basic covariates
basic_covariates_output <- left_join(basic_covariates_output, lpa_values_outofrange_rebadged_570811_ %>% select(eid, conven_lpa_corr), by = "eid")

# Merge datasets (except for proteomic)
datasets_list <- list(basic_covariates_output, cad_output, cad2_output, stroke_output, ischemic_stroke_output, pad_output, medications_output)

dataset_temp <- Reduce(full_join, datasets_list)

# Right join with proteomic data - N = 502,128 => N = 44,782
dataset <- right_join(dataset_temp, olink_filtered_imputed, by ="eid")

# Exclude prevalent CAD, stroke and PAD (2014 broad CAD, 744 stroke, 364 PAD ) N = 44,782 => N = 41,951 (2831 removed)
dataset <- dataset %>% filter(cad_ukb_prev == 0 & stroke_ukb_prev == 0 & pad_ukb_prev == 0)

# Exclude missing total cholesterol (N = 2067), direct LDL (N = 2154), HDL cholesterol (N = 5438) or lipoprotein(a) [conven_lpa_corr] (N = 3092) N = 41,951 => N = 35,269 (6682 removed)
dataset <- dataset %>% filter(!is.na(conven_tchol) & !is.na(conven_ldl) & !is.na(conven_hdl) & !is.na(conven_lpa_corr))

# Calculate necessary lipid variables
dataset <- dataset %>% 
  mutate(
    conven_rchol = conven_tchol - conven_hdl - conven_ldl,
    conven_rchol = ifelse(conven_rchol < 0, 0, conven_rchol),
    conven_lpa_corr = ifelse(conven_lpa_corr < 0, 0, conven_lpa_corr),
    conven_trlp = conven_rchol * 250,
    conven_ldlp = conven_ldl * 500,
    conven_trlp_100 = conven_trlp / 100,
    conven_ldlp_100 = conven_ldlp / 100,
    conven_lpa_corr_100 = conven_lpa_corr / 100
  )

# Generate a report of variable metadata
generate_report <- function(df) {
  report <- data.frame(
    Variable = names(df),
    Type = sapply(df, function(x) class(x)[1]),
    Missing_NA = sapply(df, function(x) sum(is.na(x))),
    stringsAsFactors = FALSE
  )
  return(report)
}

metadata <- generate_report(dataset)
#> print(metadata %>% filter(Missing_NA>0))
#Variable    Type Missing_NA
#bmi                       bmi numeric        146
#race                     race  factor        150
#smoking_status smoking_status  factor        154
#conven_tg           conven_tg numeric          3
#conven_lpa         conven_lpa numeric       5996
#conven_apob       conven_apob numeric        203
#conven_hba1c     conven_hba1c numeric       1730
#conven_hscrp     conven_hscrp numeric        119
#conven_creat     conven_creat numeric          1
#gfr_mdrd             gfr_mdrd numeric          1
#blood_wbc           blood_wbc numeric        881
#blood_rdw           blood_rdw numeric        880
#blood_lympho     blood_lympho numeric        936
#blood_neutro     blood_neutro numeric        936
#sbp_mean             sbp_mean numeric         32
#dbp_mean             dbp_mean numeric         32

## Missing "race" will be imputed with "Mixed/other category"
dataset <- dataset %>% mutate(race_imp = fct_na_value_to_level(race, level = "Mixed/Other"))

## Missing "conven_apob" will be imputed from LDL cholesterol and log(triglycerides) - All participants with missing apoB have LDL-C and TG values
summary(lm(conven_apob ~ conven_ldl + log(conven_tg+1), data = dataset))

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)        0.0419479  0.0016556   25.34   <2e-16 ***
#  conven_ldl         0.2604738  0.0004268  610.29   <2e-16 ***
#  log(conven_tg + 1) 0.0638908  0.0011692   54.65   <2e-16 ***
#  Multiple R-squared:  0.9243,	Adjusted R-squared:  0.9243 
#  F-statistic: 2.139e+05 on 2 and 35060 DF,  p-value: < 2.2e-16

dataset <- dataset %>% mutate(conven_apob_imp = ifelse(is.na(conven_apob), 0.0419479 + conven_ldl * 0.2604738 + log(conven_tg) * 0.0638908, conven_apob))

## bmi, smokign_status, conven_tg, conve3n_hba1c, gfr_mdrd, sbp_mean and. dbp_mean using k-nearest neighbors 

# Variables to impute (deduped)
impute_vars <- c("bmi", "smoking_status", "conven_hba1c", "gfr_mdrd", "sbp_mean", "dbp_mean")

# Variables to use for distance (deduped)
dist_vars <- c(
  "age", "sex", "race", 
  "conven_ldl", "conven_rchol", "conven_lpa_corr", "med_lipid", "med_hyperten", "med_diab",
  "sbp_mean", "dbp_mean"
)

# kNN imputation (k=5 by default).
set.seed(123)
dataset_imputed <- VIM::kNN(
  data     = dataset,
  variable = impute_vars,
  dist_var = dist_vars,
  k        = 5,
  imp_var  = FALSE,
  imp_suffix = TRUE,
  weightDist = FALSE
)

# Selected, rename and merge imputed variables
dataset_imputed_temp <- dataset_imputed %>%
  select(
    eid,
    bmi,
    smoking_status,
    conven_hba1c,
    gfr_mdrd,
    sbp_mean,
    dbp_mean
  ) %>%
  rename(
    bmi_imp          = bmi,
    smoking_status_imp = smoking_status,
    conven_hba1c_imp = conven_hba1c,
    gfr_mdrd_imp     = gfr_mdrd,
    sbp_mean_imp     = sbp_mean,
    dbp_mean_imp     = dbp_mean
  )

dataset <- left_join(dataset, dataset_imputed_temp, by = "eid")

## Remove OLink APOB100 and LPA from the work dataset and save them separately

dataset_side1 <- dataset %>% 
  select(eid, `APOB;Apolipoprotein B-100`, `LPA;Apolipoprotein(a)`) %>% 
  rename(APOB = `APOB;Apolipoprotein B-100`,
         LPA = `LPA;Apolipoprotein(a)`)

save(dataset_side1, file="dataset_side1.rda")

dataset <- dataset %>% select(-`APOB;Apolipoprotein B-100`, -`LPA;Apolipoprotein(a)`)

## Change protein names to short variants (gene column from name_uniprot_gene.csv)
name_uniprot_gene <- read_excel("name_uniprot_gene_pathway.xlsx")

# Create a lookup table
name_map <- setNames(name_uniprot_gene$gene,
                     name_uniprot_gene$protein_name)

# Rename columns
colnames(dataset) <- ifelse(
  colnames(dataset) %in% names(name_map),
  name_map[colnames(dataset)],
  colnames(dataset)
)

##### SAVE FINAL WORK DATASET #####
save(dataset, file = "dataset.rda")

