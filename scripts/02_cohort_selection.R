###############################################################################
# 02_cohort_selection.R
#
# Build semaglutide and tirzepatide treatment cohorts from the Helix Research
# Network (HRN) OMOP CDM database.
#
# Steps:
#   1. Extract body weight and height measurements
#   2. Extract demographics and genetic similarity (ancestry)
#   3. Extract drug exposure records (semaglutide, tirzepatide, other GLP-1 RAs)
#   4. Build prescription pattern indicators
#   5. Extract comorbidities (type 2 diabetes, malignancy, bariatric surgery,
#      pregnancy)
#   6. Construct censoring dates
#   7. Calculate baseline BMI
#   8. Apply inclusion/exclusion criteria
#   9. Save cohort files and longitudinal weight data
#
# Outputs:
#   output/semaglutide_cohort.RData
#   output/tirzepatide_cohort.RData
#   output/semaglutide_longitudinal_weights.RData
#   output/tirzepatide_longitudinal_weights.RData
#
# IMPORTANT: This script runs on Helix Research Network infrastructure.
# Database connection, query functions, and concept ID lists are specific to
# this environment and cannot be executed externally.
# See Methods for data sources and access details.
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(lubridate)
library(reticulate)
library(RNOmni)    # for RankNorm()

source("config.R")

###############################################################################
# HELIX INTERNAL: DATABASE CONNECTION
# The following block connects to the HRN OMOP CDM database via an internal
# API. This requires Helix Research Network credentials and environment.
###############################################################################

# [HELIX INTERNAL: activate environment and load internal database API]
# use_condaenv('[HELIX INTERNAL: conda environment path]')
# py_run_file('[HELIX INTERNAL: database API helper path]')
# [HELIX INTERNAL: set database role]

# query() — wrapper that executes SQL against the HRN OMOP CDM
# [HELIX INTERNAL: defined in internal API; not reproduced here]

###############################################################################
# OMOP CONCEPT ID LISTS
# Standard OMOP CDM concept IDs used to identify clinical entities.
# [HELIX INTERNAL: full concept ID lists curated from the HRN OMOP vocabulary.
#  Key concept domains used: Measurement (weight, height), Condition
#  (malignancy, T2D, pregnancy), Observation (bariatric surgery, malignancy,
#  pregnancy), Procedure (bariatric surgery).
#  See Methods for concept selection criteria.]
###############################################################################

# weight_concept_id_list      — OMOP concept IDs for body weight measurements
# height_concept_id_list      — OMOP concept IDs for body height measurements
# malignancy_concept_id_list  — OMOP concept IDs for malignancy diagnoses
# bariatric_surgery_concept_id_list — OMOP concept IDs for bariatric procedures
# pregnancy_concept_id_list   — OMOP concept IDs for pregnancy diagnoses
# t2diabetes_concept_id_list  — OMOP concept IDs for type 2 diabetes diagnoses

# [HELIX INTERNAL: load concept ID lists from internal vocabulary resource]
# source("[HELIX INTERNAL: concept_id_lists path]")

###############################################################################
# DATABASE QUERY HELPER FUNCTIONS
###############################################################################

#' Query the OMOP measurement table for specified concept IDs
query_measurement <- function(concept_ids) {
  where_clause <- paste0(
    "WHERE a.measurement_concept_id IN (", concept_ids, ")
     OR a.measurement_source_concept_id IN (", concept_ids, ")"
  )
  SQL <- paste0("
    SELECT a.measurement_concept_id,
           a.measurement_date,
           a.value_as_number,
           a.value_as_concept_id,
           a.unit_concept_id,
           a.measurement_source_value,
           a.measurement_source_concept_id,
           a.unit_source_value,
           a.unit_source_concept_id,
           a.value_source_value,
           p.person_source_value,
           p.site,
           c1.concept_name AS measurement_concept_name,
           c2.concept_name AS value_as_concept_name,
           c3.concept_name AS unit_concept_name
    FROM measurement a
    INNER JOIN person p ON a.person_id = p.person_id
    LEFT JOIN concept c1 ON a.measurement_concept_id = c1.concept_id
    LEFT JOIN concept c2 ON a.value_as_concept_id   = c2.concept_id
    LEFT JOIN concept c3 ON a.unit_concept_id        = c3.concept_id
    ", where_clause, ";"
  )
  query(SQL) %>% mutate(measurement_date = as.Date(measurement_date))
}

#' Query drug_exposure via ATC hierarchy
#' @param atc_level  ATC level (integer, typically 5 for drug-level code)
#' @param atc_codes  Character vector of ATC codes
get_drug <- function(atc_level, atc_codes) {
  atc_column   <- paste0("atc", atc_level, "_concept_code")
  atc_codes_sql <- paste0("('", paste(atc_codes, collapse = "','"), "')")

  SQL <- paste0("
    SELECT a.drug_concept_id,
           a.drug_exposure_start_date,
           a.drug_exposure_end_date,
           a.days_supply,
           a.quantity,
           a.refills,
           a.sig,
           a.drug_source_value,
           a.drug_source_concept_id,
           a.data_source,
           p.site,
           p.person_source_value,
           c1.concept_name AS drug_concept_name,
           c2.concept_name AS drug_type_concept_name,
           c3.concept_name AS route_concept_name,
           c4.concept_name AS drug_source_concept_name
    FROM drug_exposure a
    INNER JOIN person p ON a.person_id = p.person_id
    LEFT JOIN concept c1 ON a.drug_concept_id        = c1.concept_id
    LEFT JOIN concept c2 ON a.drug_type_concept_id   = c2.concept_id
    LEFT JOIN concept c3 ON a.route_concept_id       = c3.concept_id
    LEFT JOIN concept c4 ON a.drug_source_concept_id = c4.concept_id
    INNER JOIN (
      SELECT DISTINCT ds.drug_concept_id
      FROM drug_strength ds
      INNER JOIN atc_drug_reference ah ON ah.ingredient_concept_id = ds.ingredient_concept_id
      WHERE ah.", atc_column, " IN ", atc_codes_sql, "
    ) sub ON a.drug_concept_id = sub.drug_concept_id;"
  )
  query(SQL) %>%
    mutate(drug_exposure_start_date = as.Date(drug_exposure_start_date),
           drug_exposure_end_date   = as.Date(drug_exposure_end_date)) %>%
    arrange(person_source_value, drug_exposure_start_date)
}

#' Query condition_occurrence in chunks (for large concept ID lists)
query_condition_chunks <- function(concept_ids) {
  # [HELIX INTERNAL: chunk() splits large concept ID lists for query efficiency]
  chunks <- chunk(concept_ids)
  result <- data.frame()
  for (ch in chunks) {
    SQL <- paste0("
      SELECT a.condition_concept_id,
             a.condition_start_date,
             a.condition_source_value,
             a.condition_source_concept_id,
             p.person_source_value,
             c1.concept_name AS condition_concept_name,
             c2.concept_name AS condition_source_concept_name
      FROM condition_occurrence a
      INNER JOIN person p ON a.person_id = p.person_id
      LEFT JOIN concept c1 ON a.condition_concept_id        = c1.concept_id
      LEFT JOIN concept c2 ON a.condition_source_concept_id = c2.concept_id
      WHERE a.condition_concept_id IN (", ch, ")
         OR a.condition_source_concept_id IN (", ch, ");"
    )
    result <- rbind(result, query(SQL))
  }
  result %>% mutate(condition_start_date = as.Date(condition_start_date))
}

#' Query observation table in chunks
query_observation_chunks <- function(concept_ids) {
  chunks <- chunk(concept_ids)
  result <- data.frame()
  for (ch in chunks) {
    SQL <- paste0("
      SELECT a.observation_concept_id,
             a.observation_date,
             a.observation_source_value,
             a.observation_source_concept_id,
             p.person_source_value,
             c1.concept_name AS observation_concept_name,
             c2.concept_name AS observation_source_concept_name
      FROM observation a
      INNER JOIN person p ON a.person_id = p.person_id
      LEFT JOIN concept c1 ON a.observation_concept_id        = c1.concept_id
      LEFT JOIN concept c2 ON a.observation_source_concept_id = c2.concept_id
      WHERE a.observation_concept_id IN (", ch, ")
         OR a.observation_source_concept_id IN (", ch, ");"
    )
    result <- rbind(result, query(SQL))
  }
  result %>% mutate(observation_date = as.Date(observation_date))
}

#' Query procedure_occurrence
query_procedure <- function(concept_ids) {
  where_clause <- paste0(
    "WHERE a.procedure_concept_id IN (", concept_ids, ")
     OR a.procedure_source_concept_id IN (", concept_ids, ")"
  )
  SQL <- paste0("
    SELECT a.procedure_concept_id,
           a.procedure_date,
           a.procedure_source_value,
           a.procedure_source_concept_id,
           p.person_source_value,
           c1.concept_name AS procedure_concept_name
    FROM procedure_occurrence a
    INNER JOIN person p ON a.person_id = p.person_id
    LEFT JOIN concept c1 ON a.procedure_concept_id = c1.concept_id
    ", where_clause, ";"
  )
  query(SQL) %>% mutate(procedure_date = as.Date(procedure_date))
}


###############################################################################
# 1. WEIGHT AND HEIGHT EXTRACTION
###############################################################################

weight_raw <- query_measurement(weight_concept_id_list)
height_raw <- query_measurement(height_concept_id_list)

###############################################################################
# 2. DEMOGRAPHICS AND GENETIC SIMILARITY
###############################################################################

SQL_person <- "
  SELECT a.*,
         c1.concept_name AS sex_ehr,
         c2.concept_name AS race,
         c3.concept_name AS ethnicity
  FROM person a
  LEFT JOIN concept c1 ON a.gender_concept_id    = c1.concept_id
  LEFT JOIN concept c2 ON a.race_concept_id      = c2.concept_id
  LEFT JOIN concept c3 ON a.ethnicity_concept_id = c3.concept_id;"

person <- query(SQL_person) %>%
  mutate(
    birth_date = as_date(birth_datetime),
    sex_ehr = case_when(
      sex_ehr == "No matching concept" ~ "Unknown",
      TRUE ~ str_to_title(sex_ehr)
    ),
    race_eth = case_when(
      ethnicity == "Hispanic or Latino"              ~ "Hispanic",
      race %in% c("Asian", "Chinese", "Filipino")   ~ "Asian",
      race %in% c("Black or African American","Black") ~ "Black",
      race == "White"                                ~ "White",
      TRUE                                           ~ "Other/unknown"
    )
  )

# Genetic similarity (ancestry) and inferred sex from array data
SQL_metadata <- "
  SELECT d_id,
         ancestry_group,
         CAST(enrollment_date AS VARCHAR) AS enrollment_date,
         genetic_sex
  FROM sample_metadata"

sample_metadata <- query(SQL_metadata) %>%
  filter(d_id != "") %>%
  mutate(enrollment_date = as.Date(enrollment_date)) %>%
  unique()

# Merge person with sample metadata; prefer inferred sex, fall back to EHR sex
# When multiple enrollment records exist, use the earliest
person <- person %>%
  inner_join(sample_metadata, by = c("person_source_value" = "d_id")) %>%
  mutate(
    sex = case_when(
      genetic_sex %in% c("F", "M")                  ~ genetic_sex,
      sex_ehr == "Female"                            ~ "F",
      sex_ehr == "Male"                              ~ "M",
      sex_ehr == "Unknown" & genetic_sex == "UNASSIGNED" ~ "U"
    )
  ) %>%
  group_by(person_source_value) %>%
  arrange(enrollment_date) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  filter(!is.na(person_source_value) & person_source_value != "")


###############################################################################
# 3. CLEAN WEIGHT AND HEIGHT DATA
###############################################################################

# Restrict to individuals with genetic data
weight_raw <- weight_raw %>%
  inner_join(person %>% select(person_source_value), by = "person_source_value")

height_raw <- height_raw %>%
  inner_join(person %>% select(person_source_value, birth_date),
             by = "person_source_value") %>%
  mutate(age = floor(as.numeric(measurement_date - birth_date) / 365.25)) %>%
  filter(age >= 18)

# Convert all weights to pounds (US)
# Site-specific heuristics applied for records with missing unit information
# (identified during data QC; see Methods)
weight <- weight_raw %>%
  mutate(weight = case_when(
    unit_concept_name %in% c("ounce (avoirdupois)", "oz")  ~ value_as_number / 16,
    unit_concept_name %in% c("kilogram", "kg")             ~ value_as_number * 2.20462262185,
    unit_concept_name == "pound (US)"                      ~ value_as_number,
    unit_source_value %in% c("lbs", "[lb_av];kg", "LBS")  ~ value_as_number,
    unit_source_value == "kg"                              ~ value_as_number * 2.20462262185,
    unit_source_value == "oz"                              ~ value_as_number / 16,
    # Recover weight from value_source_value when value_as_number is missing
    is.na(value_as_number) & grepl("\\d", value_source_value) ~
      as.numeric(str_extract(value_source_value, "\\d+\\.?\\d*")),
    # Site-specific unit assignments for records with missing unit labels
    # [HELIX INTERNAL: replace site codes below with HRN partner site identifiers]
    unit_concept_name %in% c("", "No matching concept") &
      site %in% c("site_a", "site_b", "site_c", "site_d") &
      unit_source_value == "" & !is.na(value_as_number)        ~ value_as_number / 16,
    unit_concept_name %in% c("", "No matching concept") &
      site == "site_e" & unit_source_value == "" & !is.na(value_as_number) ~
      value_as_number * 2.20462262185,
    unit_concept_name %in% c("", "No matching concept") &
      site == "site_f" & unit_source_value == "" & !is.na(value_as_number) ~
      value_as_number,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(weight) & weight > 10 & weight < 1430)  # plausibility filter

# Convert all heights to inches
height <- height_raw %>%
  mutate(height = case_when(
    unit_concept_name %in% c("inch (US)", "in") ~ value_as_number,
    unit_concept_name %in% c("cm", "centimeter") ~ value_as_number * 0.393701,
    unit_source_value %in% c("in", "Inch(es)")   ~ value_as_number,
    unit_source_value == "cm"                    ~ value_as_number * 0.393701,
    unit_concept_name %in% c("", "No matching concept") &
      unit_source_value == "" & !is.na(value_as_number) ~ value_as_number,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(height) & height > 21 & height < 107)  # plausibility filter (inches)

# Use median height per person for BMI calculation
median_height <- height %>%
  group_by(person_source_value) %>%
  summarise(median_height = median(height), .groups = "drop")

# Compute BMI time series: weight (lbs) and height (in) -> BMI (kg/m2)
bmi <- weight %>%
  select(person_source_value, measurement_date, weight) %>%
  left_join(median_height, by = "person_source_value") %>%
  rename(bmi_date = measurement_date) %>%
  mutate(bmi = (weight / (median_height^2)) * 703) %>%  # standard BMI formula (lbs, inches)
  filter(!is.na(bmi) & bmi > 7 & bmi < 200) %>%
  arrange(person_source_value, bmi_date, bmi) %>%
  group_by(person_source_value, bmi_date) %>%
  filter(row_number() == 1) %>%  # one record per person per date
  ungroup()


###############################################################################
# 4. DRUG EXPOSURE RECORDS
###############################################################################

# Semaglutide (ATC A10BJ06): injectable (Ozempic, Wegovy) and oral (Rybelsus)
semaglutide_raw <- get_drug(atc_level = 5, atc_codes = ATC$semaglutide) %>%
  inner_join(person %>% select(person_source_value), by = "person_source_value")

# Tirzepatide (ATC A10BX16): Mounjaro, Zepbound
tirzepatide_raw <- get_drug(atc_level = 5, atc_codes = ATC$tirzepatide) %>%
  inner_join(person %>% select(person_source_value), by = "person_source_value")

# Other GLP-1 RAs: exenatide, liraglutide, lixisenatide, albiglutide,
#                  dulaglutide, beinaglutide
other_glp1_raw <- get_drug(atc_level = 5, atc_codes = ATC$other_glp1ra) %>%
  inner_join(person %>% select(person_source_value), by = "person_source_value")

# --- Semaglutide brand and dose mapping ---
# Oral doses encoded as fractions so they sort below injectable doses
# (e.g., Rybelsus 3 mg → 0.03, 7 mg → 0.3007, 14 mg → 0.6014)
semaglutide_brand_dose <- semaglutide_raw %>%
  mutate(
    brand_drug = case_when(
      drug_concept_id %in% c(793147, 793152, 37003617, 37003616, 793154,
                              793153, 741832, 780248, 780249, 741834) ~ "Ozempic",
      drug_concept_id %in% c(1537596, 1537597, 1537598, 1537599, 1537600,
                              1537601, 1537602, 1537603, 1537604, 1537605) ~ "Wegovy",
      drug_concept_id %in% c(37496746, 37496842, 37496844, 37496840,
                              37496838, 37496832, 1465681)               ~ "Rybelsus",
      drug_concept_id == 793143                                          ~ "Unknown"
    ),
    dose = case_when(
      # Ozempic doses (mg, injectable)
      drug_concept_id %in% c(793147, 793152, 741832, 741834)            ~ 0.5,
      drug_concept_id %in% c(37003617, 37003616, 793154, 793153)        ~ 1.0,
      drug_concept_id %in% c(780248, 780249)                            ~ 2.0,
      # Wegovy doses (mg, injectable)
      drug_concept_id %in% c(1537596, 1537597)                          ~ 0.25,
      drug_concept_id %in% c(1537598, 1537599)                          ~ 0.5,
      drug_concept_id %in% c(1537600, 1537601)                          ~ 1.0,
      drug_concept_id %in% c(1537602, 1537603)                          ~ 1.7,
      drug_concept_id %in% c(1537604, 1537605)                          ~ 2.4,
      # Rybelsus doses (mg, oral — encoded as fractions for sorting)
      drug_concept_id %in% c(37496842, 37496840)                        ~ 0.03,   # 3 mg
      drug_concept_id %in% c(37496746, 37496844)                        ~ 0.3007, # 7 mg
      drug_concept_id %in% c(37496838, 37496832)                        ~ 0.6014, # 14 mg
      drug_concept_id == 1465681                                         ~ 0.015   # 3 mg (alt concept)
    )
  )

# First semaglutide start date (overall, EHR, and claims separately)
sema_start <- semaglutide_raw %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date = drug_exposure_start_date)

sema_start_ehr <- semaglutide_raw %>%
  filter(data_source == "ehr") %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date_ehr = drug_exposure_start_date)

sema_start_claims <- semaglutide_raw %>%
  filter(data_source == "claims") %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date_claims = drug_exposure_start_date)

# First oral semaglutide (Rybelsus) — used as censor event
first_oral_sema <- semaglutide_brand_dose %>%
  filter(brand_drug == "Rybelsus") %>%
  arrange(person_source_value, drug_exposure_start_date) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(first_oral_semaglutide = drug_exposure_start_date)

# Censor date: last semaglutide prescription + 60 days
sema_last_start <- semaglutide_raw %>%
  arrange(person_source_value, desc(drug_exposure_start_date)) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  mutate(censor_last_drug = drug_exposure_start_date + 60)

# EHR repeat prescription indicator (≥1 refill 7–90 days after first EHR fill)
sema_ehr_repeat <- semaglutide_raw %>%
  filter(data_source == "ehr") %>%
  left_join(sema_start_ehr, by = "person_source_value") %>%
  filter(drug_exposure_start_date >= (drug_start_date_ehr + 7) &
           drug_exposure_start_date <= (drug_start_date_ehr + 90)) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(drug_ehr_repeat = 1) %>%
  select(person_source_value, drug_ehr_repeat)

# Two+ distinct start-date records (used as data quality indicator)
sema_two_records <- semaglutide_raw %>%
  group_by(person_source_value) %>%
  filter(row_number() >= 2) %>%
  distinct(person_source_value) %>%
  mutate(include_2records = 1)

###############################################################################
# Parallel logic for tirzepatide cohort
###############################################################################
# (Identical structure; substituting tirzepatide for semaglutide)

tirz_start <- tirzepatide_raw %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date = drug_exposure_start_date)

tirz_start_ehr <- tirzepatide_raw %>%
  filter(data_source == "ehr") %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date_ehr = drug_exposure_start_date)

tirz_start_claims <- tirzepatide_raw %>%
  filter(data_source == "claims") %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  rename(drug_start_date_claims = drug_exposure_start_date)

tirz_last_start <- tirzepatide_raw %>%
  arrange(person_source_value, desc(drug_exposure_start_date)) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, drug_exposure_start_date) %>%
  mutate(censor_last_drug = drug_exposure_start_date + 60)

tirz_ehr_repeat <- tirzepatide_raw %>%
  filter(data_source == "ehr") %>%
  left_join(tirz_start_ehr, by = "person_source_value") %>%
  filter(drug_exposure_start_date >= (drug_start_date_ehr + 7) &
           drug_exposure_start_date <= (drug_start_date_ehr + 90)) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(drug_ehr_repeat = 1) %>%
  select(person_source_value, drug_ehr_repeat)

tirz_two_records <- tirzepatide_raw %>%
  group_by(person_source_value) %>%
  filter(row_number() >= 2) %>%
  distinct(person_source_value) %>%
  mutate(include_2records = 1)

# Tirzepatide brand and dose mapping
# Identified via ATC code A10BX16 (tirzepatide, injectable, 1.4 mg/dose)
# Doses: 2.5 mg, 5 mg, 7.5 mg, 10 mg, 12.5 mg, 15 mg
# [HELIX INTERNAL: replace the drug_concept_id values below with the correct
#  OMOP concept IDs for each tirzepatide dose, following the same pattern as
#  the semaglutide_brand_dose mapping above]
tirzepatide_brand_dose <- tirzepatide_raw %>%
  mutate(
    brand_drug = "Mounjaro/Zepbound",  # both brands use ATC A10BX16
    dose = case_when(
      drug_concept_id %in% c()  ~ 2.5,   # [HELIX INTERNAL: concept IDs for 2.5 mg]
      drug_concept_id %in% c()  ~ 5.0,   # [HELIX INTERNAL: concept IDs for 5 mg]
      drug_concept_id %in% c()  ~ 7.5,   # [HELIX INTERNAL: concept IDs for 7.5 mg]
      drug_concept_id %in% c()  ~ 10.0,  # [HELIX INTERNAL: concept IDs for 10 mg]
      drug_concept_id %in% c()  ~ 12.5,  # [HELIX INTERNAL: concept IDs for 12.5 mg]
      drug_concept_id %in% c()  ~ 15.0,  # [HELIX INTERNAL: concept IDs for 15 mg]
      TRUE ~ NA_real_
    )
  )


###############################################################################
# 5. PRESCRIPTION PATTERN INDICATORS
###############################################################################

build_glp1_indicators <- function(drug_start, other_drug_start,
                                  other_glp1_raw, index_drug_label) {
  # Prior index drug in 12 months (crossover prior exposure)
  prior_index <- other_drug_start %>%
    inner_join(drug_start %>% rename(index_start = drug_start_date),
               by = "person_source_value") %>%
    filter(drug_start_date < index_start &
             as.numeric(index_start - drug_start_date) <= 365) %>%
    distinct(person_source_value) %>%
    mutate(prior_sema_tirz = 1)

  # Prior other GLP-1 RA in 12 months
  prior_other <- other_glp1_raw %>%
    inner_join(drug_start %>% rename(index_start = drug_start_date),
               by = "person_source_value") %>%
    filter(drug_exposure_start_date < index_start &
             as.numeric(index_start - drug_exposure_start_date) <= 365) %>%
    distinct(person_source_value) %>%
    mutate(prior_other_glp1 = 1)

  # GLP-1 switch within first 3 months (censor event)
  glp1_within_3mo <- rbind(
    other_drug_start %>% rename(drug_exposure_start_date = drug_start_date) %>%
      select(person_source_value, drug_exposure_start_date),
    other_glp1_raw %>% select(person_source_value, drug_exposure_start_date)
  ) %>%
    inner_join(drug_start %>% rename(index_start = drug_start_date),
               by = "person_source_value") %>%
    filter(drug_exposure_start_date >= index_start &
             as.numeric(drug_exposure_start_date - index_start) < 91.25) %>%
    arrange(person_source_value, drug_exposure_start_date) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(censor_glp1_within_3mo = drug_exposure_start_date) %>%
    select(person_source_value, censor_glp1_within_3mo)

  # GLP-1 switch after 3 months (censor event)
  glp1_after_3mo <- rbind(
    other_drug_start %>% rename(drug_exposure_start_date = drug_start_date) %>%
      select(person_source_value, drug_exposure_start_date),
    other_glp1_raw %>% select(person_source_value, drug_exposure_start_date)
  ) %>%
    inner_join(drug_start %>% rename(index_start = drug_start_date),
               by = "person_source_value") %>%
    filter(drug_exposure_start_date >= index_start &
             as.numeric(drug_exposure_start_date - index_start) >= 91.25) %>%
    arrange(person_source_value, drug_exposure_start_date) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(censor_glp1_after_3mo = drug_exposure_start_date) %>%
    select(person_source_value, censor_glp1_after_3mo)

  list(
    prior_index    = prior_index,
    prior_other    = prior_other,
    within_3mo     = glp1_within_3mo,
    after_3mo      = glp1_after_3mo
  )
}

sema_indicators <- build_glp1_indicators(
  drug_start       = sema_start,
  other_drug_start = tirz_start,
  other_glp1_raw   = other_glp1_raw,
  index_drug_label = "semaglutide"
)

tirz_indicators <- build_glp1_indicators(
  drug_start       = tirz_start,
  other_drug_start = sema_start,
  other_glp1_raw   = other_glp1_raw,
  index_drug_label = "tirzepatide"
)


###############################################################################
# 6. COMORBIDITIES
###############################################################################

# --- Type 2 diabetes ---
t2d_condition <- query_condition_chunks(t2diabetes_concept_id_list) %>%
  rename(t2d_date = condition_start_date) %>%
  arrange(person_source_value, t2d_date) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, t2d_date)

# --- Malignancy (condition + observation tables) ---
# Defined as any malignancy diagnosis within 1 year before index date
malignancy_condition <- query_condition_chunks(malignancy_concept_id_list)
malignancy_observation <- query_observation_chunks(malignancy_concept_id_list)

malignancy_all <- rbind(
  malignancy_condition %>%
    select(person_source_value, condition_start_date) %>%
    rename(event_date = condition_start_date),
  malignancy_observation %>%
    select(person_source_value, observation_date) %>%
    rename(event_date = observation_date)
) %>%
  rename(malignancy_date = event_date)

# --- Bariatric surgery (observation + procedure tables) ---
bariatric_observation <- query_observation_chunks(bariatric_surgery_concept_id_list) %>%
  rename(bariatric_date = observation_date) %>%
  arrange(person_source_value, bariatric_date) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, bariatric_date)

bariatric_procedure <- query_procedure(bariatric_surgery_concept_id_list) %>%
  rename(bariatric_date = procedure_date) %>%
  arrange(person_source_value, bariatric_date) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(person_source_value, bariatric_date)

bariatric <- rbind(bariatric_observation, bariatric_procedure) %>%
  arrange(person_source_value, bariatric_date) %>%
  group_by(person_source_value) %>%
  filter(row_number() == 1) %>%
  ungroup()

# --- Pregnancy (condition + observation tables) ---
# Defined as any pregnancy code within 9 months (273 days) before index date
pregnancy_condition <- query_condition_chunks(pregnancy_concept_id_list) %>%
  rename(pregnancy_date = condition_start_date) %>%
  select(person_source_value, pregnancy_date)

pregnancy_observation <- query_observation_chunks(pregnancy_concept_id_list) %>%
  rename(pregnancy_date = observation_date) %>%
  select(person_source_value, pregnancy_date)

pregnancy <- rbind(pregnancy_condition, pregnancy_observation)


###############################################################################
# 7. ASSEMBLE COHORT — FUNCTION (REUSED FOR BOTH DRUGS)
###############################################################################

assemble_cohort <- function(drug_label,
                             drug_start,
                             drug_start_ehr,
                             drug_start_claims,
                             drug_last_start,
                             drug_ehr_repeat,
                             drug_two_records,
                             drug_indicators,
                             first_oral_sema = NULL,  # semaglutide only
                             brand_dose_df) {
  cat("\n===", drug_label, "===\n")

  # --- Latest BMI within 12 months post start (primary eligibility window) ---
  latest_bmi_12mo <- bmi %>%
    inner_join(drug_start, by = "person_source_value") %>%
    filter(bmi_date > drug_start_date &
             bmi_date <= (drug_start_date + 365)) %>%
    arrange(person_source_value, desc(bmi_date)) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(latest_bmi_12mo_date = bmi_date,
           latest_bmi_12mo_days = bmi_date) %>%
    mutate(latest_bmi_12mo_days = as.numeric(latest_bmi_12mo_date - drug_start_date)) %>%
    select(person_source_value, latest_bmi_12mo_date, latest_bmi_12mo_days)

  # --- Malignancy: keep only events within 1 year before index date ---
  malignancy_indexed <- malignancy_all %>%
    inner_join(drug_start, by = "person_source_value") %>%
    filter(malignancy_date >= (drug_start_date - 365.25)) %>%
    arrange(person_source_value, malignancy_date) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(person_source_value, malignancy_date)

  # --- Pregnancy: within 9 months before index date ---
  pregnancy_indexed <- pregnancy %>%
    inner_join(drug_start, by = "person_source_value") %>%
    filter(pregnancy_date >= (drug_start_date - 273)) %>%
    arrange(person_source_value, pregnancy_date) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(person_source_value, pregnancy_date)

  # --- Censor date: earliest of multiple competing events ---
  # Events: last follow-up BMI (in window), GLP-1 switch <3mo, GLP-1 switch >3mo,
  #         last drug Rx +60 days, oral semaglutide (if applicable),
  #         bariatric surgery, pregnancy, malignancy
  censor_components <- latest_bmi_12mo %>%
    left_join(drug_start, by = "person_source_value") %>%
    left_join(drug_start_ehr,    by = "person_source_value") %>%
    left_join(drug_start_claims, by = "person_source_value") %>%
    left_join(drug_last_start %>% select(person_source_value, censor_last_drug),
              by = "person_source_value") %>%
    left_join(drug_ehr_repeat, by = "person_source_value") %>%
    left_join(drug_two_records %>% select(person_source_value, include_2records),
              by = "person_source_value") %>%
    left_join(drug_indicators$within_3mo, by = "person_source_value") %>%
    left_join(drug_indicators$after_3mo,  by = "person_source_value") %>%
    left_join(pregnancy_indexed,  by = "person_source_value") %>%
    left_join(bariatric %>% select(person_source_value, bariatric_date),
              by = "person_source_value") %>%
    left_join(malignancy_indexed, by = "person_source_value")

  if (!is.null(first_oral_sema)) {
    censor_components <- censor_components %>%
      left_join(first_oral_sema, by = "person_source_value")
  } else {
    censor_components$first_oral_semaglutide <- as.Date(NA)
  }

  censor <- censor_components %>%
    mutate(
      include_ehr_records    = ifelse(!is.na(drug_ehr_repeat) & drug_ehr_repeat == 1, 1, 0),
      include_claims_records = ifelse(!is.na(drug_start_date_claims) &
                                        drug_start_date_claims >= drug_start_date &
                                        drug_start_date_claims <= (drug_start_date + 365), 1, 0),
      censor_date = pmin(latest_bmi_12mo_date,
                         censor_glp1_within_3mo,
                         censor_glp1_after_3mo,
                         censor_last_drug,
                         first_oral_semaglutide,
                         bariatric_date,
                         pregnancy_date,
                         malignancy_date,
                         na.rm = TRUE),
      censor_days = as.numeric(censor_date - drug_start_date)
    )

  # --- Baseline BMI: most recent BMI on or before index date, within 6 months ---
  bmi_baseline <- bmi %>%
    inner_join(drug_start, by = "person_source_value") %>%
    filter(bmi_date <= drug_start_date) %>%
    arrange(person_source_value, desc(bmi_date)) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(bmi_baseline = bmi,
           weight_baseline = weight,
           bmi_baseline_date = bmi_date) %>%
    mutate(
      bmi_baseline_days   = as.numeric(bmi_baseline_date - drug_start_date),
      include_bmi_baseline_6mo = ifelse(!is.na(bmi_baseline) &
                                          bmi_baseline_days >= -182.5, 1, NA)
    )

  # --- Last BMI before censor date ---
  last_bmi_before_censor <- censor %>%
    select(person_source_value, censor_date) %>%
    left_join(bmi %>% select(person_source_value, bmi_date),
              by = "person_source_value") %>%
    filter(bmi_date <= censor_date) %>%
    arrange(person_source_value, desc(bmi_date)) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(last_bmi_before_censor_date = bmi_date) %>%
    inner_join(drug_start, by = "person_source_value") %>%
    mutate(last_bmi_before_censor_days =
             as.numeric(last_bmi_before_censor_date - drug_start_date))

  # --- Dose category (highest dose reached through follow-up window) ---
  # Oral doses decoded back to mg (0.015/0.03 -> 3, 0.3007 -> 7, 0.6014 -> 14)
  dose_category_df <- brand_dose_df %>%
    inner_join(last_bmi_before_censor %>%
                 select(person_source_value, last_bmi_before_censor_date),
               by = "person_source_value") %>%
    filter(drug_exposure_start_date < last_bmi_before_censor_date) %>%
    mutate(dose_decoded = case_when(
      dose %in% c(0.015, 0.03) ~ 3,
      dose == 0.3007 ~ 7,
      dose == 0.6014 ~ 14,
      TRUE ~ dose
    )) %>%
    arrange(person_source_value, desc(dose_decoded)) %>%
    group_by(person_source_value) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(dose_category = factor(case_when(
      dose_decoded == 2.4  ~ "2.4 mg",
      dose_decoded == 2.0  ~ "2 mg",
      dose_decoded == 1.7  ~ "1.7 mg",
      dose_decoded == 1.0  ~ "1 mg",
      dose_decoded == 14   ~ "Oral 14 mg",
      dose_decoded %in% c(0.25, 0.5) ~ "0.25 or 0.5 mg",
      dose_decoded %in% c(3, 7)      ~ "Oral 3 or 7 mg",
      dose_decoded %in% c(2.5, 5, 7.5, 10, 12.5, 15) ~
        as.character(dose_decoded),  # tirzepatide doses
      TRUE ~ "Unknown"
    ))) %>%
    select(person_source_value, dose_category, dose_decoded)

  # --- T2D flag (any T2D diagnosis on or before index date) ---
  t2d_flag <- t2d_condition %>%
    inner_join(drug_start, by = "person_source_value") %>%
    filter(t2d_date <= drug_start_date) %>%
    distinct(person_source_value) %>%
    mutate(diabetes_type2 = 1)

  # --- Age at index date ---
  age_df <- person %>%
    select(person_source_value, birth_date, sex, ancestry_group, race_eth) %>%
    inner_join(drug_start, by = "person_source_value") %>%
    mutate(age = floor(as.numeric(drug_start_date - birth_date) / 365.25))

  # --- Inclusion/exclusion criteria ---
  cohort <- drug_start %>%
    left_join(drug_start_ehr,    by = "person_source_value") %>%
    left_join(drug_start_claims, by = "person_source_value") %>%
    left_join(drug_ehr_repeat,   by = "person_source_value") %>%
    left_join(drug_two_records %>% select(person_source_value, include_2records),
              by = "person_source_value") %>%
    left_join(censor %>% select(person_source_value, censor_date, censor_days,
                                 include_ehr_records, include_claims_records),
              by = "person_source_value") %>%
    left_join(bmi_baseline %>% select(person_source_value, bmi_baseline,
                                       weight_baseline, bmi_baseline_date,
                                       bmi_baseline_days, include_bmi_baseline_6mo),
              by = "person_source_value") %>%
    left_join(last_bmi_before_censor %>% select(person_source_value,
                                                  last_bmi_before_censor_date,
                                                  last_bmi_before_censor_days),
              by = "person_source_value") %>%
    left_join(dose_category_df %>% select(person_source_value, dose_category),
              by = "person_source_value") %>%
    left_join(drug_indicators$prior_index %>% select(person_source_value, prior_sema_tirz),
              by = "person_source_value") %>%
    left_join(drug_indicators$prior_other %>% select(person_source_value, prior_other_glp1),
              by = "person_source_value") %>%
    left_join(t2d_flag, by = "person_source_value") %>%
    left_join(age_df %>% select(person_source_value, age, sex,
                                 ancestry_group, race_eth),
              by = "person_source_value") %>%
    left_join(bariatric %>% select(person_source_value, bariatric_date),
              by = "person_source_value") %>%
    left_join(pregnancy_indexed, by = "person_source_value") %>%
    left_join(malignancy_indexed, by = "person_source_value") %>%
    mutate(
      drug              = drug_label,
      prior_sema_tirz   = ifelse(is.na(prior_sema_tirz), 0, prior_sema_tirz),
      prior_other_glp1  = ifelse(is.na(prior_other_glp1), 0, prior_other_glp1),
      diabetes_type2    = ifelse(is.na(diabetes_type2), 0, diabetes_type2),
      followup_duration = censor_days,

      # Inclusion flag (see Methods for full criteria)
      include = as.integer(
        # At least one valid prescription record (EHR or claims)
        (include_ehr_records == 1 | include_claims_records == 1) &
          # Baseline BMI available within 6 months
          !is.na(include_bmi_baseline_6mo) &
          # At least one follow-up BMI (any point post-index, before censor)
          !is.na(last_bmi_before_censor_days) & last_bmi_before_censor_days > 0 &
          # No bariatric surgery before index date
          (is.na(bariatric_date)  | bariatric_date  > drug_start_date) &
          # No pregnancy within 9 months before index date
          (is.na(pregnancy_date) | pregnancy_date  > drug_start_date) &
          # No malignancy within 1 year before index date
          (is.na(malignancy_date) | malignancy_date > drug_start_date) &
          # No oral semaglutide before (or at) index date (sema cohort only)
          (is.na(first_oral_semaglutide) |
             first_oral_semaglutide > drug_start_date)
      )
    )

  cat(sprintf("  Total eligible: %d\n",
              sum(cohort$include == 1, na.rm = TRUE)))
  cohort
}


###############################################################################
# 8. BUILD BOTH COHORTS
###############################################################################

sema_cohort <- assemble_cohort(
  drug_label        = "semaglutide",
  drug_start        = sema_start,
  drug_start_ehr    = sema_start_ehr,
  drug_start_claims = sema_start_claims,
  drug_last_start   = sema_last_start,
  drug_ehr_repeat   = sema_ehr_repeat,
  drug_two_records  = sema_two_records,
  drug_indicators   = sema_indicators,
  first_oral_sema   = first_oral_sema,
  brand_dose_df     = semaglutide_brand_dose
)

tirz_cohort <- assemble_cohort(
  drug_label        = "tirzepatide",
  drug_start        = tirz_start,
  drug_start_ehr    = tirz_start_ehr,
  drug_start_claims = tirz_start_claims,
  drug_last_start   = tirz_last_start,
  drug_ehr_repeat   = tirz_ehr_repeat,
  drug_two_records  = tirz_two_records,
  drug_indicators   = tirz_indicators,
  first_oral_sema   = NULL,
  brand_dose_df     = tirzepatide_brand_dose
)


###############################################################################
# 9. LONGITUDINAL WEIGHT DATA
# Saved separately; used in 03_outcome_definition.R to define the
# 6–12 month minimum weight outcome.
###############################################################################

build_longitudinal <- function(drug_start, cohort_included, drug_label) {
  bmi %>%
    inner_join(
      cohort_included %>%
        filter(include == 1) %>%
        select(person_source_value, drug_start_date, bmi_baseline,
               weight_baseline) %>%
        inner_join(drug_start, by = "person_source_value"),
      by = "person_source_value"
    ) %>%
    mutate(
      days_from_start   = as.numeric(bmi_date - drug_start_date),
      weight_pct_change = ((weight - weight_baseline) / weight_baseline) * 100
    ) %>%
    filter(days_from_start > 0) %>%
    select(person_source_value, bmi_date, days_from_start,
           weight, bmi, weight_pct_change) %>%
    arrange(person_source_value, bmi_date)
}

sema_longitudinal <- build_longitudinal(sema_start, sema_cohort, "semaglutide")
tirz_longitudinal <- build_longitudinal(tirz_start, tirz_cohort, "tirzepatide")


###############################################################################
# NOTE ON INDIVIDUALS APPEARING IN BOTH COHORTS
# Some individuals have records for both semaglutide and tirzepatide. This
# script retains them in both cohorts. Crossover handling is performed in
# 05_primary_analysis.R:
#   - Overall analysis: if a carrier appears in both cohorts, only their
#     first treatment episode is used (earlier drug_start_date).
#   - Drug-specific analyses: each cohort is used independently.
###############################################################################

###############################################################################
# 10. SAVE
###############################################################################

save(sema_cohort, file = PATHS$sema_cohort)
save(tirz_cohort, file = PATHS$tirz_cohort)
save(sema_longitudinal, file = PATHS$sema_longitudinal)
save(tirz_longitudinal, file = PATHS$tirz_longitudinal)

cat("\nSaved:\n",
    PATHS$sema_cohort, "\n",
    PATHS$tirz_cohort, "\n",
    PATHS$sema_longitudinal, "\n",
    PATHS$tirz_longitudinal, "\n")
