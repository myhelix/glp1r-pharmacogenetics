###############################################################################
# 03_outcome_definition.R
#
# Apply cohort inclusion filters, classify GLP1R/GIPR rare variants, and
# define the primary outcome: minimum % weight change from baseline in the
# 6–12 month window post treatment initiation.
#
# Steps:
#   1. Load cohort and longitudinal weight data (from 02_cohort_selection.R)
#   2. Classify GLP1R rare variants (hierarchical: pLoF > damaging missense >
#      rs146868158 > uncertain missense > benign missense)
#   3. Apply cohort inclusion filters (BMI ≥25, ≥6 months follow-up)
#   4. Calculate primary outcome (minimum weight % change, 6–12 months)
#   5. Assemble analysis-ready datasets (overall + drug-specific)
#   6. Save
#
# Outputs:
#   output/analysis_overall.rds       — first treatment per person
#   output/analysis_semaglutide.rds   — semaglutide cohort
#   output/analysis_tirzepatide.rds   — tirzepatide cohort
#   output/variant_details_glp1r.rds  — variant-level annotation table
#   output/variant_details_gipr.rds
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(RNOmni)

source("config.R")

###############################################################################
# 1. LOAD COHORT DATA
###############################################################################

load(PATHS$sema_cohort)       # -> sema_cohort
load(PATHS$tirz_cohort)       # -> tirz_cohort
load(PATHS$sema_longitudinal) # -> sema_longitudinal
load(PATHS$tirz_longitudinal) # -> tirz_longitudinal

# 2-record quality filter: retain only individuals with confirmed repeated
# prescriptions meeting EHR or claims criteria
# [HELIX INTERNAL: these flag files are produced by Helix QC pipelines]
sema_2rec <- read_tsv(PATHS$sema_2record_flag, show_col_types = FALSE)
tirz_2rec <- read_tsv(PATHS$tirz_2record_flag, show_col_types = FALSE)

sema_cohort <- sema_cohort %>%
  inner_join(sema_2rec %>% filter(meets_any_7_90_within_90d_index == 1) %>%
               select(person_source_value), by = "person_source_value")

tirz_cohort <- tirz_cohort %>%
  inner_join(tirz_2rec %>% filter(meets_any_7_90_within_90d_index == 1) %>%
               select(person_source_value), by = "person_source_value")

# Principal components (Europeans; used in regenie sensitivity analysis)
pcs <- read_tsv(PATHS$pcs, show_col_types = FALSE) %>%
  select(person_source_value = FID, PC1:PC10)

# Binary indicator: individual has genetic data from HRN array
has_genetic <- read_tsv(PATHS$has_genetic, show_col_types = FALSE) %>%
  select(person_source_value)


###############################################################################
# 2. VARIANT CLASSIFICATION
#
# Input: carrier TSVs produced by 01_variant_classification.py (Hail).
# Each row is one carrier-variant pair with annotation fields:
#   variant_id, MANE_CT_consequence, REVEL, MANE_CT_am_pathogenicity,
#   MANE_CT_lof
#
# Output: person-level flags saved to output/variant_details_glp1r.rds
#         and output/variant_details_gipr.rds; referenced in all downstream
#         scripts.
###############################################################################

# --- GLP1R carriers (rare variants: gnomAD MAF <0.1%) ---
glp1r_carriers <- read_tsv(PATHS$glp1r_carriers, show_col_types = FALSE) %>%
  rename(person_source_value = d_id) %>%
  mutate(
    REVEL = replace_na(REVEL, 0),
    AM    = replace_na(MANE_CT_am_pathogenicity, 0),
    is_lof = MANE_CT_lof == TRUE |
      MANE_CT_consequence %in% c("stop_gained", "frameshift_variant",
                                 "splice_donor_variant", "splice_acceptor_variant"),
    is_damaging_missense = !is_lof &
      (REVEL >= PARAMS$revel_damaging_threshold | AM >= PARAMS$am_damaging_threshold),
    is_uncertain_missense = !is_lof & !is_damaging_missense &
      ((AM >= PARAMS$am_uncertain_lower & AM < PARAMS$am_damaging_threshold) |
         (REVEL >= PARAMS$revel_uncertain_lower & REVEL < PARAMS$revel_damaging_threshold)),
    is_benign_missense = !is_lof & !is_damaging_missense & !is_uncertain_missense &
      (AM > 0 | REVEL > 0)
  )

# --- rs146868158 (chr6:39085942:C:T) ---
# Queried separately because this variant's gnomAD MAF slightly exceeds the
# 0.1% threshold in Finnish populations; included as a functionally validated
# category (see Methods).
carriers_rs146 <- read_tsv(PATHS$rs146868158_carriers, show_col_types = FALSE) %>%
  rename(person_source_value = d_id) %>%
  select(person_source_value, n_alt_alleles) %>%
  mutate(has_rs146 = 1)

# Collapse GLP1R to person level; assign highest-severity category
# Hierarchy: pLoF > damaging > rs146868158 > uncertain > benign
glp1r_variant_flags <- glp1r_carriers %>%
  group_by(person_source_value) %>%
  summarise(
    has_lof               = as.integer(any(is_lof)),
    has_damaging_missense = as.integer(any(is_damaging_missense)),
    has_uncertain_missense = as.integer(any(is_uncertain_missense)),
    has_benign_missense   = as.integer(any(is_benign_missense)),
    max_REVEL             = max(REVEL, na.rm = TRUE),
    max_AM                = max(AM,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(carriers_rs146 %>% select(person_source_value, has_rs146),
            by = "person_source_value") %>%
  mutate(
    has_rs146 = replace_na(has_rs146, 0),
    # Hierarchical primary classification
    variant_class = case_when(
      has_lof == 1               ~ "pLoF",
      has_damaging_missense == 1 ~ "damaging_missense",
      has_rs146 == 1             ~ "rs146868158",
      has_uncertain_missense == 1 ~ "uncertain_missense",
      has_benign_missense == 1   ~ "benign_missense",
      TRUE                       ~ "non_carrier"
    ),
    carrier_primary      = as.integer(variant_class %in%
      c("pLoF", "damaging_missense", "rs146868158", "uncertain_missense")),
    carrier_lof_damaging = as.integer(variant_class %in% c("pLoF", "damaging_missense")),
    carrier_uncertain    = as.integer(variant_class == "uncertain_missense"),
    carrier_benign       = as.integer(variant_class == "benign_missense"),
    any_glp1r_variant    = as.integer(variant_class != "non_carrier")
  )

cat("GLP1R variant counts:\n")
print(table(glp1r_variant_flags$variant_class))

saveRDS(glp1r_variant_flags,
        file.path(PATHS$output_dir, "variant_details_glp1r.rds"))


# --- GIPR carriers ---
# Four categories: pLoF, damaging missense, uncertain missense, benign missense
# (no functionally validated variant equivalent to rs146868158)
gipr_carriers <- read_tsv(PATHS$gipr_carriers, show_col_types = FALSE) %>%
  rename(person_source_value = d_id) %>%
  mutate(
    REVEL = replace_na(REVEL, 0),
    AM    = replace_na(MANE_CT_am_pathogenicity, 0),
    is_lof = MANE_CT_lof == TRUE |
      MANE_CT_consequence %in% c("stop_gained", "frameshift_variant",
                                 "splice_donor_variant", "splice_acceptor_variant"),
    is_damaging_missense = !is_lof &
      (REVEL >= PARAMS$revel_damaging_threshold | AM >= PARAMS$am_damaging_threshold),
    is_uncertain_missense = !is_lof & !is_damaging_missense &
      ((AM >= PARAMS$am_uncertain_lower & AM < PARAMS$am_damaging_threshold) |
         (REVEL >= PARAMS$revel_uncertain_lower & REVEL < PARAMS$revel_damaging_threshold)),
    is_benign_missense = !is_lof & !is_damaging_missense & !is_uncertain_missense &
      (AM > 0 | REVEL > 0)
  )

gipr_variant_flags <- gipr_carriers %>%
  group_by(person_source_value) %>%
  summarise(
    has_gipr_lof               = as.integer(any(is_lof)),
    has_gipr_damaging_missense = as.integer(any(is_damaging_missense)),
    has_gipr_uncertain_missense = as.integer(any(is_uncertain_missense)),
    has_gipr_benign_missense   = as.integer(any(is_benign_missense)),
    max_REVEL                  = max(REVEL, na.rm = TRUE),
    max_AM                     = max(AM,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    gipr_variant_class = case_when(
      has_gipr_lof == 1                ~ "pLoF",
      has_gipr_damaging_missense == 1  ~ "damaging_missense",
      has_gipr_uncertain_missense == 1 ~ "uncertain_missense",
      has_gipr_benign_missense == 1    ~ "benign_missense",
      TRUE                             ~ "non_carrier"
    ),
    carrier_gipr_primary     = as.integer(gipr_variant_class %in%
                                 c("pLoF", "damaging_missense", "uncertain_missense")),
    carrier_gipr_lof_damaging = as.integer(gipr_variant_class %in%
                                 c("pLoF", "damaging_missense")),
    carrier_gipr_benign      = as.integer(gipr_variant_class == "benign_missense")
  )

cat("GIPR variant counts:\n")
print(table(gipr_variant_flags$gipr_variant_class))

saveRDS(gipr_variant_flags,
        file.path(PATHS$output_dir, "variant_details_gipr.rds"))


###############################################################################
# 3. APPLY COHORT INCLUSION FILTERS
###############################################################################

#' Apply standard inclusion filters and dose categorisation
#' @param cohort  Data frame from 02_cohort_selection.R
#' @param drug    "semaglutide" or "tirzepatide"
prepare_cohort <- function(cohort, drug) {
  df <- cohort %>%
    filter(include == 1) %>%
    # BMI ≥25 kg/m² at baseline (see Methods)
    filter(bmi_baseline >= PARAMS$min_bmi) %>%
    # ≥6 months follow-up (182.5 days)
    filter(followup_duration >= PARAMS$min_followup_days) %>%
    mutate(drug = drug)
    # NOTE: Outlier inspection (weight change beyond the threshold reported in
    # the Methods) is performed manually after matching, within the matched
    # sample only. See 04_matching.R.

  cat(sprintf("  %s after filters: %d\n", drug, nrow(df)))
  df
}

sema_filtered <- prepare_cohort(sema_cohort, "semaglutide")
tirz_filtered <- prepare_cohort(tirz_cohort, "tirzepatide")

# Require genetic data
sema_filtered <- sema_filtered %>%
  inner_join(has_genetic, by = "person_source_value")
tirz_filtered <- tirz_filtered %>%
  inner_join(has_genetic, by = "person_source_value")


###############################################################################
# 4. PRIMARY OUTCOME: MINIMUM WEIGHT % CHANGE AT 6–12 MONTHS
###############################################################################

#' Calculate minimum % weight change in the 6–12 month window
#' (minimum = maximum weight loss, i.e. most negative value)
#' Requires at least one weight measurement between 182.5 and 365 days
#'
#' @param longitudinal  Longitudinal weight data frame (from 02_cohort_selection.R)
#' @param cohort        Filtered cohort data frame
define_outcome_6to12mo <- function(longitudinal, cohort) {
  outcome <- longitudinal %>%
    # Restrict to individuals in the filtered cohort
    inner_join(cohort %>% select(person_source_value),
               by = "person_source_value") %>%
    # 6–12 month window
    filter(days_from_start >= PARAMS$min_followup_days &
             days_from_start <= PARAMS$max_followup_days) %>%
    group_by(person_source_value) %>%
    summarise(
      outcome_weight_pct_change = min(weight_pct_change, na.rm = TRUE),
      n_measurements_6to12mo    = n(),
      .groups = "drop"
    )

  cat(sprintf("  Individuals with ≥1 measurement (6–12 mo): %d\n",
              nrow(outcome)))

  # Merge back; keep only those with outcome data
  cohort %>%
    left_join(outcome, by = "person_source_value") %>%
    filter(!is.na(outcome_weight_pct_change))
}

cat("Semaglutide:\n")
sema_outcome <- define_outcome_6to12mo(sema_longitudinal, sema_filtered)
cat("Tirzepatide:\n")
tirz_outcome <- define_outcome_6to12mo(tirz_longitudinal, tirz_filtered)


###############################################################################
# 5. ASSEMBLE ANALYSIS-READY DATASETS
###############################################################################

add_variant_flags <- function(df) {
  df %>%
    left_join(glp1r_variant_flags, by = "person_source_value") %>%
    left_join(gipr_variant_flags,  by = "person_source_value") %>%
    mutate(across(where(is.integer), ~replace_na(., 0L)),
           variant_class      = replace_na(variant_class,      "non_carrier"),
           gipr_variant_class = replace_na(gipr_variant_class, "non_carrier"),
           carrier_gipr_lof_damaging = replace_na(carrier_gipr_lof_damaging, 0L),
           carrier_gipr_benign       = replace_na(carrier_gipr_benign, 0L))
}

sema_analysis <- sema_outcome %>%
  left_join(pcs, by = "person_source_value") %>%
  add_variant_flags() %>%
  mutate(ancestry_group = factor(ancestry_group,
                                     levels = c("European","African","Americas",
                                                "EastAsian","SouthAsian","Other")))

tirz_analysis <- tirz_outcome %>%
  left_join(pcs, by = "person_source_value") %>%
  add_variant_flags() %>%
  mutate(ancestry_group = factor(ancestry_group,
                                     levels = c("European","African","Americas",
                                                "EastAsian","SouthAsian","Other")))

# Overall dataset: when an individual appears in both cohorts, use only their
# first treatment episode (earlier drug_start_date)
overall_analysis <- bind_rows(sema_analysis, tirz_analysis) %>%
  group_by(person_source_value) %>%
  arrange(drug_start_date) %>%
  slice(1) %>%
  ungroup()

cat(sprintf(
  "\nFinal cohort sizes:\n  Overall:      %d (carriers: %d)\n  Semaglutide:  %d\n  Tirzepatide:  %d\n",
  nrow(overall_analysis),  sum(overall_analysis$carrier_primary, na.rm = TRUE),
  nrow(sema_analysis),
  nrow(tirz_analysis)
))


###############################################################################
# 6. SAVE
###############################################################################

saveRDS(overall_analysis, file.path(PATHS$output_dir, "analysis_overall.rds"))
saveRDS(sema_analysis,    file.path(PATHS$output_dir, "analysis_semaglutide.rds"))
saveRDS(tirz_analysis,    file.path(PATHS$output_dir, "analysis_tirzepatide.rds"))

cat("\nSaved analysis datasets to", PATHS$output_dir, "\n")
