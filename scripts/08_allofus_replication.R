###############################################################################
# 08_allofus_replication.R
#
# Replication of primary GLP1R and GIPR findings in the All of Us (AoU)
# Research Program.
#
# IMPORTANT: This script runs on the All of Us Researcher Workbench
# (cloud-based; Google Cloud / BigQuery). Paths and dataset identifiers
# marked [AOU INTERNAL] require an approved AoU data access tier and cannot
# be executed externally. See Methods for access details.
#
# Data:
#   - Phenotype / EHR: All of Us OMOP CDM (BigQuery; accessed via bigrquery)
#   - Genomic: All of Us short-read WGS (GRCh38), extracted via Hail on
#     Dataproc; see note in Section 2 below.
#
# Steps:
#   1. Connect to AoU BigQuery CDR; define concept IDs
#   2. Load GLP1R / GIPR carrier data (Hail WGS extraction)
#   3. Query semaglutide / tirzepatide drug exposure
#   4. Query longitudinal weight measurements
#   5. Apply cohort inclusion criteria (same as primary HRN analysis)
#   6. Define primary outcome (minimum % weight change, 6–12 months)
#   7. Assemble analysis-ready dataset; variant classification
#   8. Perform 1:10 matching (same specification as 04_matching.R)
#   9. Run weighted regression + drug × carrier interaction test
#  10. Save results
#
# Outputs (saved to output/allofus/):
#   aou_matched_glp1r.rds       — matched GLP1R dataset
#   aou_matched_gipr.rds        — matched GIPR dataset
#   aou_balance_summary.rds     — covariate balance SMDs
#   aou_results_primary.rds     — regression results
#   aou_results_interaction.rds — drug × carrier interaction
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(bigrquery)
library(DBI)
library(MatchIt)
library(sandwich)
library(lmtest)

source("config.R")

aou_dir <- file.path(PATHS$output_dir, "allofus")
dir.create(aou_dir, showWarnings = FALSE)


###############################################################################
# 1. BIGQUERY CONNECTION AND CONCEPT IDS
###############################################################################

# All of Us CDR dataset (set by the Researcher Workbench environment)
# [AOU INTERNAL: WORKSPACE_CDR and GOOGLE_PROJECT are environment variables
#  populated automatically in the AoU Researcher Workbench]
CDR     <- Sys.getenv("WORKSPACE_CDR")
PROJECT <- Sys.getenv("GOOGLE_PROJECT")

con <- dbConnect(
  bigrquery::bigquery(),
  project = PROJECT,
  dataset = CDR,
  billing = PROJECT
)

# Helper: run BigQuery SQL and return a tibble
bq_query <- function(sql) {
  dbGetQuery(con, sql) %>% as_tibble()
}

# --- Concept IDs ---
# Standard OMOP concept IDs verified in the AoU CDR for this analysis.
# [AOU INTERNAL: confirm concept IDs against the AoU CDR vocabulary tables
#  before use — vocabulary versions may differ across CDR releases]

# Drug ingredient concept IDs (used to identify GLP-1 RA prescriptions)
CONCEPT <- list(
  # GLP-1 RA ingredients
  semaglutide = "[AOU INTERNAL: semaglutide ingredient concept_id]",
  tirzepatide = "[AOU INTERNAL: tirzepatide ingredient concept_id]",

  # Weight measurements
  body_weight_kg  = "[AOU INTERNAL: body weight (kg) concept_id]",
  body_weight_lbs = "[AOU INTERNAL: body weight (lbs) concept_id]",

  # Height
  body_height_cm  = "[AOU INTERNAL: body height (cm) concept_id]",
  body_height_in  = "[AOU INTERNAL: body height (inches) concept_id]",

  # BMI
  bmi             = "[AOU INTERNAL: BMI concept_id]",

  # Conditions
  diabetes_type2  = "[AOU INTERNAL: type 2 diabetes mellitus concept_id or descendant set]",
  bariatric_surgery = "[AOU INTERNAL: bariatric surgery concept_id or descendant set]",
  malignancy        = "[AOU INTERNAL: malignancy concept_id or descendant set]",
  pregnancy         = "[AOU INTERNAL: pregnancy concept_id or descendant set]"
)


###############################################################################
# 2. GLP1R / GIPR CARRIER DATA FROM HAIL WGS EXTRACTION
#
# Rare coding variant carriers were extracted from the All of Us short-read
# WGS data using a Hail script analogous to 01_variant_classification.py,
# applied to the AoU WGS VDS (Variant Dataset) on a Dataproc cluster.
#
# Filtering criteria (identical to HRN primary analysis):
#   - MANE Select transcript coding variants
#   - MAF < 0.1% in gnomAD v4.1 (each ancestry group)
#   - Unrelated individuals per AoU relatedness estimates
#
# The Hail extraction script is provided as a code comment block below.
# Run it on a Dataproc cluster before executing this R script.
#
# --- BEGIN HAIL EXTRACTION (Python; run on AoU Dataproc cluster) ---
#
# import hail as hl
# import pandas as pd
#
# hl.init(default_reference='GRCh38',
#         tmp_dir='[AOU INTERNAL: Hail temporary directory on GCS]')
#
# # AoU WGS VDS path
# # [AOU INTERNAL: path to AoU WGS VDS on GCS — provided via Researcher Workbench]
# vds = hl.vds.read_vds('[AOU INTERNAL: path to AoU WGS VDS]')
# mt  = hl.vds.to_dense_mt(vds)
#
# # Filter to unrelated individuals
# # [AOU INTERNAL: unrelated sample list from AoU relatedness pipeline]
# unrelated = hl.import_table('[AOU INTERNAL: path to unrelated sample list]',
#                              key='person_id')
# mt = mt.filter_cols(hl.is_defined(unrelated[mt.col_key[0]]))
#
# # Apply the same extraction logic as scripts/01_variant_classification.py:
# # filter to gene, construct variant_id, filter n_alt_alleles > 0, export.
# # Annotation fields: REVEL, AlphaMissense (MANE_CT_am_pathogenicity),
# # LOFTEE (MANE_CT_lof), VEP consequence (MANE_CT_consequence).
# #
# # [AOU INTERNAL: annotation fields from AoU VEP annotation tables may use
# #  different column names — verify against the AoU WGS annotation schema]
#
# for gene in ['GLP1R', 'GIPR']:
#     mt_gene = mt.filter_rows(mt.gene_symbol == gene)
#     mt_gene = mt_gene.filter_rows(mt_gene.gnomad_maf < 0.001)
#     mt_gene = mt_gene.annotate_entries(
#         n_alt_alleles = mt_gene.GT.n_alt_alleles()
#     )
#     mt_gene = mt_gene.filter_entries(mt_gene.n_alt_alleles > 0)
#     entries = mt_gene.entries()
#     entries.export(f'[AOU INTERNAL: GCS output path]/{gene.lower()}_aou_carriers.tsv')
#
# --- END HAIL EXTRACTION ---

# Load carrier TSVs written by the Hail extraction above
# [AOU INTERNAL: adjust paths to GCS or local paths after localizing files]
aou_glp1r_carriers <- read_tsv(
  "[AOU INTERNAL: path to glp1r_aou_carriers.tsv]",
  show_col_types = FALSE
) %>% rename(person_id = person_id)   # AoU uses integer person_id

aou_gipr_carriers <- read_tsv(
  "[AOU INTERNAL: path to gipr_aou_carriers.tsv]",
  show_col_types = FALSE
)

# AoU uses integer person_id; convert to character for join consistency
aou_glp1r_carriers <- aou_glp1r_carriers %>%
  mutate(person_id = as.character(person_id))
aou_gipr_carriers  <- aou_gipr_carriers  %>%
  mutate(person_id = as.character(person_id))


###############################################################################
# 3. VARIANT CLASSIFICATION
# Same logic as 03_outcome_definition.R; thresholds from config.R
###############################################################################

classify_glp1r <- function(carriers, rs146_carriers) {
  variant_flags <- carriers %>%
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
    ) %>%
    group_by(person_id) %>%
    summarise(
      has_lof               = as.integer(any(is_lof)),
      has_damaging_missense = as.integer(any(is_damaging_missense)),
      has_uncertain_missense = as.integer(any(is_uncertain_missense)),
      has_benign_missense   = as.integer(any(is_benign_missense)),
      max_REVEL             = max(REVEL, na.rm = TRUE),
      max_AM                = max(AM,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(rs146_carriers %>% select(person_id, has_rs146),
              by = "person_id") %>%
    mutate(
      has_rs146    = replace_na(has_rs146, 0L),
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

  cat("AoU GLP1R variant counts:\n")
  print(table(variant_flags$variant_class))
  variant_flags
}

classify_gipr <- function(carriers) {
  carriers %>%
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
           (REVEL >= PARAMS$revel_uncertain_lower & REVEL < PARAMS$revel_damaging_threshold))
    ) %>%
    group_by(person_id) %>%
    summarise(
      has_gipr_lof               = as.integer(any(is_lof)),
      has_gipr_damaging_missense = as.integer(any(is_damaging_missense)),
      has_gipr_uncertain_missense = as.integer(any(is_uncertain_missense)),
      max_REVEL                  = max(REVEL, na.rm = TRUE),
      max_AM                     = max(AM,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      gipr_variant_class   = case_when(
        has_gipr_lof == 1                ~ "pLoF",
        has_gipr_damaging_missense == 1  ~ "damaging_missense",
        has_gipr_uncertain_missense == 1 ~ "uncertain_missense",
        TRUE                             ~ "non_carrier"
      ),
      carrier_gipr_primary = as.integer(gipr_variant_class != "non_carrier")
    )
}

# rs146868158 — queried separately (see Methods; MAF slightly exceeds 0.1% threshold in Finnish)
# [AOU INTERNAL: rs146868158 carrier list extracted from AoU WGS or genotype data]
aou_rs146_carriers <- read_tsv(
  "[AOU INTERNAL: path to rs146868158_aou_carriers.tsv]",
  show_col_types = FALSE
) %>%
  mutate(person_id = as.character(person_id), has_rs146 = 1L)

aou_glp1r_flags <- classify_glp1r(aou_glp1r_carriers, aou_rs146_carriers)
aou_gipr_flags  <- classify_gipr(aou_gipr_carriers)


###############################################################################
# 4. DRUG EXPOSURE QUERY
###############################################################################

# Query semaglutide drug exposures from AoU OMOP CDM
# drug_exposure rows are deduplicated to one row per person per drug_start_date
sema_sql <- sprintf("
  SELECT
    CAST(de.person_id AS STRING)                         AS person_id,
    de.drug_concept_id,
    c.concept_name                                       AS drug_name,
    MIN(de.drug_exposure_start_date)                     AS drug_start_date,
    MAX(COALESCE(de.drug_exposure_end_date,
         DATE_ADD(de.drug_exposure_start_date,
                  INTERVAL de.days_supply DAY),
         DATE_ADD(de.drug_exposure_start_date, INTERVAL 30 DAY)))
                                                         AS drug_end_date
  FROM `%s.drug_exposure` de
  JOIN `%s.concept_ancestor` ca
    ON de.drug_concept_id = ca.descendant_concept_id
  JOIN `%s.concept` c
    ON de.drug_concept_id = c.concept_id
  WHERE ca.ancestor_concept_id = %s   -- semaglutide ingredient
  GROUP BY de.person_id, de.drug_concept_id, c.concept_name
",
CDR, CDR, CDR, CONCEPT$semaglutide)

tirz_sql <- sprintf("
  SELECT
    CAST(de.person_id AS STRING)                         AS person_id,
    de.drug_concept_id,
    c.concept_name                                       AS drug_name,
    MIN(de.drug_exposure_start_date)                     AS drug_start_date,
    MAX(COALESCE(de.drug_exposure_end_date,
         DATE_ADD(de.drug_exposure_start_date,
                  INTERVAL de.days_supply DAY),
         DATE_ADD(de.drug_exposure_start_date, INTERVAL 30 DAY)))
                                                         AS drug_end_date
  FROM `%s.drug_exposure` de
  JOIN `%s.concept_ancestor` ca
    ON de.drug_concept_id = ca.descendant_concept_id
  JOIN `%s.concept` c
    ON de.drug_concept_id = c.concept_id
  WHERE ca.ancestor_concept_id = %s   -- tirzepatide ingredient
  GROUP BY de.person_id, de.drug_concept_id, c.concept_name
",
CDR, CDR, CDR, CONCEPT$tirzepatide)

cat("Querying semaglutide drug exposures...\n")
aou_sema_rx <- bq_query(sema_sql) %>%
  mutate(drug = "semaglutide", drug_start_date = as.Date(drug_start_date))

cat("Querying tirzepatide drug exposures...\n")
aou_tirz_rx <- bq_query(tirz_sql) %>%
  mutate(drug = "tirzepatide", drug_start_date = as.Date(drug_start_date))

cat(sprintf("Drug exposures: sema=%d rows, tirz=%d rows\n",
            nrow(aou_sema_rx), nrow(aou_tirz_rx)))


###############################################################################
# 5. WEIGHT AND HEIGHT QUERIES
###############################################################################

weight_sql <- sprintf("
  SELECT
    CAST(m.person_id AS STRING)     AS person_id,
    m.measurement_date              AS measurement_date,
    m.value_as_number               AS value_raw,
    m.unit_concept_id,
    uc.concept_name                 AS unit_name,
    m.measurement_concept_id
  FROM `%s.measurement` m
  LEFT JOIN `%s.concept` uc ON m.unit_concept_id = uc.concept_id
  WHERE m.measurement_concept_id IN (%s, %s)
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0
",
CDR, CDR, CONCEPT$body_weight_kg, CONCEPT$body_weight_lbs)

height_sql <- sprintf("
  SELECT
    CAST(m.person_id AS STRING)     AS person_id,
    m.measurement_date              AS measurement_date,
    m.value_as_number               AS value_raw,
    m.measurement_concept_id,
    uc.concept_name                 AS unit_name
  FROM `%s.measurement` m
  LEFT JOIN `%s.concept` uc ON m.unit_concept_id = uc.concept_id
  WHERE m.measurement_concept_id IN (%s, %s)
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0
",
CDR, CDR, CONCEPT$body_height_cm, CONCEPT$body_height_in)

cat("Querying weight measurements...\n")
aou_weights_raw <- bq_query(weight_sql) %>%
  mutate(measurement_date = as.Date(measurement_date))

cat("Querying height measurements...\n")
aou_heights_raw <- bq_query(height_sql) %>%
  mutate(measurement_date = as.Date(measurement_date))

# Convert all weights to kg
# [AOU INTERNAL: confirm unit_concept_id values for kg vs lbs in the CDR]
aou_weights <- aou_weights_raw %>%
  mutate(
    weight_kg = case_when(
      measurement_concept_id == as.integer(CONCEPT$body_weight_lbs) ~ value_raw * 0.453592,
      TRUE                                                           ~ value_raw
    )
  ) %>%
  filter(weight_kg >= 30, weight_kg <= 400) %>%   # plausible range filter
  select(person_id, measurement_date, weight_kg)

# Convert heights to cm; take modal height per person (stable trait)
aou_heights <- aou_heights_raw %>%
  mutate(
    height_cm = case_when(
      measurement_concept_id == as.integer(CONCEPT$body_height_in) ~ value_raw * 2.54,
      TRUE                                                          ~ value_raw
    )
  ) %>%
  filter(height_cm >= 100, height_cm <= 250) %>%
  group_by(person_id) %>%
  summarise(height_cm = median(height_cm, na.rm = TRUE), .groups = "drop")


###############################################################################
# 6. ASSEMBLE COHORTS
###############################################################################

# Person-level demographics from AoU
demo_sql <- sprintf("
  SELECT
    CAST(p.person_id AS STRING)     AS person_id,
    p.birth_datetime,
    CASE p.sex_at_birth_concept_id
      WHEN 45878463 THEN 'M'
      WHEN 45880669 THEN 'F'
      ELSE 'Unknown'
    END                             AS sex,
    -- AoU uses HARE for ancestry; map to match primary analysis categories
    -- [AOU INTERNAL: ancestry grouping from AoU HARE or PC-based assignment]
    '[AOU INTERNAL: ancestry column]'  AS ancestry_group
  FROM `%s.person` p
",
CDR)

cat("Querying demographics...\n")
aou_demo <- bq_query(demo_sql) %>%
  mutate(birth_datetime = as.Date(birth_datetime))

# T2D flag
t2d_sql <- sprintf("
  SELECT DISTINCT CAST(person_id AS STRING) AS person_id, 1 AS diabetes_type2
  FROM `%s.condition_occurrence`
  WHERE condition_concept_id IN
    (SELECT descendant_concept_id FROM `%s.concept_ancestor`
     WHERE ancestor_concept_id = %s)
",
CDR, CDR, CONCEPT$diabetes_type2)

cat("Querying T2D...\n")
aou_t2d <- bq_query(t2d_sql)

# Bariatric exclusion flag
bariatric_sql <- sprintf("
  SELECT DISTINCT CAST(person_id AS STRING) AS person_id, 1 AS had_bariatric
  FROM `%s.procedure_occurrence`
  WHERE procedure_concept_id IN
    (SELECT descendant_concept_id FROM `%s.concept_ancestor`
     WHERE ancestor_concept_id = %s)
",
CDR, CDR, CONCEPT$bariatric_surgery)
aou_bariatric <- bq_query(bariatric_sql)

#' Build cohort for one drug
#' Returns one row per person per treatment episode; applies same inclusion
#' criteria as 02_cohort_selection.R / 03_outcome_definition.R
build_aou_cohort <- function(rx_df, drug_label) {
  # Take earliest drug start per person
  rx_first <- rx_df %>%
    group_by(person_id) %>%
    arrange(drug_start_date) %>%
    slice(1) %>%
    ungroup()

  # Baseline weight: closest measurement before drug_start (up to 180 days prior)
  wt_baseline <- aou_weights %>%
    inner_join(rx_first %>% select(person_id, drug_start_date), by = "person_id") %>%
    filter(measurement_date >= drug_start_date - 180,
           measurement_date <= drug_start_date) %>%
    group_by(person_id) %>%
    arrange(desc(measurement_date)) %>%
    slice(1) %>%
    ungroup() %>%
    rename(weight_baseline_kg = weight_kg, baseline_date = measurement_date)

  # Age at drug start
  cohort <- rx_first %>%
    inner_join(aou_demo,    by = "person_id") %>%
    inner_join(wt_baseline, by = c("person_id", "drug_start_date")) %>%
    inner_join(aou_heights, by = "person_id") %>%
    left_join(aou_t2d,      by = "person_id") %>%
    left_join(aou_bariatric,by = "person_id") %>%
    mutate(
      drug            = drug_label,
      age             = as.numeric(drug_start_date - birth_datetime) / 365.25,
      bmi_baseline    = weight_baseline_kg / (height_cm / 100)^2,
      diabetes_type2  = replace_na(diabetes_type2, 0L),
      had_bariatric   = replace_na(had_bariatric,  0L)
    ) %>%
    # Inclusion filters (same as primary analysis)
    filter(bmi_baseline >= PARAMS$min_bmi) %>%
    filter(had_bariatric == 0) %>%
    select(person_id, drug, drug_start_date, age, sex, bmi_baseline,
           weight_baseline_kg, height_cm, diabetes_type2, ancestry_group)

  cat(sprintf("  AoU %s after inclusion filters: N=%d\n", drug_label, nrow(cohort)))
  cohort
}

cat("\nBuilding AoU cohorts...\n")
aou_sema_cohort <- build_aou_cohort(aou_sema_rx, "semaglutide")
aou_tirz_cohort <- build_aou_cohort(aou_tirz_rx, "tirzepatide")


###############################################################################
# 7. PRIMARY OUTCOME: MINIMUM % WEIGHT CHANGE AT 6–12 MONTHS
###############################################################################

define_aou_outcome <- function(cohort) {
  aou_weights %>%
    inner_join(cohort %>% select(person_id, drug_start_date, weight_baseline_kg),
               by = "person_id") %>%
    mutate(days_from_start = as.numeric(measurement_date - drug_start_date),
           weight_pct_change = (weight_kg - weight_baseline_kg) / weight_baseline_kg * 100) %>%
    filter(days_from_start >= PARAMS$min_followup_days,
           days_from_start <= PARAMS$max_followup_days) %>%
    group_by(person_id) %>%
    summarise(
      outcome_weight_pct_change = min(weight_pct_change, na.rm = TRUE),
      n_measurements_6to12mo    = n(),
      followup_duration         = max(days_from_start),
      .groups = "drop"
    ) %>%
    inner_join(cohort, by = "person_id") %>%
    filter(!is.na(outcome_weight_pct_change))
}

cat("\nDefining outcomes...\n")
aou_sema_outcome <- define_aou_outcome(aou_sema_cohort)
aou_tirz_outcome <- define_aou_outcome(aou_tirz_cohort)

cat(sprintf("  With outcome: sema=%d, tirz=%d\n",
            nrow(aou_sema_outcome), nrow(aou_tirz_outcome)))

# Add variant flags; overall = first treatment per person
add_aou_variant_flags <- function(df) {
  df %>%
    mutate(person_id = as.character(person_id)) %>%
    left_join(aou_glp1r_flags, by = "person_id") %>%
    left_join(aou_gipr_flags,  by = "person_id") %>%
    mutate(
      across(c(carrier_primary, carrier_lof_damaging, carrier_uncertain,
               carrier_benign, any_glp1r_variant, carrier_gipr_primary),
             ~replace_na(., 0L)),
      variant_class      = replace_na(variant_class,      "non_carrier"),
      gipr_variant_class = replace_na(gipr_variant_class, "non_carrier"),
      # Placeholder for prior GLP-1 RA use
      # [AOU INTERNAL: compute prior_sema_tirz and prior_other_glp1 from
      #  drug_exposure table using same logic as 02_cohort_selection.R]
      prior_sema_tirz  = 0L,
      prior_other_glp1 = 0L,
      dose_category    = "standard",
      ancestry_group = factor(ancestry_group,
                                  levels = c("European","African","Americas",
                                             "EastAsian","SouthAsian","Other"))
    )
}

aou_sema_analysis <- add_aou_variant_flags(aou_sema_outcome)
aou_tirz_analysis <- add_aou_variant_flags(aou_tirz_outcome)

aou_overall <- bind_rows(aou_sema_analysis, aou_tirz_analysis) %>%
  group_by(person_id) %>%
  arrange(drug_start_date) %>%
  slice(1) %>%
  ungroup()

cat(sprintf(
  "\nAoU final cohort: overall N=%d (GLP1R carriers=%d; GIPR carriers=%d)\n",
  nrow(aou_overall),
  sum(aou_overall$carrier_primary,      na.rm = TRUE),
  sum(aou_overall$carrier_gipr_primary, na.rm = TRUE)
))


###############################################################################
# 8. MATCHING (same specification as 04_matching.R)
###############################################################################

perform_aou_matching <- function(data, carrier_var, label,
                                 ratio = PARAMS$match_ratio) {
  cat(sprintf("\n--- %s ---\n", label))
  n_carriers <- sum(data[[carrier_var]] == 1, na.rm = TRUE)
  cat(sprintf("  Input: N=%d | carriers=%d\n", nrow(data), n_carriers))
  if (n_carriers < 1) { cat("  Skipped: no carriers\n"); return(NULL) }

  exact_vars <- PARAMS$match_exact
  exact_vars <- exact_vars[exact_vars %in% names(data)]
  exact_vars <- exact_vars[sapply(exact_vars, function(v) {
    length(unique(na.omit(data[[v]]))) >= 2 &&
      length(unique(na.omit(data[[v]][data[[carrier_var]] == 1]))) >= 1
  })]
  mahal_vars <- PARAMS$match_mahal[PARAMS$match_mahal %in% names(data)]
  all_vars   <- unique(c(exact_vars, mahal_vars))

  match_formula <- as.formula(paste(carrier_var, "~", paste(all_vars, collapse = " + ")))
  exact_formula <- if (length(exact_vars) > 0)
    as.formula(paste("~", paste(exact_vars, collapse = " + "))) else NULL

  calipers <- c(bmi_baseline = PARAMS$caliper_bmi_kgm2,
                followup_duration = PARAMS$caliper_fu_days)
  calipers <- calipers[names(calipers) %in% names(data)]

  matched_obj <- tryCatch(
    matchit(match_formula, data = data,
            method = PARAMS$match_method, distance = PARAMS$match_distance,
            exact = exact_formula, caliper = calipers, std.caliper = FALSE,
            ratio = ratio, replace = FALSE, estimand = "ATT"),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
  if (is.null(matched_obj)) return(NULL)

  matched_df <- match.data(matched_obj) %>%
    group_by(subclass) %>%
    mutate(
      n_controls_in_set = sum(.data[[carrier_var]] == 0),
      weights_att = if_else(.data[[carrier_var]] == 1, 1, 1 / n_controls_in_set)
    ) %>%
    ungroup()

  n_c <- sum(matched_df[[carrier_var]] == 1)
  n_k <- sum(matched_df[[carrier_var]] == 0)
  cat(sprintf("  Matched: %d carriers, %d controls (ratio %.1f:1)\n",
              n_c, n_k, n_k / n_c))
  list(matched_data = matched_df, matchit_object = matched_obj,
       n_carriers = n_c, n_controls = n_k, carrier_var = carrier_var)
}

no_variant_aou <- aou_overall %>%
  filter(any_glp1r_variant == 0 & carrier_gipr_primary == 0)

aou_match_glp1r <- perform_aou_matching(
  bind_rows(aou_overall %>% filter(carrier_primary == 1), no_variant_aou) %>%
    mutate(ancestry_group = droplevels(ancestry_group)),
  "carrier_primary", "AoU GLP1R primary — overall"
)

aou_match_gipr <- perform_aou_matching(
  bind_rows(aou_overall %>% filter(carrier_gipr_primary == 1), no_variant_aou) %>%
    mutate(ancestry_group = droplevels(ancestry_group)),
  "carrier_gipr_primary", "AoU GIPR primary — overall"
)


###############################################################################
# 9. REGRESSION
###############################################################################

adj_vars <- c("age", "bmi_baseline", "followup_duration")

run_aou_regression <- function(df, carrier_var, label) {
  if (is.null(df)) return(NULL)
  avars <- adj_vars[adj_vars %in% names(df)]
  adj_str <- if (length(avars) > 0) paste("+", paste(avars, collapse = " + ")) else ""
  formula_obj <- as.formula(paste("outcome_weight_pct_change ~", carrier_var, adj_str))

  fit <- lm(formula_obj, data = df, weights = df$weights_att)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = PARAMS$vcov_type))
  if (!carrier_var %in% rownames(ct)) return(NULL)

  coef_row <- ct[carrier_var, ]
  est <- coef_row["Estimate"]; se <- coef_row["Std. Error"]; pval <- coef_row["Pr(>|t|)"]
  carriers <- df[df[[carrier_var]] == 1, ]
  controls <- df[df[[carrier_var]] == 0, ]

  tibble(
    analysis     = label,
    n_carriers   = nrow(carriers),
    n_controls   = nrow(controls),
    mean_carrier = round(weighted.mean(carriers$outcome_weight_pct_change,
                                       carriers$weights_att, na.rm = TRUE), 2),
    mean_control = round(weighted.mean(controls$outcome_weight_pct_change,
                                       controls$weights_att, na.rm = TRUE), 2),
    estimate     = round(est, 2),
    se           = round(se,  3),
    ci_lo        = round(est - 1.96 * se, 2),
    ci_hi        = round(est + 1.96 * se, 2),
    p_value      = pval
  )
}

aou_results <- bind_rows(
  run_aou_regression(aou_match_glp1r$matched_data, "carrier_primary",
                     "AoU GLP1R primary — overall"),
  run_aou_regression(
    aou_match_glp1r$matched_data %>% filter(drug == "semaglutide"),
    "carrier_primary", "AoU GLP1R primary — semaglutide"),
  run_aou_regression(
    aou_match_glp1r$matched_data %>% filter(drug == "tirzepatide"),
    "carrier_primary", "AoU GLP1R primary — tirzepatide"),
  run_aou_regression(aou_match_gipr$matched_data, "carrier_gipr_primary",
                     "AoU GIPR primary — overall")
)

cat("\nAoU replication results:\n")
print(
  aou_results %>%
    mutate(p_fmt = formatC(p_value, digits = 3, format = "g"),
           ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)) %>%
    select(analysis, n_carriers, n_controls, mean_carrier, mean_control,
           estimate, ci, p_fmt),
  n = Inf, width = 120
)


###############################################################################
# 10. SAVE
###############################################################################

if (!is.null(aou_match_glp1r))
  saveRDS(aou_match_glp1r$matched_data, file.path(aou_dir, "aou_matched_glp1r.rds"))
if (!is.null(aou_match_gipr))
  saveRDS(aou_match_gipr$matched_data,  file.path(aou_dir, "aou_matched_gipr.rds"))
saveRDS(aou_results, file.path(aou_dir, "aou_results_primary.rds"))

cat("\nAoU replication results saved to", aou_dir, "\n")
