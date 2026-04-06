###############################################################################
# config.R
# Configuration file for GLP1R pharmacogenetics analysis
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
#
# USAGE: source("config.R") at the top of each analysis script.
#
# Pat. Pending, US Pat App. 19/002,020.
#
# IMPORTANT: This code was developed on Helix Research Network (HRN)
# infrastructure and the All of Us Research Program. Paths marked
# [HELIX INTERNAL] and [AOU INTERNAL] require access to those respective
# environments and cannot be run externally.
###############################################################################

###############################################################################
# INPUT DATA PATHS
###############################################################################

# --- Genetic data (from Hail extraction; see scripts/01_variant_classification.py) ---
PATHS <- list(

  # GLP1R carrier file — one row per carrier-variant pair
  # (output of 01_variant_classification.py; variants with gnomAD MAF <0.1%)
  glp1r_carriers = "output/glp1r_rare_coding_carriers.tsv",

  # GIPR carrier file
  gipr_carriers  = "output/gipr_rare_coding_carriers.tsv",

  # rs146868158 (chr6:39085942:C:T) carrier file — queried separately because
  # this variant's gnomAD MAF slightly exceeds the 0.1% threshold in Finnish
  # populations; included as a functionally validated category (see Methods).
  # One row per carrier.
  # [HELIX INTERNAL: extracted from HRN genotype data]
  rs146868158_carriers = "output/rs146868158_carriers.tsv",

  # Phenotype data (output of 02_cohort_selection.R)
  sema_cohort       = "output/semaglutide_cohort.RData",
  tirz_cohort       = "output/tirzepatide_cohort.RData",
  sema_longitudinal = "output/semaglutide_longitudinal_weights.RData",
  tirz_longitudinal = "output/tirzepatide_longitudinal_weights.RData",

  # Duplicate-record flag files (from internal QC; see Methods)
  # [HELIX INTERNAL: produced by Helix QC pipelines]
  sema_2record_flag = "output/sema_2record_persons.tsv",
  tirz_2record_flag = "output/tirz_2record_persons.tsv",

  # Principal components (first 10 PCs from HRN array data)
  # [HELIX INTERNAL]
  pcs = "output/hrn_pcs.tsv",

  # Binary indicator: individual has genetic data
  # [HELIX INTERNAL]
  has_genetic = "output/has_genetic_data.tsv",

  # Output directory
  output_dir = "output/"
)

###############################################################################
# ANALYSIS PARAMETERS
# Do NOT modify — these match the published analysis
###############################################################################

PARAMS <- list(

  # --- Variant classification thresholds ---
  # AlphaMissense pathogenicity thresholds
  am_damaging_threshold  = 0.564,   # >= this: damaging missense
  am_uncertain_lower     = 0.340,   # >= this (and < damaging): uncertain
  # REVEL thresholds
  revel_damaging_threshold = 0.75,  # >= this: damaging missense
  revel_uncertain_lower    = 0.25,  # >= this (and < damaging): uncertain
  # Variant classified as damaging if AM >= am_damaging_threshold OR REVEL >= revel_damaging_threshold
  # Variant classified as uncertain if:
  #   (AM >= am_uncertain_lower OR REVEL >= revel_uncertain_lower) AND not damaging

  # Functionally validated variant (rs146868158)
  # chr6:39085942:C>T (GRCh38); forms its own category between pLoF and damaging missense
  rs146868158_id = "chr6:39085942:C:T",

  # --- Cohort inclusion criteria ---
  min_bmi           = 25,           # BMI >= 25 kg/m2 required at baseline
  min_followup_days = 182.5,        # >= 6 months follow-up required (182.5 days)
  max_followup_days = 365,          # Weight outcome window: up to 12 months (365 days)

  # --- Matching parameters ---
  match_ratio = 10,                 # 1 carrier : 10 controls (nearest-neighbor)
  match_method = "nearest",         # MatchIt method (nearest-neighbor; supports calipers)
  match_distance = "mahalanobis",   # Distance metric

  # Variables matched exactly
  match_exact = c("drug", "dose_category", "prior_sema_tirz", "prior_other_glp1",
                  "ancestry_group", "sex", "diabetes_type2"),

  # Variables in Mahalanobis distance
  match_mahal = c("age", "bmi_baseline", "followup_duration"),

  # Calipers
  caliper_bmi_kgm2  = 3,           # BMI: +/- 3 kg/m2
  caliper_fu_days   = 30,           # Follow-up duration: +/- 30 days

  # Weights: carriers get weight 1; each control gets 1/n_controls_in_matched_set
  carrier_weight   = 1,

  # --- Regression ---
  vcov_type = "HC1",                # Heteroskedasticity-consistent SE type

  # --- regenie ---
  regenie_version = "2.2.4",
  regenie_ancestry = "EUR",         # Europeans only for regenie sensitivity

  # --- Quality control exclusions ---
  # A small number of individuals were excluded after manual QC review (see Methods).
  # One individual with a variant annotation discrepancy was reclassified
  # from uncertain_missense to benign_missense based on corrected REVEL/AM scores.

  # --- Significance threshold ---
  alpha = 0.05
)

###############################################################################
# ATC DRUG CODES (OMOP CDM / WHO ATC hierarchy)
###############################################################################

ATC <- list(
  semaglutide  = "A10BJ06",
  tirzepatide  = "A10BX16",
  # Other GLP-1 RAs (used for prior GLP-1 RA flag and crossover censoring)
  other_glp1ra = c("A10BJ01",   # exenatide
                   "A10BJ02",   # liraglutide
                   "A10BJ03",   # lixisenatide
                   "A10BJ04",   # albiglutide
                   "A10BJ05",   # dulaglutide
                   "A10BJ07")   # beinaglutide
)

###############################################################################
# COLOUR PALETTE (used in all figures)
###############################################################################

COLORS <- list(
  pink    = "#f45b83",
  yellow  = "#f8be08",
  blue    = "#3f4c77",
  green   = "#7acb71",

  # Shaded variants (lighter = higher index)
  pink_shades   = c("#f45b83", "#f57693", "#f791a4", "#f8acb4",
                    "#fac7c5", "#fbe2d5", "#fdfde6"),
  yellow_shades = c("#f8be08", "#f9c833", "#fad35e", "#fbdd89",
                    "#fce7b4", "#fdf1df", "#fffef0"),
  blue_shades   = c("#3f4c77", "#566283", "#6c788f", "#838e9b",
                    "#9aa4a7", "#b1bab3", "#c8d0bf"),
  green_shades  = c("#7acb71", "#88d280", "#97d98f", "#a6e09e",
                    "#b5e7ad", "#c4eebc", "#d3f5cb"),

  # Drug-specific colours used in figures
  semaglutide = "#f8be08",   # yellow
  tirzepatide = "#7acb71",   # green

  # Variant category colours (used in forest plots)
  plof               = "#3f4c77",   # blue
  damaging_missense  = "#f45b83",   # pink
  rs146868158        = "#f8be08",   # yellow
  uncertain_missense = "#984EA3",   # purple
  benign_missense    = "#c8d0bf"    # light grey
)
