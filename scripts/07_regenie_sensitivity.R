###############################################################################
# 07_regenie_sensitivity.R
#
# Prepare regenie input files and parse regenie output for the European-only
# sensitivity analysis.
#
# Regenie tests the GLP1R and GIPR rare variant burden against weight %
# change in European-ancestry individuals, adjusting for population structure
# (PCs 1–10), age, BMI, sex, drug, and prior GLP-1 RA use.  This complements
# the matched analysis by explicitly accounting for relatedness and ancestry
# via whole-genome ridge regression (regenie step 1).
#
# Steps:
#   1. Load semaglutide + tirzepatide analysis datasets; restrict to Europeans
#   2. Rank-based inverse normal transform (RINT) of the outcome
#   3. Write phenotype file and covariate file for regenie
#   4. Write variant annotation, set-list, and mask files for burden testing
#   5. Parse regenie step-2 output (if present)
#   6. Save formatted results
#
# Shell commands:
#   See scripts/07_regenie_sensitivity.sh for the regenie step-1 and step-2
#   invocations.  Run that script first, then re-source this file to parse
#   the output.
#
# Outputs (written to output/regenie/):
#   regenie_pheno.txt       — phenotype file (IID, RINT weight change)
#   regenie_covariates.txt  — covariate file (age, BMI, sex, drug, PCs 1–10)
#   regenie_anno.txt        — variant annotation file (for burden masks)
#   regenie_setlist.txt     — set-list file
#   regenie_masks.txt       — mask definitions
#   results_regenie.rds     — parsed regenie results (if output present)
#
# Requires: regenie v2.2.4 (see scripts/07_regenie_sensitivity.sh)
#           run on Helix Research Network infrastructure (see Methods)
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(RNOmni)

source("config.R")

regenie_dir <- file.path(PATHS$output_dir, "regenie")
dir.create(regenie_dir, showWarnings = FALSE)


###############################################################################
# 1. LOAD AND FILTER TO EUROPEANS
###############################################################################

sema_analysis <- readRDS(file.path(PATHS$output_dir, "analysis_semaglutide.rds"))
tirz_analysis <- readRDS(file.path(PATHS$output_dir, "analysis_tirzepatide.rds"))

# Overall dataset: first treatment per person (same logic as 03_outcome_definition.R)
overall_analysis <- bind_rows(sema_analysis, tirz_analysis) %>%
  group_by(person_source_value) %>%
  arrange(drug_start_date) %>%
  slice(1) %>%
  ungroup()

# Restrict to European-ancestry individuals
# (ancestry_group assigned from HRN ancestry inference; see Methods)
eur_analysis <- overall_analysis %>%
  filter(ancestry_group == "European")

cat(sprintf(
  "European-ancestry subset: N=%d (GLP1R primary carriers: %d; GIPR primary carriers: %d)\n",
  nrow(eur_analysis),
  sum(eur_analysis$carrier_primary,      na.rm = TRUE),
  sum(eur_analysis$carrier_gipr_primary, na.rm = TRUE)
))


###############################################################################
# 2. RANK-BASED INVERSE NORMAL TRANSFORMATION (RINT)
#
# The primary outcome (minimum % weight change at 6–12 months) is right-skewed.
# regenie sensitivity analysis uses RINT to better satisfy normality assumptions.
# Primary matched analyses use the untransformed outcome (see 05_primary_analysis.R).
###############################################################################

eur_analysis <- eur_analysis %>%
  mutate(outcome_rint = RankNorm(outcome_weight_pct_change))

cat(sprintf(
  "Outcome: raw mean=%.2f%%, SD=%.2f%%; RINT mean=%.4f, SD=%.4f\n",
  mean(eur_analysis$outcome_weight_pct_change, na.rm = TRUE),
  sd(eur_analysis$outcome_weight_pct_change,   na.rm = TRUE),
  mean(eur_analysis$outcome_rint, na.rm = TRUE),
  sd(eur_analysis$outcome_rint,   na.rm = TRUE)
))


###############################################################################
# 3. WRITE PHENOTYPE AND COVARIATE FILES
###############################################################################

# --- Phenotype file ---
# Format: FID IID phenotype  (space-separated; NA coded as NA)
pheno_df <- eur_analysis %>%
  select(FID = person_source_value,
         IID = person_source_value,
         outcome_weight_pct_change_rint = outcome_rint)

write_delim(pheno_df,
            file.path(regenie_dir, "regenie_pheno.txt"),
            delim = " ", na = "NA")
cat(sprintf("Wrote phenotype file: %d individuals\n", nrow(pheno_df)))

# --- Covariate file ---
# Age, BMI, sex (binary: 1=male, 2=female per PLINK convention),
# drug (binary: 1=semaglutide, 0=tirzepatide), prior_other_glp1,
# and the first 10 principal components.
#
# dose_category and diabetes_type2 are included to mirror the primary
# matched analysis adjustment; regenie encodes categorical variables
# as dummy indicators (or they can be included as binary flags if already coded).

pc_cols <- paste0("PC", 1:10)
pc_cols_present <- pc_cols[pc_cols %in% names(eur_analysis)]

covar_df <- eur_analysis %>%
  mutate(
    sex_code  = case_when(sex == "M" ~ 1L, sex == "F" ~ 2L, TRUE ~ NA_integer_),
    drug_sema = as.integer(drug == "semaglutide"),
    t2d       = as.integer(diabetes_type2 == 1)
  ) %>%
  select(FID = person_source_value,
         IID = person_source_value,
         age, bmi_baseline, sex_code, drug_sema,
         prior_other_glp1, t2d,
         all_of(pc_cols_present))

write_delim(covar_df,
            file.path(regenie_dir, "regenie_covariates.txt"),
            delim = " ", na = "NA")
cat(sprintf("Wrote covariate file: %d individuals, %d covariates\n",
            nrow(covar_df), ncol(covar_df) - 2))


###############################################################################
# 4. VARIANT ANNOTATION, SET-LIST, AND MASK FILES
#
# regenie burden / SKAT tests require three auxiliary files:
#   --anno-file   : variant ID → annotation category
#   --set-list    : gene → list of variant IDs
#   --mask-def    : mask name → comma-separated list of annotation categories
#
# These are built from the GLP1R and GIPR carrier TSVs (output of
# 01_variant_classification.py), cross-referenced with the variant
# classification in 03_outcome_definition.R.
###############################################################################

glp1r_carriers <- read_tsv(PATHS$glp1r_carriers, show_col_types = FALSE)
gipr_carriers  <- read_tsv(PATHS$gipr_carriers,  show_col_types = FALSE)

# Assign annotation categories matching the primary analysis hierarchy
classify_anno <- function(df) {
  df %>%
    mutate(
      REVEL = replace_na(REVEL, 0),
      AM    = replace_na(MANE_CT_am_pathogenicity, 0),
      is_lof = MANE_CT_lof == TRUE |
        MANE_CT_consequence %in% c("stop_gained", "frameshift_variant",
                                   "splice_donor_variant", "splice_acceptor_variant"),
      is_damaging = !is_lof &
        (REVEL >= PARAMS$revel_damaging_threshold | AM >= PARAMS$am_damaging_threshold),
      is_uncertain = !is_lof & !is_damaging &
        ((AM >= PARAMS$am_uncertain_lower & AM < PARAMS$am_damaging_threshold) |
           (REVEL >= PARAMS$revel_uncertain_lower & REVEL < PARAMS$revel_damaging_threshold)),
      anno_category = case_when(
        is_lof      ~ "pLoF",
        is_damaging ~ "damaging_missense",
        is_uncertain ~ "uncertain_missense",
        TRUE        ~ "benign_missense"
      )
    ) %>%
    select(variant_id, MANE_CT_gene, anno_category) %>%
    distinct()
}

glp1r_anno <- classify_anno(glp1r_carriers)
gipr_anno  <- classify_anno(gipr_carriers)

# Add rs146868158 as its own annotation category
rs146_anno <- tibble(
  variant_id    = PARAMS$rs146868158_id,
  MANE_CT_gene  = "GLP1R",
  anno_category = "rs146868158"
)
glp1r_anno <- bind_rows(glp1r_anno, rs146_anno) %>% distinct()

all_anno <- bind_rows(glp1r_anno, gipr_anno)

# anno-file: variant_id  gene  annotation
write_delim(all_anno %>% select(variant_id, MANE_CT_gene, anno_category),
            file.path(regenie_dir, "regenie_anno.txt"),
            delim = " ", col_names = FALSE)

# set-list: gene  chr  pos  variant_id1,variant_id2,...
# chr and pos are extracted from variant_id (chr:pos:ref:alt)
set_list <- all_anno %>%
  mutate(
    chrom = sub(":.*", "", variant_id),
    pos   = as.integer(sub("^[^:]+:([^:]+):.*", "\\1", variant_id))
  ) %>%
  group_by(gene = MANE_CT_gene) %>%
  summarise(
    chrom    = first(chrom),
    pos      = min(pos),
    variants = paste(variant_id, collapse = ","),
    .groups  = "drop"
  )
write_delim(set_list,
            file.path(regenie_dir, "regenie_setlist.txt"),
            delim = " ", col_names = FALSE)

# mask-def: mask_name  annotation1,annotation2,...
# Three masks mirroring the primary analysis groupings:
#   M1: pLoF only
#   M2: pLoF + damaging missense
#   M3: pLoF + rs146868158 + damaging missense + uncertain missense (= carrier_primary)
masks <- tribble(
  ~mask,  ~categories,
  "M1_pLoF",            "pLoF",
  "M2_pLoF_damaging",   "pLoF,damaging_missense",
  "M3_primary",         "pLoF,rs146868158,damaging_missense,uncertain_missense"
)
write_delim(masks,
            file.path(regenie_dir, "regenie_masks.txt"),
            delim = " ", col_names = FALSE)

cat(sprintf(
  "Wrote annotation files: %d GLP1R variants, %d GIPR variants, 3 masks\n",
  nrow(glp1r_anno), nrow(gipr_anno)
))


###############################################################################
# 5. PARSE REGENIE OUTPUT
#
# Expected after running 07_regenie_sensitivity.sh.
# regenie step-2 output files are named:
#   regenie_step2_glp1r_outcome_weight_pct_change_rint.regenie
#   regenie_step2_gipr_outcome_weight_pct_change_rint.regenie
###############################################################################

parse_regenie_output <- function(filepath, label) {
  if (!file.exists(filepath)) {
    cat(sprintf("  regenie output not found: %s\n  Run 07_regenie_sensitivity.sh first.\n",
                filepath))
    return(NULL)
  }
  read_table(filepath, show_col_types = FALSE) %>%
    mutate(
      analysis = label,
      p_value  = 10^(-LOG10P)
    ) %>%
    select(analysis, CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ,
           N, BETA, SE, CHISQ, LOG10P, p_value, TEST, EXTRA)
}

regenie_glp1r <- parse_regenie_output(
  file.path(regenie_dir, "regenie_step2_glp1r_outcome_weight_pct_change_rint.regenie"),
  "GLP1R — regenie burden"
)
regenie_gipr  <- parse_regenie_output(
  file.path(regenie_dir, "regenie_step2_gipr_outcome_weight_pct_change_rint.regenie"),
  "GIPR — regenie burden"
)

results_regenie <- bind_rows(regenie_glp1r, regenie_gipr)

if (!is.null(results_regenie) && nrow(results_regenie) > 0) {
  cat("\nRegenie results:\n")
  print(
    results_regenie %>%
      mutate(p_fmt = formatC(p_value, digits = 3, format = "g")) %>%
      select(analysis, ID, TEST, N, BETA, SE, LOG10P, p_fmt),
    n = Inf
  )
  saveRDS(results_regenie, file.path(PATHS$output_dir, "results_regenie.rds"))
  cat("\nRegenie results saved to", PATHS$output_dir, "\n")
} else {
  cat("\nNo regenie output to parse. Saved NULL placeholder.\n")
  saveRDS(NULL, file.path(PATHS$output_dir, "results_regenie.rds"))
}

cat("\nRegenie input files written to", regenie_dir, "\n")
cat("Run scripts/07_regenie_sensitivity.sh to execute regenie.\n")
