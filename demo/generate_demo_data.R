###############################################################################
# demo/generate_demo_data.R
#
# Generate fully synthetic demo data that mimics the structure and column
# layout expected by scripts 03–09.  No real patient data is included.
#
# All values are simulated.  The synthetic GLP1R carrier effect (attenuated
# weight loss with semaglutide, preserved with tirzepatide) is imposed by
# design so that downstream scripts produce interpretable output.
#
# Run this script once before running demo/run_demo.R.
#
# Outputs (written to output/ as defined in config.R):
#   glp1r_rare_coding_carriers.tsv
#   gipr_rare_coding_carriers.tsv
#   rs146868158_carriers.tsv
#   semaglutide_cohort.RData          (object: sema_cohort)
#   tirzepatide_cohort.RData          (object: tirz_cohort)
#   semaglutide_longitudinal_weights.RData (object: sema_longitudinal)
#   tirzepatide_longitudinal_weights.RData (object: tirz_longitudinal)
#   sema_2record_persons.tsv
#   tirz_2record_persons.tsv
#   hrn_pcs.tsv
#   has_genetic_data.tsv
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)

source("config.R")
set.seed(20250101)

dir.create(PATHS$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating synthetic demo data...\n")


###############################################################################
# PARAMETERS
###############################################################################

N_SEMA    <- 10000   # semaglutide cohort size: must be large because overall_analysis only
                     # uses first-treatment-per-person, halving the effective sema control
                     # pool (vs tirz-specific which uses all tirz_analysis directly)
N_TIRZ    <- 6000    # tirzepatide cohort size (large enough for tirz-specific matching)
N_OVERLAP <- 400     # individuals appearing in both cohorts (crossover)

# Carrier counts per drug arm — increased substantially so that after strict
# 1:10 nearest-neighbour matching, every clinical subgroup (T2D/no-T2D,
# male/female, BMI <40/≥40, GLP-1 naive/experienced, high dose) has enough
# carriers to produce visible forest-plot rows in all three drug strata.
N_GLP1R_LOF_SEMA       <- 20;  N_GLP1R_LOF_TIRZ       <- 18
N_GLP1R_RS146_SEMA     <- 20;  N_GLP1R_RS146_TIRZ      <- 18
N_GLP1R_DAMAGING_SEMA  <- 22;  N_GLP1R_DAMAGING_TIRZ   <- 20
N_GLP1R_UNCERTAIN_SEMA <- 22;  N_GLP1R_UNCERTAIN_TIRZ  <- 20
N_GLP1R_BENIGN_SEMA    <- 12;  N_GLP1R_BENIGN_TIRZ     <- 10
N_GIPR_LOF_SEMA        <- 15;  N_GIPR_LOF_TIRZ         <- 15
N_GIPR_DAMAGING_SEMA   <- 22;  N_GIPR_DAMAGING_TIRZ    <- 20
N_GIPR_UNCERTAIN_SEMA  <- 22;  N_GIPR_UNCERTAIN_TIRZ   <- 20
N_GIPR_BENIGN_SEMA     <- 12;  N_GIPR_BENIGN_TIRZ      <- 10

# Weight change parameters (mean %, SD) by group
# Negative = weight loss
WC <- list(
  sema_noncarrier  = list(mu = -8.0,  sd = 5.5),
  tirz_noncarrier  = list(mu = -12.0, sd = 6.0),
  # GLP1R: attenuated response with sema, normal with tirz
  glp1r_sema       = list(mu = -3.5,  sd = 5.0),  # attenuated
  glp1r_tirz       = list(mu = -11.8, sd = 6.0),  # same as non-carrier
  # GIPR: attenuated response with tirz
  gipr_sema        = list(mu = -7.8,  sd = 5.5),  # same as non-carrier
  gipr_tirz        = list(mu = -6.0,  sd = 5.5)   # attenuated
)


###############################################################################
# 1. SYNTHETIC VARIANT TABLES
###############################################################################

# GLP1R gene region: chr6:39,016,000-39,086,000 (GRCh38)
# GIPR gene region:  chr19:46,154,000-46,174,000 (GRCh38)

make_variant_id <- function(chrom, pos, ref, alt) {
  sprintf("%s:%d:%s:%s", chrom, pos, ref, alt)
}

# --- Synthetic GLP1R variants ---
glp1r_variants <- tribble(
  ~variant_id,                          ~MANE_CT_consequence,    ~MANE_CT_HGVSP_custom, ~REVEL, ~MANE_CT_am_pathogenicity, ~MANE_CT_lof,
  "chr6:39021445:C:T",  "stop_gained",          "p.Arg176Ter",   NA,   NA,  TRUE,
  "chr6:39031882:GA:G", "frameshift_variant",   "p.Leu247fs",    NA,   NA,     TRUE,
  "chr6:39044213:C:A",  "splice_donor_variant", "p.?",           NA,   NA,     TRUE,
  "chr6:39019876:A:G",  "stop_gained",          "p.Trp112Ter",   NA,   NA,     TRUE,
  "chr6:39085942:C:T",  "missense_variant",     "p.Arg68Cys",    0.62, 0.48,   FALSE,  # rs146868158
  "chr6:39085942:C:T",  "missense_variant",     "p.Arg68Cys",    0.62, 0.48,   FALSE,  # rs146868158 dup row for 2nd carrier
  "chr6:39022100:G:A",  "missense_variant",     "p.Gly186Arg",   0.81, 0.71,   FALSE,
  "chr6:39027543:T:C",  "missense_variant",     "p.Ile241Thr",   0.78, 0.68,   FALSE,
  "chr6:39033291:A:T",  "missense_variant",     "p.Asn298Ile",   0.76, 0.65,   FALSE,
  "chr6:39040118:C:G",  "missense_variant",     "p.Pro347Ala",   0.77, 0.72,   FALSE,
  "chr6:39018774:G:C",  "missense_variant",     "p.Asp114His",   0.79, 0.69,   FALSE,
  "chr6:39046892:T:A",  "missense_variant",     "p.Leu385His",   0.75, 0.58,   FALSE,
  "chr6:39050341:C:T",  "missense_variant",     "p.Pro412Leu",   0.58, 0.42,   FALSE,
  "chr6:39023654:A:G",  "missense_variant",     "p.Ile201Val",   0.52, 0.37,   FALSE,
  "chr6:39029877:G:A",  "missense_variant",     "p.Ala261Thr",   0.49, 0.35,   FALSE,
  "chr6:39035100:T:C",  "missense_variant",     "p.Tyr321His",   0.44, 0.34,   FALSE,
  "chr6:39037826:G:T",  "missense_variant",     "p.Glu341Asp",   0.30, 0.22,   FALSE,
  "chr6:39042019:C:A",  "missense_variant",     "p.Arg367Ser",   0.28, 0.18,   FALSE,
  "chr6:39048733:A:G",  "missense_variant",     "p.Met398Val",   0.22, 0.15,   FALSE,
  "chr6:39053410:T:C",  "missense_variant",     "p.Lys422Arg",   0.19, 0.12,   FALSE
) %>%
  mutate(MANE_CT_gene = "GLP1R", n_alt_alleles = 1L)

# --- Synthetic GIPR variants ---
gipr_variants <- tribble(
  ~variant_id,                          ~MANE_CT_consequence,    ~MANE_CT_HGVSP_custom, ~REVEL, ~MANE_CT_am_pathogenicity, ~MANE_CT_lof,
  "chr19:46157823:C:T", "stop_gained",       "p.Arg143Ter",   NA,   NA,     TRUE,
  "chr19:46159441:AT:A","frameshift_variant","p.Thr201fs",    NA,   NA,     TRUE,
  "chr19:46162114:G:A", "missense_variant",  "p.Val278Met",   0.80, 0.67,   FALSE,
  "chr19:46163890:T:C", "missense_variant",  "p.Ile319Thr",   0.77, 0.63,   FALSE,
  "chr19:46165201:A:G", "missense_variant",  "p.Lys356Glu",   0.76, 0.61,   FALSE,
  "chr19:46166734:C:T", "missense_variant",  "p.Arg398Cys",   0.75, 0.59,   FALSE,
  "chr19:46168042:G:C", "missense_variant",  "p.Asp432His",   0.78, 0.65,   FALSE,
  "chr19:46160381:T:A", "missense_variant",  "p.Tyr226Asn",   0.55, 0.40,   FALSE,
  "chr19:46161792:G:T", "missense_variant",  "p.Ala271Ser",   0.51, 0.36,   FALSE,
  "chr19:46163100:A:C", "missense_variant",  "p.Ile311Leu",   0.47, 0.35,   FALSE,
  "chr19:46164500:T:G", "missense_variant",  "p.Phe352Val",   0.44, 0.34,   FALSE,
  "chr19:46165900:C:A", "missense_variant",  "p.Pro396Gln",   0.42, 0.33,   FALSE,
  "chr19:46167200:G:A", "missense_variant",  "p.Ala432Thr",   0.40, 0.32,   FALSE,
  # Benign missense (REVEL < 0.25 AND AM < 0.34)
  "chr19:46158900:G:T", "missense_variant",  "p.Ala192Val",   0.22, 0.20,   FALSE,
  "chr19:46160100:C:A", "missense_variant",  "p.Pro246Gln",   0.19, 0.17,   FALSE,
  "chr19:46161400:T:G", "missense_variant",  "p.Ser268Ala",   0.18, 0.15,   FALSE,
  "chr19:46162800:A:C", "missense_variant",  "p.Leu305Pro",   0.15, 0.12,   FALSE
) %>%
  mutate(MANE_CT_gene = "GIPR", n_alt_alleles = 1L)


###############################################################################
# 2. ASSIGN CARRIERS TO SYNTHETIC INDIVIDUALS
###############################################################################

# All individuals across both cohorts (unique)
all_ids <- sprintf("DEMO_%05d", seq_len(N_SEMA + N_TIRZ - N_OVERLAP))
sema_ids <- all_ids[seq_len(N_SEMA)]
tirz_ids <- c(all_ids[(N_SEMA - N_OVERLAP + 1):N_SEMA],         # overlap
              all_ids[(N_SEMA + 1):(N_SEMA + N_TIRZ - N_OVERLAP)])

# Split pool into sema-only, tirz-only, and overlap so we can guarantee
# representation in both drug arms for every variant class.
sema_only_ids <- setdiff(sema_ids, tirz_ids)
tirz_only_ids <- setdiff(tirz_ids, sema_ids)

# Draw n_s carriers from sema-only pool and n_t from tirz-only pool.
# excluded: character vector of already-assigned IDs to avoid duplicates.
draw_both_arms <- function(n_s, n_t, exclude = character(0)) {
  s <- sample(setdiff(sema_only_ids, exclude), n_s)
  t <- sample(setdiff(tirz_only_ids, c(exclude, s)), n_t)
  c(s, t)
}

used <- character(0)
glp1r_lof_carriers   <- draw_both_arms(N_GLP1R_LOF_SEMA,       N_GLP1R_LOF_TIRZ,       used); used <- c(used, glp1r_lof_carriers)
glp1r_rs146_carriers <- draw_both_arms(N_GLP1R_RS146_SEMA,     N_GLP1R_RS146_TIRZ,     used); used <- c(used, glp1r_rs146_carriers)
glp1r_dam_carriers   <- draw_both_arms(N_GLP1R_DAMAGING_SEMA,  N_GLP1R_DAMAGING_TIRZ,  used); used <- c(used, glp1r_dam_carriers)
glp1r_unc_carriers   <- draw_both_arms(N_GLP1R_UNCERTAIN_SEMA, N_GLP1R_UNCERTAIN_TIRZ, used); used <- c(used, glp1r_unc_carriers)
glp1r_ben_carriers   <- draw_both_arms(N_GLP1R_BENIGN_SEMA,    N_GLP1R_BENIGN_TIRZ,    used); used <- c(used, glp1r_ben_carriers)
gipr_lof_carriers    <- draw_both_arms(N_GIPR_LOF_SEMA,        N_GIPR_LOF_TIRZ,        used); used <- c(used, gipr_lof_carriers)
gipr_dam_carriers    <- draw_both_arms(N_GIPR_DAMAGING_SEMA,   N_GIPR_DAMAGING_TIRZ,   used); used <- c(used, gipr_dam_carriers)
gipr_unc_carriers    <- draw_both_arms(N_GIPR_UNCERTAIN_SEMA,  N_GIPR_UNCERTAIN_TIRZ,  used); used <- c(used, gipr_unc_carriers)
gipr_ben_carriers    <- draw_both_arms(N_GIPR_BENIGN_SEMA,     N_GIPR_BENIGN_TIRZ,     used); used <- c(used, gipr_ben_carriers)

# Build carrier-variant pair tables
assign_carriers_to_variants <- function(carrier_ids, variant_rows) {
  # Cycle through variants when there are more carriers than variant rows
  n_v <- nrow(variant_rows)
  map_dfr(seq_along(carrier_ids), function(i) {
    variant_rows[((i - 1) %% n_v) + 1, ] %>%
      mutate(d_id = carrier_ids[i])
  })
}

glp1r_carrier_rows <- bind_rows(
  assign_carriers_to_variants(glp1r_lof_carriers,   glp1r_variants %>% filter(MANE_CT_lof)),
  assign_carriers_to_variants(glp1r_rs146_carriers, glp1r_variants %>% filter(variant_id == "chr6:39085942:C:T") %>% slice(1)),
  assign_carriers_to_variants(glp1r_dam_carriers,   glp1r_variants %>% filter(!MANE_CT_lof, REVEL >= 0.75) %>% slice(1:6)),
  assign_carriers_to_variants(glp1r_unc_carriers,   glp1r_variants %>% filter(!MANE_CT_lof, !is.na(REVEL), REVEL >= 0.25, REVEL < 0.75)),
  assign_carriers_to_variants(glp1r_ben_carriers,   glp1r_variants %>% filter(!MANE_CT_lof, !is.na(REVEL), REVEL < 0.25))
) %>%
  select(d_id, variant_id, n_alt_alleles, MANE_CT_gene,
         MANE_CT_consequence, MANE_CT_HGVSP_custom,
         REVEL, MANE_CT_am_pathogenicity, MANE_CT_lof)

gipr_carrier_rows <- bind_rows(
  assign_carriers_to_variants(gipr_lof_carriers, gipr_variants %>% filter(MANE_CT_lof)),
  assign_carriers_to_variants(gipr_dam_carriers, gipr_variants %>% filter(!MANE_CT_lof, REVEL >= 0.75)),
  assign_carriers_to_variants(gipr_unc_carriers, gipr_variants %>% filter(!MANE_CT_lof, !is.na(REVEL), REVEL >= 0.25, REVEL < 0.75)),
  assign_carriers_to_variants(gipr_ben_carriers, gipr_variants %>% filter(!MANE_CT_lof, !is.na(REVEL), REVEL < 0.25))
) %>%
  select(d_id, variant_id, n_alt_alleles, MANE_CT_gene,
         MANE_CT_consequence, MANE_CT_HGVSP_custom,
         REVEL, MANE_CT_am_pathogenicity, MANE_CT_lof)

rs146_rows <- tibble(
  d_id          = glp1r_rs146_carriers,
  n_alt_alleles = 1L
)

write_tsv(glp1r_carrier_rows, file.path(PATHS$output_dir, "glp1r_rare_coding_carriers.tsv"))
write_tsv(gipr_carrier_rows,  file.path(PATHS$output_dir, "gipr_rare_coding_carriers.tsv"))
write_tsv(rs146_rows,         file.path(PATHS$output_dir, "rs146868158_carriers.tsv"))
cat("  Written: carrier TSVs\n")


###############################################################################
# 3. COHORT DATA FRAMES
#
# Column layout must match what 03_outcome_definition.R expects after
# 02_cohort_selection.R (include flag, bmi_baseline, followup_duration,
# drug_start_date, dose_category, prior_sema_tirz, prior_other_glp1,
# ancestry_group, sex, diabetes_type2).
###############################################################################

gen_ancestry <- function(n) {
  sample(c("European","African","Americas","EastAsian","SouthAsian","Other"),
         n, replace = TRUE,
         prob = c(0.62, 0.14, 0.12, 0.06, 0.04, 0.02))
}

gen_cohort <- function(ids, drug_label) {
  n <- length(ids)
  tibble(
    person_source_value = ids,
    drug_start_date     = as.Date("2020-01-01") +
                          sample(0:1460, n, replace = TRUE),
    age                 = round(rnorm(n, mean = 52, sd = 11), 1),
    sex                 = sample(c("M","F"), n, replace = TRUE, prob = c(0.38, 0.62)),
    bmi_baseline        = round(pmax(25.5, rnorm(n, mean = 36, sd = 7)), 1),
    followup_duration   = round(runif(n, min = 185, max = 365)),
    dose_category       = sample(c("standard","high"), n, replace = TRUE,
                                 prob = c(0.65, 0.35)),
    prior_sema_tirz     = as.integer(rbinom(n, 1, prob = 0.05)),
    prior_other_glp1    = as.integer(rbinom(n, 1, prob = 0.12)),
    diabetes_type2      = as.integer(rbinom(n, 1, prob = 0.45)),
    ancestry_group  = gen_ancestry(n),
    include             = 1L
  )
}

sema_cohort <- gen_cohort(sema_ids, "semaglutide")
tirz_cohort <- gen_cohort(tirz_ids, "tirzepatide")

cat("  Generated cohorts: sema n=", nrow(sema_cohort),
    " tirz n=", nrow(tirz_cohort), "\n")


###############################################################################
# 3b. BALANCE CARRIER DEMOGRAPHICS
#
# Override the randomly-assigned demographics of carrier individuals to ensure
# each clinical subgroup (T2D/no-T2D, male/female, BMI <40/≥40,
# GLP-1 naive/experienced, high dose) has enough carriers in all three drug
# strata to produce visible forest-plot rows.
#
# Key decisions:
#   - followup_duration set to ≥275 days so all carriers pass the 6–12-month
#     measurement filter in 03_outcome_definition.R.
#   - Proportions chosen to match the overall cohort distribution so that the
#     carrier pool is matchable to the large non-carrier control pool.
#   - Prior-GLP-1 RA split: 0% prior_sema_tirz=1, 10% prior_other_glp1=1,
#     remaining naive (both are exact-match variables).
###############################################################################

# Allocate exactly round(n*p) individuals to each binary subgroup category
set_balanced <- function(n, p) {
  n1 <- max(1L, as.integer(round(n * p)))
  n0 <- n - n1
  sample(c(rep(1L, n1), rep(0L, n0)))
}

balance_carrier_demographics <- function(cohort_df, carrier_ids_in_arm) {
  idx <- which(cohort_df$person_source_value %in% carrier_ids_in_arm)
  n   <- length(idx)
  if (n == 0) return(cohort_df)

  # T2D: 45%
  cohort_df$diabetes_type2[idx]  <- set_balanced(n, 0.45)

  # Sex: 40% male
  cohort_df$sex[idx]              <- ifelse(set_balanced(n, 0.40) == 1L, "M", "F")

  # BMI: 35% ≥40 kg/m².
  # Draw from N(33, 3.5) for low-BMI carriers (centred near the control
  # distribution peak) and N(42.5, 1.5) for high-BMI carriers.  Add targeted
  # jitter to any values that land exactly at a boundary (27.0 or 38.9 for low;
  # 40.5 for high) to prevent identical-BMI spikes: carriers with identical BMI
  # compete for the exact same caliper-window controls, causing later carriers
  # to find zero controls in 1:10 nearest-neighbour matching.
  is_bmi_high <- set_balanced(n, 0.35) == 1L
  n_high <- sum(is_bmi_high)
  n_low  <- n - n_high
  bmi_low  <- round(pmin(38.9, pmax(27.0, rnorm(max(1L, n_low),  33.0, 3.5))), 1)
  at_ceil  <- bmi_low >= 38.9
  if (any(at_ceil)) bmi_low[at_ceil] <- round(runif(sum(at_ceil), 37.5, 38.8), 1)
  at_flr   <- bmi_low <= 27.0
  if (any(at_flr))  bmi_low[at_flr]  <- round(runif(sum(at_flr),  27.1, 28.5), 1)
  bmi_high <- round(pmax(40.5, rnorm(max(1L, n_high), 42.5, 1.5)), 1)
  at_flr_h <- bmi_high <= 40.5
  if (any(at_flr_h)) bmi_high[at_flr_h] <- round(runif(sum(at_flr_h), 40.6, 41.0), 1)
  all_bmi <- numeric(n)
  if (n_low  > 0) all_bmi[!is_bmi_high] <- bmi_low
  if (n_high > 0) all_bmi[is_bmi_high]  <- bmi_high
  cohort_df$bmi_baseline[idx] <- all_bmi

  # Prior GLP-1 RA:
  #   prior_sema_tirz = 0 for ALL carriers — the 5% population rate produces
  #   so few controls with prior sema/tirz that carriers in this strata almost
  #   never find 10 matching controls (the exact-match constraint kills the pool).
  #   prior_other_glp1 = 10% of carriers (needed for the GLP-1-experienced
  #   subgroup row in the ED forest plot).  gen_cohort sets the background rate
  #   at 12% so there are sufficient controls available in each caliper window.
  pog <- as.integer(set_balanced(n, 0.10) == 1L)
  cohort_df$prior_sema_tirz[idx]  <- 0L
  cohort_df$prior_other_glp1[idx] <- pog

  # Dose: 35% high dose
  cohort_df$dose_category[idx] <- ifelse(set_balanced(n, 0.35) == 1L,
                                          "high", "standard")

  # Ancestry: 90% European to maximise matchability
  cohort_df$ancestry_group[idx] <- sample(
    c("European", "African"),
    n, replace = TRUE, prob = c(0.90, 0.10)
  )

  # Follow-up: ≥275 days so all carriers reach the day-270 measurement
  # (the earliest measurement in the 6–12 month outcome window)
  cohort_df$followup_duration[idx] <- round(runif(n, 275, 365))

  cohort_df
}

# Identify sema-arm and tirz-arm primary carriers
all_glp1r_primary_ids <- c(glp1r_lof_carriers, glp1r_rs146_carriers,
                             glp1r_dam_carriers, glp1r_unc_carriers)
all_gipr_primary_ids  <- c(gipr_lof_carriers, gipr_dam_carriers, gipr_unc_carriers)

glp1r_sema_primary <- intersect(all_glp1r_primary_ids, sema_only_ids)
glp1r_tirz_primary <- intersect(all_glp1r_primary_ids, tirz_only_ids)
gipr_sema_primary  <- intersect(all_gipr_primary_ids,  sema_only_ids)
gipr_tirz_primary  <- intersect(all_gipr_primary_ids,  tirz_only_ids)

# Benign carriers also need balanced demographics so they have valid follow-up
# measurements and land in matchable strata (without balance they have random
# followup 185-365 days, and ~47% miss the day-270 measurement window).
glp1r_ben_sema <- intersect(glp1r_ben_carriers, sema_only_ids)
glp1r_ben_tirz <- intersect(glp1r_ben_carriers, tirz_only_ids)
gipr_ben_sema  <- intersect(gipr_ben_carriers,  sema_only_ids)
gipr_ben_tirz  <- intersect(gipr_ben_carriers,  tirz_only_ids)

sema_cohort <- balance_carrier_demographics(sema_cohort, glp1r_sema_primary)
tirz_cohort <- balance_carrier_demographics(tirz_cohort, glp1r_tirz_primary)
sema_cohort <- balance_carrier_demographics(sema_cohort, gipr_sema_primary)
tirz_cohort <- balance_carrier_demographics(tirz_cohort, gipr_tirz_primary)
sema_cohort <- balance_carrier_demographics(sema_cohort, glp1r_ben_sema)
tirz_cohort <- balance_carrier_demographics(tirz_cohort, glp1r_ben_tirz)
sema_cohort <- balance_carrier_demographics(sema_cohort, gipr_ben_sema)
tirz_cohort <- balance_carrier_demographics(tirz_cohort, gipr_ben_tirz)


###############################################################################
# 4. LONGITUDINAL WEIGHT DATA
#
# Simulate measurements at days 90, 180, 270, 365.
# The primary outcome uses the minimum value in days 182.5–365.
# GLP1R carriers: attenuated weight loss with sema; normal with tirz.
# GIPR carriers: attenuated weight loss with tirz; normal with sema.
###############################################################################

# Classify carrier status per individual
glp1r_primary_ids <- all_glp1r_primary_ids
gipr_primary_ids  <- all_gipr_primary_ids

gen_longitudinal <- function(cohort_df, drug_label) {
  timepoints <- c(90, 180, 270, 365)

  map_dfr(seq_len(nrow(cohort_df)), function(i) {
    pid   <- cohort_df$person_source_value[i]
    fu    <- cohort_df$followup_duration[i]

    # Determine weight change trajectory parameters
    if (pid %in% glp1r_primary_ids && drug_label == "semaglutide") {
      # GLP1R carrier on sema — attenuated effect
      final_wc <- rnorm(1, WC$glp1r_sema$mu, WC$glp1r_sema$sd)
    } else if (pid %in% gipr_primary_ids && drug_label == "tirzepatide") {
      # GIPR carrier on tirz — attenuated effect
      final_wc <- rnorm(1, WC$gipr_tirz$mu, WC$gipr_tirz$sd)
    } else if (drug_label == "semaglutide") {
      final_wc <- rnorm(1, WC$sema_noncarrier$mu, WC$sema_noncarrier$sd)
    } else {
      final_wc <- rnorm(1, WC$tirz_noncarrier$mu, WC$tirz_noncarrier$sd)
    }

    # Generate measurements up to followup_duration
    valid_tp <- timepoints[timepoints <= fu]
    if (length(valid_tp) == 0) return(NULL)

    # Weight change increases (more negative) over time with noise
    wc_at_tp <- final_wc * (valid_tp / 365) +
                rnorm(length(valid_tp), 0, 1.5)

    tibble(
      person_source_value = pid,
      days_from_start     = valid_tp,
      weight_pct_change   = round(wc_at_tp, 2)
    )
  })
}

cat("  Generating longitudinal weight data (this may take a moment)...\n")
sema_longitudinal <- gen_longitudinal(sema_cohort, "semaglutide")
tirz_longitudinal <- gen_longitudinal(tirz_cohort, "tirzepatide")

cat(sprintf("  Longitudinal rows: sema=%d, tirz=%d\n",
            nrow(sema_longitudinal), nrow(tirz_longitudinal)))


###############################################################################
# 5. ANCILLARY FILES
###############################################################################

# 2-record flag files (all pass in demo)
sema_2rec <- sema_cohort %>%
  select(person_source_value) %>%
  mutate(meets_any_7_90_within_90d_index = 1L)
tirz_2rec <- tirz_cohort %>%
  select(person_source_value) %>%
  mutate(meets_any_7_90_within_90d_index = 1L)

write_tsv(sema_2rec, file.path(PATHS$output_dir, "sema_2record_persons.tsv"))
write_tsv(tirz_2rec, file.path(PATHS$output_dir, "tirz_2record_persons.tsv"))

# Principal components (random; used only in regenie sensitivity)
all_unique_ids <- unique(c(sema_ids, tirz_ids))
pcs_df <- tibble(FID = all_unique_ids, IID = all_unique_ids) %>%
  bind_cols(
    as_tibble(matrix(rnorm(length(all_unique_ids) * 10, 0, 1),
                     ncol = 10,
                     dimnames = list(NULL, paste0("PC", 1:10))))
  )
write_tsv(pcs_df, file.path(PATHS$output_dir, "hrn_pcs.tsv"))

# has_genetic_data flag (all individuals have genetic data in demo)
has_genetic <- tibble(person_source_value = all_unique_ids)
write_tsv(has_genetic, file.path(PATHS$output_dir, "has_genetic_data.tsv"))

cat("  Written: ancillary files\n")


###############################################################################
# 6. SAVE .RData FILES
###############################################################################

save(sema_cohort,      file = file.path(PATHS$output_dir, "semaglutide_cohort.RData"))
save(tirz_cohort,      file = file.path(PATHS$output_dir, "tirzepatide_cohort.RData"))
save(sema_longitudinal, file = file.path(PATHS$output_dir, "semaglutide_longitudinal_weights.RData"))
save(tirz_longitudinal, file = file.path(PATHS$output_dir, "tirzepatide_longitudinal_weights.RData"))
cat("  Written: .RData cohort and longitudinal files\n")


###############################################################################
# SUMMARY
###############################################################################

cat("\nDemo data summary:\n")
cat(sprintf(
  "  Individuals:    sema=%d, tirz=%d, unique=%d\n",
  nrow(sema_cohort), nrow(tirz_cohort), length(all_unique_ids)
))
cat(sprintf(
  "  GLP1R carriers: pLoF=%d  rs146=%d  damaging=%d  uncertain=%d  benign=%d\n",
  length(glp1r_lof_carriers), length(glp1r_rs146_carriers),
  length(glp1r_dam_carriers), length(glp1r_unc_carriers), length(glp1r_ben_carriers)
))
cat(sprintf(
  "  GIPR carriers:  pLoF=%d  damaging=%d  uncertain=%d  benign=%d\n",
  length(gipr_lof_carriers), length(gipr_dam_carriers),
  length(gipr_unc_carriers),  length(gipr_ben_carriers)
))
cat(sprintf(
  "  Crossover (in both cohorts): %d\n", N_OVERLAP
))
cat(sprintf(
  "  GLP1R primary (sema/tirz arm): %d / %d\n",
  length(glp1r_sema_primary), length(glp1r_tirz_primary)
))
cat(sprintf(
  "  GIPR primary  (sema/tirz arm): %d / %d\n",
  length(gipr_sema_primary), length(gipr_tirz_primary)
))
cat("\nAll demo data written to:", PATHS$output_dir, "\n")
cat("Run demo/run_demo.R to execute the analysis pipeline.\n")
