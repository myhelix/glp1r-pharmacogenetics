###############################################################################
# 05_primary_analysis.R
#
# Primary weighted regression analyses on matched datasets.
#
# Steps:
#   1. Load matched datasets (from 04_matching.R)
#   2. Weighted regression: outcome ~ carrier + covariates (HC1 SEs)
#   3. Drug-stratified analyses (semaglutide / tirzepatide)
#   4. Drug × carrier interaction test
#   5. Negative control (GLP1R benign missense)
#   6. GIPR primary analyses
#   7. Save results
#
# Statistical approach:
#   - Weighted linear regression; ATT weights (carrier = 1; each control =
#     1 / n_controls_in_matched_set) as assigned in 04_matching.R
#   - HC1 heteroskedasticity-consistent SEs (sandwich package)
#   - Regression adjusts for age, bmi_baseline, followup_duration
#     (Mahalanobis matching variables; included for precision; exact-match
#     variables are already balanced by design)
#   - Crossover: overall analysis uses each person's first treatment only,
#     as assembled in 03_outcome_definition.R; within the overall matched
#     dataset drug is an exact-match variable, so filtering to a single
#     drug yields self-contained matched sets with valid weights
#
# Outputs (saved to output/):
#   results_primary.rds      — carrier effect estimates across all analyses
#   results_interaction.rds  — drug × carrier interaction test results
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(sandwich)
library(lmtest)

source("config.R")


###############################################################################
# 1. LOAD MATCHED DATASETS
###############################################################################

match_glp1r_overall <- readRDS(file.path(PATHS$output_dir, "matched_glp1r_overall.rds"))
match_glp1r_tirz    <- readRDS(file.path(PATHS$output_dir, "matched_glp1r_tirzepatide.rds"))
match_glp1r_benign  <- readRDS(file.path(PATHS$output_dir, "matched_glp1r_benign.rds"))
match_gipr_overall  <- readRDS(file.path(PATHS$output_dir, "matched_gipr_overall.rds"))
match_gipr_tirz     <- readRDS(file.path(PATHS$output_dir, "matched_gipr_tirzepatide.rds"))
match_gipr_benign   <- readRDS(file.path(PATHS$output_dir, "matched_gipr_benign.rds"))

cat(sprintf(
  "Loaded: GLP1R overall n=%d, GLP1R tirz n=%d, GLP1R benign n=%d, GIPR overall n=%d, GIPR tirz n=%d, GIPR benign n=%d\n",
  nrow(match_glp1r_overall), nrow(match_glp1r_tirz), nrow(match_glp1r_benign),
  nrow(match_gipr_overall),  nrow(match_gipr_tirz),  nrow(match_gipr_benign)
))


###############################################################################
# 2. REGRESSION FUNCTIONS
###############################################################################

adj_vars <- c("age", "bmi_baseline", "followup_duration")

#' Weighted linear regression: carrier effect on % weight change
#'
#' @param df          Matched data frame (with weights_att column)
#' @param carrier_var Name of binary carrier column
#' @param label       Label for output row
#' @return Tibble with n, weighted means, estimate, HC1 SE, 95% CI, p-value
run_regression <- function(df, carrier_var, label) {
  avars <- adj_vars[adj_vars %in% names(df)]
  adj_str <- if (length(avars) > 0) paste("+", paste(avars, collapse = " + ")) else ""
  formula_obj <- as.formula(
    paste("outcome_weight_pct_change ~", carrier_var, adj_str)
  )

  fit <- lm(formula_obj, data = df, weights = df$weights_att)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = PARAMS$vcov_type))

  if (!carrier_var %in% rownames(ct)) {
    cat(sprintf("  WARNING: %s not found in model for %s\n", carrier_var, label))
    return(NULL)
  }
  coef_row <- ct[carrier_var, ]
  est  <- coef_row["Estimate"]
  se   <- coef_row["Std. Error"]
  pval <- coef_row["Pr(>|t|)"]

  carriers <- df[df[[carrier_var]] == 1, ]
  controls <- df[df[[carrier_var]] == 0, ]
  mean_c <- weighted.mean(carriers$outcome_weight_pct_change, carriers$weights_att, na.rm = TRUE)
  mean_k <- weighted.mean(controls$outcome_weight_pct_change, controls$weights_att, na.rm = TRUE)

  tibble(
    analysis     = label,
    n_carriers   = nrow(carriers),
    n_controls   = nrow(controls),
    mean_carrier = round(mean_c, 2),
    mean_control = round(mean_k, 2),
    estimate     = round(est,    2),
    se           = round(se,     3),
    ci_lo        = round(est - 1.96 * se, 2),
    ci_hi        = round(est + 1.96 * se, 2),
    p_value      = pval
  )
}

#' Drug × carrier interaction test within a matched dataset
#'
#' Tests whether the carrier effect differs between semaglutide and tirzepatide.
#' Semaglutide is the reference level (alphabetically first).
#' The coefficient of interest is carrier:drugtirzepatide, representing the
#' additional carrier effect in tirzepatide-treated individuals relative to
#' semaglutide-treated individuals.
#'
#' @param df          Matched data frame containing both sema and tirz
#' @param carrier_var Name of binary carrier column
#' @param label       Label
run_interaction <- function(df, carrier_var, label) {
  if (length(unique(df$drug)) < 2) {
    cat(sprintf("  Skipped (%s): fewer than 2 drug levels in dataset\n", label))
    return(NULL)
  }

  avars <- adj_vars[adj_vars %in% names(df)]
  adj_str <- if (length(avars) > 0) paste("+", paste(avars, collapse = " + ")) else ""
  formula_obj <- as.formula(
    paste("outcome_weight_pct_change ~", carrier_var, "* drug", adj_str)
  )

  fit <- lm(formula_obj, data = df, weights = df$weights_att)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = PARAMS$vcov_type))

  # Identify interaction term (drug levels sorted alphabetically; sema = ref)
  interaction_term <- rownames(ct)[grepl(paste0(carrier_var, ":"), rownames(ct))]
  if (length(interaction_term) == 0) {
    cat(sprintf("  WARNING: no interaction term found for %s\n", label))
    return(NULL)
  }
  interaction_term <- interaction_term[1]

  coef_row <- ct[interaction_term, ]
  est  <- coef_row["Estimate"]
  se   <- coef_row["Std. Error"]
  pval <- coef_row["Pr(>|t|)"]

  tibble(
    analysis         = label,
    interaction_term = interaction_term,
    estimate         = round(est,  2),
    se               = round(se,   3),
    ci_lo            = round(est - 1.96 * se, 2),
    ci_hi            = round(est + 1.96 * se, 2),
    p_value          = pval
  )
}


###############################################################################
# 3. PRIMARY REGRESSION ANALYSES
###############################################################################

cat("\n=== PRIMARY REGRESSION ANALYSES ===\n")

# --- GLP1R primary: overall (first treatment per person) ---
res_glp1r_overall <- run_regression(
  match_glp1r_overall, "carrier_primary",
  "GLP1R primary — overall"
)

# --- GLP1R primary: semaglutide-treated (within overall matched dataset)
# drug is an exact-match variable, so all individuals in a matched set share
# the same drug; filtering to semaglutide yields self-contained matched sets.
res_glp1r_sema <- run_regression(
  match_glp1r_overall %>% filter(drug == "semaglutide"),
  "carrier_primary",
  "GLP1R primary — semaglutide"
)

# --- GLP1R primary: tirzepatide-treated (within overall matched dataset) ---
res_glp1r_tirz_within <- run_regression(
  match_glp1r_overall %>% filter(drug == "tirzepatide"),
  "carrier_primary",
  "GLP1R primary — tirzepatide (within overall)"
)

# --- GLP1R primary: tirzepatide-specific matched dataset
# (carriers and controls drawn exclusively from tirzepatide-treated cohort;
#  provides tirzepatide-specific estimate with maximised control pool)
res_glp1r_tirz_specific <- run_regression(
  match_glp1r_tirz, "carrier_primary",
  "GLP1R primary — tirzepatide-specific"
)

# --- GLP1R benign missense: negative control (overall, sema, tirz) ---
res_glp1r_benign <- run_regression(
  match_glp1r_benign, "carrier_benign",
  "GLP1R benign missense — overall (negative control)"
)
res_glp1r_benign_sema <- run_regression(
  match_glp1r_benign %>% filter(drug == "semaglutide"), "carrier_benign",
  "GLP1R benign missense — semaglutide (negative control)"
)
res_glp1r_benign_tirz <- run_regression(
  match_glp1r_benign %>% filter(drug == "tirzepatide"), "carrier_benign",
  "GLP1R benign missense — tirzepatide (negative control)"
)

# --- GIPR primary: overall ---
res_gipr_overall <- run_regression(
  match_gipr_overall, "carrier_gipr_primary",
  "GIPR primary — overall"
)

# --- GIPR primary: semaglutide-treated (within overall matched dataset) ---
res_gipr_sema <- run_regression(
  match_gipr_overall %>% filter(drug == "semaglutide"),
  "carrier_gipr_primary",
  "GIPR primary — semaglutide"
)

# --- GIPR primary: tirzepatide-specific ---
res_gipr_tirz_specific <- run_regression(
  match_gipr_tirz, "carrier_gipr_primary",
  "GIPR primary — tirzepatide-specific"
)

# --- GIPR benign missense: negative control (overall, sema, tirz) ---
res_gipr_benign <- run_regression(
  match_gipr_benign, "carrier_gipr_benign",
  "GIPR benign missense — overall (negative control)"
)
res_gipr_benign_sema <- run_regression(
  match_gipr_benign %>% filter(drug == "semaglutide"), "carrier_gipr_benign",
  "GIPR benign missense — semaglutide (negative control)"
)
res_gipr_benign_tirz <- run_regression(
  match_gipr_benign %>% filter(drug == "tirzepatide"), "carrier_gipr_benign",
  "GIPR benign missense — tirzepatide (negative control)"
)

results_primary <- bind_rows(
  res_glp1r_overall,
  res_glp1r_sema,
  res_glp1r_tirz_within,
  res_glp1r_tirz_specific,
  res_glp1r_benign,
  res_glp1r_benign_sema,
  res_glp1r_benign_tirz,
  res_gipr_overall,
  res_gipr_sema,
  res_gipr_tirz_specific,
  res_gipr_benign,
  res_gipr_benign_sema,
  res_gipr_benign_tirz
)

cat("\nPrimary results (estimate = carrier vs. control weight % change, pp):\n")
print(
  results_primary %>%
    mutate(
      p_fmt = formatC(p_value, digits = 3, format = "g"),
      ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
    ) %>%
    select(analysis, n_carriers, n_controls,
           mean_carrier, mean_control, estimate, ci, p_fmt),
  n = Inf, width = 120
)


###############################################################################
# 4. DRUG × CARRIER INTERACTION TEST
###############################################################################

cat("\n=== DRUG x CARRIER INTERACTION ===\n")
cat("Tests whether the carrier effect differs between semaglutide and tirzepatide.\n")
cat("semaglutide = reference; interaction term = additional effect in tirzepatide.\n\n")

# GLP1R: uses overall matched dataset (contains sema and tirz, exact-matched on drug)
interaction_glp1r <- run_interaction(
  match_glp1r_overall, "carrier_primary",
  "GLP1R carrier × drug interaction"
)

# GIPR: same approach
interaction_gipr <- run_interaction(
  match_gipr_overall, "carrier_gipr_primary",
  "GIPR carrier × drug interaction"
)

results_interaction <- bind_rows(interaction_glp1r, interaction_gipr)

cat("\nInteraction results:\n")
print(
  results_interaction %>%
    mutate(
      p_fmt = formatC(p_value, digits = 3, format = "g"),
      ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
    ) %>%
    select(analysis, interaction_term, estimate, ci, p_fmt),
  n = Inf, width = 120
)


###############################################################################
# 5. SAVE
###############################################################################

saveRDS(results_primary,     file.path(PATHS$output_dir, "results_primary.rds"))
saveRDS(results_interaction, file.path(PATHS$output_dir, "results_interaction.rds"))

cat("\nResults saved to", PATHS$output_dir, "\n")
