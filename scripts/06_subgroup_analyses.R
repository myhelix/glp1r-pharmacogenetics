###############################################################################
# 06_subgroup_analyses.R
#
# Subgroup and sensitivity analyses.
#
# Steps:
#   1. Load matched datasets (from 04_matching.R)
#   2. GLP1R variant class subgroups
#        — pLoF + damaging missense (most functional)
#        — rs146868158 specifically
#        — uncertain missense
#   3. Clinical subgroup interaction tests (sex, T2D)
#   4. Sensitivity: treatment-naive (no prior sema or tirz)
#   5. GIPR variant class subgroups
#   6. Save
#
# Notes on subsetting matched data:
#   - drug, dose_category, prior_sema_tirz, prior_other_glp1, sex,
#     diabetes_type2, and ancestry_group are exact-match variables.
#     Filtering on these preserves matched-set integrity (all individuals
#     in a set share the same value).
#   - Variant class subgroups are extracted by identifying matched sets
#     (subclass) whose carrier belongs to the target class, then retaining
#     all controls in those sets.
#   - BMI and age are Mahalanobis variables; subgrouping on these would
#     break matched sets, so they are tested as regression interaction terms.
#
# Outputs (saved to output/):
#   results_subgroup.rds      — all subgroup / sensitivity estimates
#   results_subgroup_interaction.rds — clinical subgroup interaction tests
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
match_gipr_overall  <- readRDS(file.path(PATHS$output_dir, "matched_gipr_overall.rds"))
match_gipr_tirz     <- readRDS(file.path(PATHS$output_dir, "matched_gipr_tirzepatide.rds"))


###############################################################################
# HELPER FUNCTIONS
###############################################################################

adj_vars <- c("age", "bmi_baseline", "followup_duration")

#' Weighted regression on a (possibly subsetted) matched data frame.
#' Reuses the existing weights_att without recalculation.
run_regression <- function(df, carrier_var, label) {
  avars <- adj_vars[adj_vars %in% names(df)]
  adj_str <- if (length(avars) > 0) paste("+", paste(avars, collapse = " + ")) else ""
  formula_obj <- as.formula(
    paste("outcome_weight_pct_change ~", carrier_var, adj_str)
  )

  if (sum(df[[carrier_var]] == 1, na.rm = TRUE) < 1) {
    cat(sprintf("  Skipped (%s): no carriers\n", label)); return(NULL)
  }

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

  tibble(
    analysis     = label,
    n_carriers   = nrow(carriers),
    n_controls   = nrow(controls),
    mean_carrier = round(weighted.mean(carriers$outcome_weight_pct_change,
                                       carriers$weights_att, na.rm = TRUE), 2),
    mean_control = round(weighted.mean(controls$outcome_weight_pct_change,
                                       controls$weights_att, na.rm = TRUE), 2),
    estimate     = round(est,  2),
    se           = round(se,   3),
    ci_lo        = round(est - 1.96 * se, 2),
    ci_hi        = round(est + 1.96 * se, 2),
    p_value      = pval
  )
}

#' Subset a matched dataset to sets whose carrier belongs to a variant class.
#' @param df           Matched data frame (contains subclass, carrier_var)
#' @param carrier_var  Primary carrier indicator (used to identify carriers)
#' @param class_var    Variable that flags the target subclass (0/1 or logical)
subset_by_carrier_class <- function(df, carrier_var, class_var) {
  target_subclasses <- df %>%
    filter(.data[[carrier_var]] == 1, .data[[class_var]] == 1) %>%
    pull(subclass) %>%
    unique()
  df %>% filter(subclass %in% target_subclasses)
}

#' Clinical subgroup interaction test: outcome ~ carrier * subgroup_var
run_interaction_subgroup <- function(df, carrier_var, subgroup_var, label) {
  avars <- adj_vars[adj_vars %in% names(df)]
  adj_str <- if (length(avars) > 0) paste("+", paste(avars, collapse = " + ")) else ""
  formula_obj <- as.formula(
    paste("outcome_weight_pct_change ~",
          carrier_var, "*", subgroup_var, adj_str)
  )

  if (length(unique(na.omit(df[[subgroup_var]]))) < 2) {
    cat(sprintf("  Skipped (%s): %s has < 2 levels\n", label, subgroup_var))
    return(NULL)
  }
  if (sum(df[[carrier_var]] == 1, na.rm = TRUE) < 1) {
    cat(sprintf("  Skipped (%s): no carriers\n", label))
    return(NULL)
  }

  fit <- tryCatch(
    lm(formula_obj, data = df, weights = df$weights_att),
    error = function(e) { cat(sprintf("  ERROR (%s): %s\n", label, e$message)); NULL }
  )
  if (is.null(fit)) return(NULL)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = PARAMS$vcov_type))

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
    estimate         = round(est, 2),
    se               = round(se,  3),
    ci_lo            = round(est - 1.96 * se, 2),
    ci_hi            = round(est + 1.96 * se, 2),
    p_value          = pval
  )
}


###############################################################################
# 2. GLP1R VARIANT CLASS SUBGROUPS
#
# Within the overall matched GLP1R dataset, subsets are defined by identifying
# matched sets (subclass) whose carrier belongs to a specific variant class.
# All 10 controls per such carrier are retained; weights are unchanged.
###############################################################################

cat("\n=== GLP1R VARIANT CLASS SUBGROUPS ===\n")

# carrier_lof_damaging: pLoF + damaging missense (highest-confidence functional)
glp1r_lod_overall <- subset_by_carrier_class(
  match_glp1r_overall, "carrier_primary", "carrier_lof_damaging"
)
res_lod_overall <- run_regression(
  glp1r_lod_overall, "carrier_primary",
  "GLP1R pLoF+damaging — overall"
)
res_lod_sema <- run_regression(
  glp1r_lod_overall %>% filter(drug == "semaglutide"),
  "carrier_primary",
  "GLP1R pLoF+damaging — semaglutide"
)
glp1r_lod_tirz <- subset_by_carrier_class(
  match_glp1r_tirz, "carrier_primary", "carrier_lof_damaging"
)
res_lod_tirz <- run_regression(
  glp1r_lod_tirz, "carrier_primary",
  "GLP1R pLoF+damaging — tirzepatide-specific"
)

# pLoF only
glp1r_lof_overall <- subset_by_carrier_class(
  match_glp1r_overall, "carrier_primary", "has_lof"
)
res_lof_overall <- run_regression(
  glp1r_lof_overall, "carrier_primary",
  "GLP1R pLoF — overall"
)
res_lof_sema <- run_regression(
  glp1r_lof_overall %>% filter(drug == "semaglutide"),
  "carrier_primary",
  "GLP1R pLoF — semaglutide"
)
glp1r_lof_tirz <- subset_by_carrier_class(
  match_glp1r_tirz, "carrier_primary", "has_lof"
)
res_lof_tirz <- run_regression(
  glp1r_lof_tirz, "carrier_primary",
  "GLP1R pLoF — tirzepatide-specific"
)

# Damaging missense only
glp1r_dam_overall <- subset_by_carrier_class(
  match_glp1r_overall, "carrier_primary", "has_damaging_missense"
)
res_dam_overall <- run_regression(
  glp1r_dam_overall, "carrier_primary",
  "GLP1R damaging missense — overall"
)
res_dam_sema <- run_regression(
  glp1r_dam_overall %>% filter(drug == "semaglutide"),
  "carrier_primary",
  "GLP1R damaging missense — semaglutide"
)
glp1r_dam_tirz <- subset_by_carrier_class(
  match_glp1r_tirz, "carrier_primary", "has_damaging_missense"
)
res_dam_tirz <- run_regression(
  glp1r_dam_tirz, "carrier_primary",
  "GLP1R damaging missense — tirzepatide-specific"
)

# rs146868158: matched sets where the carrier has variant_class == "rs146868158"
# (identified by highest-severity classification assigned in 03_outcome_definition.R)
if ("variant_class" %in% names(match_glp1r_overall)) {
  glp1r_rs146_sets <- match_glp1r_overall %>%
    filter(carrier_primary == 1, variant_class == "rs146868158") %>%
    pull(subclass) %>% unique()
  glp1r_rs146_df <- match_glp1r_overall %>% filter(subclass %in% glp1r_rs146_sets)
  res_rs146 <- run_regression(
    glp1r_rs146_df, "carrier_primary",
    "GLP1R rs146868158 — overall"
  )
  res_rs146_sema <- run_regression(
    glp1r_rs146_df %>% filter(drug == "semaglutide"),
    "carrier_primary",
    "GLP1R rs146868158 — semaglutide"
  )
} else {
  cat("  variant_class column absent from matched data; rs146868158 subgroup skipped\n")
  res_rs146      <- NULL
  res_rs146_sema <- NULL
}
if ("variant_class" %in% names(match_glp1r_tirz)) {
  glp1r_rs146_tirz_sets <- match_glp1r_tirz %>%
    filter(carrier_primary == 1, variant_class == "rs146868158") %>%
    pull(subclass) %>% unique()
  glp1r_rs146_tirz_df <- match_glp1r_tirz %>% filter(subclass %in% glp1r_rs146_tirz_sets)
  res_rs146_tirz <- run_regression(
    glp1r_rs146_tirz_df, "carrier_primary",
    "GLP1R rs146868158 — tirzepatide-specific"
  )
} else {
  res_rs146_tirz <- NULL
}

# carrier_uncertain: uncertain missense only
glp1r_unc_overall <- subset_by_carrier_class(
  match_glp1r_overall, "carrier_primary", "carrier_uncertain"
)
res_unc_overall <- run_regression(
  glp1r_unc_overall, "carrier_primary",
  "GLP1R uncertain missense — overall"
)

res_unc_sema <- run_regression(
  glp1r_unc_overall %>% filter(drug == "semaglutide"),
  "carrier_primary",
  "GLP1R uncertain missense — semaglutide"
)

glp1r_unc_tirz <- subset_by_carrier_class(
  match_glp1r_tirz, "carrier_primary", "carrier_uncertain"
)
res_unc_tirz <- run_regression(
  glp1r_unc_tirz, "carrier_primary",
  "GLP1R uncertain missense — tirzepatide-specific"
)


###############################################################################
# 3. CLINICAL SUBGROUP INTERACTION TESTS
#
# Test whether the GLP1R carrier effect is modified by sex or T2D status.
# These use interaction terms rather than subsetting, because sex and
# diabetes_type2 are exact-match variables (all individuals in a matched set
# share the same value), meaning each stratum would be an isolated matched
# set — the interaction approach uses all data while allowing effect modification.
###############################################################################

cat("\n=== CLINICAL SUBGROUP INTERACTIONS ===\n")

int_sex <- run_interaction_subgroup(
  match_glp1r_overall, "carrier_primary", "sex",
  "GLP1R primary × sex"
)

int_t2d <- run_interaction_subgroup(
  match_glp1r_overall, "carrier_primary", "diabetes_type2",
  "GLP1R primary × diabetes_type2"
)

results_subgroup_interaction <- bind_rows(int_sex, int_t2d)

cat("\nClinical subgroup interaction results:\n")
print(
  results_subgroup_interaction %>%
    mutate(
      p_fmt = formatC(p_value, digits = 3, format = "g"),
      ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
    ) %>%
    select(analysis, interaction_term, estimate, ci, p_fmt),
  n = Inf
)


###############################################################################
# 4. SENSITIVITY: TREATMENT-NAIVE ONLY
#
# prior_sema_tirz is an exact-match variable. Filtering to prior_sema_tirz == 0
# retains only matched sets composed entirely of GLP-1 RA-naive individuals;
# matched-set integrity is preserved and weights are unchanged.
###############################################################################

cat("\n=== SENSITIVITY: TREATMENT-NAIVE ===\n")

res_naive_overall <- run_regression(
  match_glp1r_overall %>% filter(prior_sema_tirz == 0),
  "carrier_primary",
  "GLP1R primary — treatment-naive, overall"
)

res_naive_sema <- run_regression(
  match_glp1r_overall %>% filter(prior_sema_tirz == 0, drug == "semaglutide"),
  "carrier_primary",
  "GLP1R primary — treatment-naive, semaglutide"
)

res_naive_tirz <- run_regression(
  match_glp1r_tirz %>% filter(prior_sema_tirz == 0),
  "carrier_primary",
  "GLP1R primary — treatment-naive, tirzepatide-specific"
)


###############################################################################
# 5. GIPR VARIANT CLASS SUBGROUPS
###############################################################################

cat("\n=== GIPR VARIANT CLASS SUBGROUPS ===\n")

# carrier_gipr_lof_damaging flag is pre-computed in 03_outcome_definition.R
if ("carrier_gipr_lof_damaging" %in% names(match_gipr_overall)) {

  gipr_lod_overall <- subset_by_carrier_class(
    match_gipr_overall, "carrier_gipr_primary", "carrier_gipr_lof_damaging"
  )
  res_gipr_lod_overall <- run_regression(
    gipr_lod_overall, "carrier_gipr_primary",
    "GIPR pLoF+damaging — overall"
  )
  res_gipr_lod_sema <- run_regression(
    gipr_lod_overall %>% filter(drug == "semaglutide"),
    "carrier_gipr_primary",
    "GIPR pLoF+damaging — semaglutide"
  )
  gipr_lod_tirz <- subset_by_carrier_class(
    match_gipr_tirz, "carrier_gipr_primary", "carrier_gipr_lof_damaging"
  )
  res_gipr_lod_tirz <- run_regression(
    gipr_lod_tirz, "carrier_gipr_primary",
    "GIPR pLoF+damaging — tirzepatide-specific"
  )

} else {
  cat("  carrier_gipr_lof_damaging absent; skipping GIPR pLoF+damaging subgroups\n")
  res_gipr_lod_overall <- NULL
  res_gipr_lod_sema    <- NULL
  res_gipr_lod_tirz    <- NULL
}

# GIPR pLoF only
if ("has_gipr_lof" %in% names(match_gipr_overall)) {
  gipr_lof_overall <- subset_by_carrier_class(
    match_gipr_overall, "carrier_gipr_primary", "has_gipr_lof"
  )
  res_gipr_lof_overall <- run_regression(
    gipr_lof_overall, "carrier_gipr_primary", "GIPR pLoF — overall"
  )
  res_gipr_lof_sema <- run_regression(
    gipr_lof_overall %>% filter(drug == "semaglutide"),
    "carrier_gipr_primary", "GIPR pLoF — semaglutide"
  )
  gipr_lof_tirz <- subset_by_carrier_class(
    match_gipr_tirz, "carrier_gipr_primary", "has_gipr_lof"
  )
  res_gipr_lof_tirz <- run_regression(
    gipr_lof_tirz, "carrier_gipr_primary", "GIPR pLoF — tirzepatide-specific"
  )
} else {
  res_gipr_lof_overall <- NULL; res_gipr_lof_sema <- NULL; res_gipr_lof_tirz <- NULL
}

# GIPR damaging missense only
if ("has_gipr_damaging_missense" %in% names(match_gipr_overall)) {
  gipr_dam_overall <- subset_by_carrier_class(
    match_gipr_overall, "carrier_gipr_primary", "has_gipr_damaging_missense"
  )
  res_gipr_dam_overall <- run_regression(
    gipr_dam_overall, "carrier_gipr_primary", "GIPR damaging missense — overall"
  )
  res_gipr_dam_sema <- run_regression(
    gipr_dam_overall %>% filter(drug == "semaglutide"),
    "carrier_gipr_primary", "GIPR damaging missense — semaglutide"
  )
  gipr_dam_tirz <- subset_by_carrier_class(
    match_gipr_tirz, "carrier_gipr_primary", "has_gipr_damaging_missense"
  )
  res_gipr_dam_tirz <- run_regression(
    gipr_dam_tirz, "carrier_gipr_primary",
    "GIPR damaging missense — tirzepatide-specific"
  )
} else {
  res_gipr_dam_overall <- NULL; res_gipr_dam_sema <- NULL; res_gipr_dam_tirz <- NULL
}

# GIPR uncertain missense
if ("has_gipr_uncertain_missense" %in% names(match_gipr_overall)) {

  gipr_unc_overall <- subset_by_carrier_class(
    match_gipr_overall, "carrier_gipr_primary", "has_gipr_uncertain_missense"
  )
  res_gipr_unc_overall <- run_regression(
    gipr_unc_overall, "carrier_gipr_primary",
    "GIPR uncertain missense — overall"
  )
  res_gipr_unc_sema <- run_regression(
    gipr_unc_overall %>% filter(drug == "semaglutide"),
    "carrier_gipr_primary",
    "GIPR uncertain missense — semaglutide"
  )

  gipr_unc_tirz <- subset_by_carrier_class(
    match_gipr_tirz, "carrier_gipr_primary", "has_gipr_uncertain_missense"
  )
  res_gipr_unc_tirz <- run_regression(
    gipr_unc_tirz, "carrier_gipr_primary",
    "GIPR uncertain missense — tirzepatide-specific"
  )

} else {
  cat("  has_gipr_uncertain_missense absent; skipping GIPR uncertain missense subgroups\n")
  res_gipr_unc_overall <- NULL
  res_gipr_unc_sema    <- NULL
  res_gipr_unc_tirz    <- NULL
}


###############################################################################
# 6. CLINICAL SUBGROUP ANALYSES (GLP1R AND GIPR PRIMARY CARRIERS)
#
# Estimates the carrier effect within clinical subgroups defined by:
#   - T2D status (exact-match variable; matched-set integrity preserved)
#   - Sex (exact-match variable; matched-set integrity preserved)
#   - Prior GLP-1 RA use (exact-match variable; matched-set integrity preserved)
#   - BMI <40 vs ≥40 (Mahalanobis variable; approximate subgroup analysis)
#   - High-dose therapy (dose_category; exact-match variable)
#
# Each subgroup is run for: overall, semaglutide, tirzepatide-specific.
###############################################################################

cat("\n=== CLINICAL SUBGROUP ANALYSES ===\n")

#' Run clinical subgroup analyses for one gene across all drug strata
#' @param df_overall   Overall matched dataset
#' @param df_tirz      Tirzepatide-specific matched dataset
#' @param carrier_var  Carrier indicator column name
#' @param gene         Gene label ("GLP1R" or "GIPR")
run_clinical_subgroups <- function(df_overall, df_tirz, carrier_var, gene) {

  results <- list()

  # Helper: run regression with a filter applied to df_overall (both drugs)
  # and the same filter applied to df_tirz for tirzepatide-specific estimates.
  run_sg <- function(filter_expr, suffix) {
    df_o <- tryCatch(df_overall %>% filter({{ filter_expr }}), error = function(e) NULL)
    df_s <- tryCatch(df_overall %>% filter(drug == "semaglutide") %>%
                       filter({{ filter_expr }}), error = function(e) NULL)
    df_t <- tryCatch(df_tirz   %>% filter({{ filter_expr }}), error = function(e) NULL)
    list(
      overall = if (!is.null(df_o)) run_regression(df_o, carrier_var, paste0(gene, " primary — ", suffix, ", overall"))    else NULL,
      sema    = if (!is.null(df_s)) run_regression(df_s, carrier_var, paste0(gene, " primary — ", suffix, ", semaglutide")) else NULL,
      tirz    = if (!is.null(df_t)) run_regression(df_t, carrier_var, paste0(gene, " primary — ", suffix, ", tirzepatide-specific")) else NULL
    )
  }

  # T2D
  r <- run_sg(diabetes_type2 == 1, "T2D");    results <- c(results, r)
  r <- run_sg(diabetes_type2 == 0, "No T2D"); results <- c(results, r)

  # Sex
  r <- run_sg(sex == "M", "Male");   results <- c(results, r)
  r <- run_sg(sex == "F", "Female"); results <- c(results, r)

  # GLP-1 RA experience
  if (all(c("prior_sema_tirz", "prior_other_glp1") %in% names(df_overall))) {
    r <- run_sg(prior_sema_tirz == 0 & prior_other_glp1 == 0, "GLP-1 naive");
    results <- c(results, r)
    r <- run_sg(prior_sema_tirz == 1 | prior_other_glp1 == 1, "GLP-1 experienced");
    results <- c(results, r)
  }

  # BMI subgroups (approximate; bmi_baseline is Mahalanobis, not exact-match)
  if ("bmi_baseline" %in% names(df_overall)) {
    r <- run_sg(bmi_baseline <  40, "BMI <40");  results <- c(results, r)
    r <- run_sg(bmi_baseline >= 40, "BMI \u226540"); results <- c(results, r)
  }

  # High dose
  if ("dose_category" %in% names(df_overall)) {
    r <- run_sg(grepl("high", dose_category, ignore.case = TRUE), "high dose");
    results <- c(results, r)
  }

  bind_rows(results)
}

results_clinical_glp1r <- run_clinical_subgroups(
  match_glp1r_overall, match_glp1r_tirz, "carrier_primary", "GLP1R"
)
cat("\nGLP1R clinical subgroup results:\n")
print(results_clinical_glp1r %>%
        select(analysis, n_carriers, n_controls, estimate, p_value) %>%
        mutate(p_value = formatC(p_value, digits = 3, format = "g")),
      n = Inf, width = 120)

results_clinical_gipr <- run_clinical_subgroups(
  match_gipr_overall, match_gipr_tirz, "carrier_gipr_primary", "GIPR"
)
cat("\nGIPR clinical subgroup results:\n")
print(results_clinical_gipr %>%
        select(analysis, n_carriers, n_controls, estimate, p_value) %>%
        mutate(p_value = formatC(p_value, digits = 3, format = "g")),
      n = Inf, width = 120)


###############################################################################
# 7. COMPILE AND PRINT
###############################################################################

results_subgroup <- bind_rows(
  res_lod_overall,
  res_lod_sema,
  res_lod_tirz,
  res_lof_overall,
  res_lof_sema,
  res_lof_tirz,
  res_dam_overall,
  res_dam_sema,
  res_dam_tirz,
  res_rs146,
  res_rs146_sema,
  res_rs146_tirz,
  res_unc_overall,
  res_unc_sema,
  res_unc_tirz,
  res_naive_overall,
  res_naive_sema,
  res_naive_tirz,
  res_gipr_lod_overall,
  res_gipr_lod_sema,
  res_gipr_lod_tirz,
  res_gipr_lof_overall,
  res_gipr_lof_sema,
  res_gipr_lof_tirz,
  res_gipr_dam_overall,
  res_gipr_dam_sema,
  res_gipr_dam_tirz,
  res_gipr_unc_overall,
  res_gipr_unc_sema,
  res_gipr_unc_tirz,
  results_clinical_glp1r,
  results_clinical_gipr
)

cat("\nSubgroup results summary:\n")
print(
  results_subgroup %>%
    mutate(
      p_fmt = formatC(p_value, digits = 3, format = "g"),
      ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
    ) %>%
    select(analysis, n_carriers, n_controls, estimate, ci, p_fmt),
  n = Inf, width = 120
)


###############################################################################
# 8. SAVE
###############################################################################

saveRDS(results_subgroup,             file.path(PATHS$output_dir, "results_subgroup.rds"))
saveRDS(results_subgroup_interaction, file.path(PATHS$output_dir, "results_subgroup_interaction.rds"))

cat("\nSubgroup results saved to", PATHS$output_dir, "\n")
