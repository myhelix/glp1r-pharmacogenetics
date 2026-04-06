###############################################################################
# 04_matching.R
#
# Match GLP1R/GIPR rare variant carriers to controls using 1:10 nearest-
# neighbor matching with Mahalanobis distance, exact matching, and
# BMI/follow-up calipers.
#
# Steps:
#   1. Load analysis-ready datasets (from 03_outcome_definition.R)
#   2. Outlier inspection in matched sample (see note below)
#   3. Check for crossover carriers (individuals in both drug cohorts)
#   4. Perform matching for each variant category
#   5. Assess covariate balance (standardized mean differences)
#   6. Save matched datasets
#
# Outputs (saved to output/):
#   matched_glp1r_overall.rds      — GLP1R primary carriers, overall
#   matched_glp1r_tirzepatide.rds  — GLP1R primary carriers, tirzepatide only
#   matched_glp1r_benign.rds       — GLP1R benign missense (negative control)
#   matched_gipr_overall.rds       — GIPR primary carriers, overall
#   matched_gipr_tirzepatide.rds   — GIPR primary carriers, tirzepatide only
#   matched_gipr_benign.rds        — GIPR benign missense (negative control)
#   balance_summary.rds            — SMD table for all matched datasets
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

library(tidyverse)
library(MatchIt)
library(sandwich)
library(lmtest)

source("config.R")

###############################################################################
# 1. LOAD DATA
###############################################################################

overall_analysis  <- readRDS(file.path(PATHS$output_dir, "analysis_overall.rds"))
sema_analysis     <- readRDS(file.path(PATHS$output_dir, "analysis_semaglutide.rds"))
tirz_analysis     <- readRDS(file.path(PATHS$output_dir, "analysis_tirzepatide.rds"))

cat(sprintf("Loaded: overall n=%d, sema n=%d, tirz n=%d\n",
            nrow(overall_analysis), nrow(sema_analysis), nrow(tirz_analysis)))


###############################################################################
# 2. OUTLIER INSPECTION
# Per the Methods, weight changes beyond the reported threshold are inspected
# manually within the matched sample before finalising results.
# Flag extreme values here for review; do not remove automatically.
###############################################################################

flag_outliers <- function(df, outcome_var = "outcome_weight_pct_change",
                          threshold = 60) {
  df %>%
    mutate(outlier_flag = abs(.data[[outcome_var]]) > threshold)
}

overall_analysis <- flag_outliers(overall_analysis)
sema_analysis    <- flag_outliers(sema_analysis)
tirz_analysis    <- flag_outliers(tirz_analysis)

cat(sprintf(
  "Outlier flags (|weight change| > 60%%): overall=%d, sema=%d, tirz=%d\n",
  sum(overall_analysis$outlier_flag, na.rm = TRUE),
  sum(sema_analysis$outlier_flag,    na.rm = TRUE),
  sum(tirz_analysis$outlier_flag,    na.rm = TRUE)
))
# Review flagged individuals before proceeding; remove if confirmed errors


###############################################################################
# 3. CROSSOVER CARRIER CHECK
# Some carriers appear in both sema and tirz cohorts. The overall analysis
# uses each person's first treatment only (handled in 03_outcome_definition.R).
# For the tirzepatide-specific analysis, we use the tirz_analysis dataset
# directly (independent of crossover status).
###############################################################################

check_crossover <- function(sema_df, tirz_df, carrier_var) {
  sema_ids <- sema_df %>% filter(.data[[carrier_var]] == 1) %>%
    pull(person_source_value)
  tirz_ids <- tirz_df %>% filter(.data[[carrier_var]] == 1) %>%
    pull(person_source_value)
  crossover <- intersect(sema_ids, tirz_ids)
  cat(sprintf("  %s: sema=%d, tirz=%d, crossover=%d\n",
              carrier_var, length(sema_ids), length(tirz_ids), length(crossover)))
  crossover
}

cat("\nCrossover carrier check:\n")
check_crossover(sema_analysis, tirz_analysis, "carrier_primary")
check_crossover(sema_analysis, tirz_analysis, "carrier_benign")
check_crossover(sema_analysis, tirz_analysis, "carrier_gipr_primary")
check_crossover(sema_analysis, tirz_analysis, "carrier_gipr_benign")


###############################################################################
# 4. MATCHING FUNCTION
###############################################################################

#' Perform 1:10 nearest-neighbor matching
#'
#' @param data         Data frame containing carriers and potential controls
#' @param carrier_var  Name of binary carrier indicator column (0/1)
#' @param label        Label for printed output
#' @param ratio        Carrier:control ratio (default 10)
#'
#' Matching specification (see Methods and config.R):
#'   - Method: nearest-neighbor without replacement (MatchIt method="nearest")
#'   - Distance: Mahalanobis on age, bmi_baseline, followup_duration
#'   - Exact: drug, dose_category, prior_sema_tirz, prior_other_glp1,
#'            ancestry_group, sex, diabetes_type2
#'   - Calipers: BMI ±3 kg/m², follow-up ±30 days (raw, not standardised)
#'   - Weights: carrier = 1; each control = 1 / n_controls_in_matched_set

perform_matching <- function(data, carrier_var, label,
                             ratio = PARAMS$match_ratio) {
  cat(sprintf("\n--- %s ---\n", label))

  n_carriers <- sum(data[[carrier_var]] == 1, na.rm = TRUE)
  cat(sprintf("  Input: N=%d | carriers=%d\n", nrow(data), n_carriers))

  if (n_carriers < 1) {
    cat("  Skipped: no carriers\n")
    return(NULL)
  }

  # Drop exact-match variables that are constant or absent in this subset
  exact_vars <- PARAMS$match_exact
  exact_vars <- exact_vars[exact_vars %in% names(data)]
  exact_vars <- exact_vars[sapply(exact_vars, function(v) {
    n_levels <- length(unique(na.omit(data[[v]])))
    carrier_n_levels <- length(unique(na.omit(
      data[[v]][data[[carrier_var]] == 1]
    )))
    n_levels >= 2 && carrier_n_levels >= 1
  })]

  mahal_vars   <- PARAMS$match_mahal[PARAMS$match_mahal %in% names(data)]
  all_vars     <- unique(c(exact_vars, mahal_vars))

  match_formula <- as.formula(
    paste(carrier_var, "~", paste(all_vars, collapse = " + "))
  )
  exact_formula <- if (length(exact_vars) > 0)
    as.formula(paste("~", paste(exact_vars, collapse = " + ")))
  else NULL

  calipers <- c(
    bmi_baseline      = PARAMS$caliper_bmi_kgm2,
    followup_duration = PARAMS$caliper_fu_days
  )
  calipers <- calipers[names(calipers) %in% names(data)]

  cat(sprintf("  Exact: %s\n",  paste(exact_vars,  collapse = ", ")))
  cat(sprintf("  Mahal: %s\n",  paste(mahal_vars,  collapse = ", ")))
  cat(sprintf("  Calipers: %s\n",
              paste(names(calipers), calipers, sep = "=±", collapse = ", ")))

  matched_obj <- tryCatch(
    matchit(
      match_formula,
      data       = data,
      method     = PARAMS$match_method,
      distance   = PARAMS$match_distance,
      exact      = exact_formula,
      caliper    = calipers,
      std.caliper = FALSE,
      ratio      = ratio,
      replace    = FALSE,
      estimand   = "ATT"
    ),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )

  if (is.null(matched_obj)) return(NULL)

  matched_df <- match.data(matched_obj) %>%
    group_by(subclass) %>%
    mutate(
      n_controls_in_set = sum(.data[[carrier_var]] == 0),
      weights_att = if_else(.data[[carrier_var]] == 1,
                            1,
                            1 / n_controls_in_set)
    ) %>%
    ungroup()

  n_c <- sum(matched_df[[carrier_var]] == 1)
  n_k <- sum(matched_df[[carrier_var]] == 0)
  cat(sprintf("  Matched: %d carriers, %d controls (ratio %.1f:1)\n",
              n_c, n_k, n_k / n_c))

  list(matched_data    = matched_df,
       matchit_object  = matched_obj,
       n_carriers      = n_c,
       n_controls      = n_k,
       carrier_var     = carrier_var)
}


###############################################################################
# 5. RUN MATCHING
###############################################################################

# Controls pool: individuals with no GLP1R/GIPR variant of any kind
# Exclude GIPR benign carriers so they are not used as controls in any analysis
no_variant <- overall_analysis %>%
  filter(any_glp1r_variant == 0 & carrier_gipr_primary == 0 &
           carrier_gipr_benign == 0)

# --- GLP1R PRIMARY (overall: first treatment per person) ---
match_glp1r_overall <- perform_matching(
  data        = bind_rows(
    overall_analysis %>% filter(carrier_primary == 1),
    no_variant
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_primary",
  label       = "GLP1R primary — overall"
)

# --- GLP1R PRIMARY (tirzepatide-specific) ---
match_glp1r_tirz <- perform_matching(
  data        = bind_rows(
    tirz_analysis %>% filter(carrier_primary == 1),
    tirz_analysis %>% filter(any_glp1r_variant == 0 & carrier_gipr_primary == 0)
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_primary",
  label       = "GLP1R primary — tirzepatide"
)

# --- GLP1R BENIGN MISSENSE (negative control; overall) ---
match_glp1r_benign <- perform_matching(
  data        = bind_rows(
    overall_analysis %>% filter(carrier_benign == 1),
    no_variant
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_benign",
  label       = "GLP1R benign missense — overall"
)

# --- GIPR PRIMARY (overall) ---
match_gipr_overall <- perform_matching(
  data        = bind_rows(
    overall_analysis %>% filter(carrier_gipr_primary == 1),
    no_variant
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_gipr_primary",
  label       = "GIPR primary — overall"
)

# --- GIPR PRIMARY (tirzepatide-specific) ---
match_gipr_tirz <- perform_matching(
  data        = bind_rows(
    tirz_analysis %>% filter(carrier_gipr_primary == 1),
    tirz_analysis %>% filter(any_glp1r_variant == 0 & carrier_gipr_primary == 0)
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_gipr_primary",
  label       = "GIPR primary — tirzepatide"
)

# --- GIPR BENIGN MISSENSE (negative control; overall) ---
match_gipr_benign <- perform_matching(
  data        = bind_rows(
    overall_analysis %>% filter(carrier_gipr_benign == 1),
    no_variant
  ) %>% mutate(ancestry_group = droplevels(ancestry_group)),
  carrier_var = "carrier_gipr_benign",
  label       = "GIPR benign missense — overall"
)


###############################################################################
# 6. COVARIATE BALANCE (STANDARDIZED MEAN DIFFERENCES)
###############################################################################

#' Calculate weighted SMD for one variable
calc_smd <- function(data, var, carrier_var, weight_var = "weights_att") {
  carriers <- data[data[[carrier_var]] == 1, ]
  controls <- data[data[[carrier_var]] == 0, ]

  if (is.numeric(data[[var]])) {
    c_mean <- weighted.mean(carriers[[var]], carriers[[weight_var]], na.rm = TRUE)
    k_mean <- weighted.mean(controls[[var]], controls[[weight_var]], na.rm = TRUE)
    c_var  <- weighted.mean((carriers[[var]] - c_mean)^2,
                            carriers[[weight_var]], na.rm = TRUE)
    k_var  <- weighted.mean((controls[[var]] - k_mean)^2,
                            controls[[weight_var]], na.rm = TRUE)
    pooled_sd <- sqrt((c_var + k_var) / 2)
    if (pooled_sd == 0) return(0)
    return((c_mean - k_mean) / pooled_sd)
  } else {
    lvls <- sort(unique(na.omit(data[[var]])))
    if (length(lvls) < 2) return(0)
    ref <- lvls[1]
    c_prop <- weighted.mean(carriers[[var]] == ref,
                            carriers[[weight_var]], na.rm = TRUE)
    k_prop <- weighted.mean(controls[[var]] == ref,
                            controls[[weight_var]], na.rm = TRUE)
    pooled_sd <- sqrt((c_prop * (1 - c_prop) + k_prop * (1 - k_prop)) / 2)
    if (pooled_sd == 0) return(0)
    return((c_prop - k_prop) / pooled_sd)
  }
}

balance_vars <- c("age", "bmi_baseline", "followup_duration",
                  "sex", "diabetes_type2", "prior_sema_tirz",
                  "prior_other_glp1", "ancestry_group")

compute_balance <- function(match_result, label) {
  if (is.null(match_result)) return(NULL)
  df  <- match_result$matched_data
  cvar <- match_result$carrier_var
  smds <- sapply(balance_vars[balance_vars %in% names(df)],
                 calc_smd, data = df, carrier_var = cvar)
  tibble(analysis   = label,
         variable   = names(smds),
         smd        = smds,
         balanced   = abs(smds) < 0.10)
}

balance_summary <- bind_rows(
  compute_balance(match_glp1r_overall, "GLP1R overall"),
  compute_balance(match_glp1r_tirz,    "GLP1R tirzepatide"),
  compute_balance(match_glp1r_benign,  "GLP1R benign"),
  compute_balance(match_gipr_overall,  "GIPR overall"),
  compute_balance(match_gipr_tirz,     "GIPR tirzepatide"),
  compute_balance(match_gipr_benign,   "GIPR benign")
)

cat("\nCovariate balance (SMD; <0.10 = well-balanced):\n")
print(balance_summary %>%
        mutate(smd = round(smd, 3)) %>%
        pivot_wider(names_from = analysis, values_from = c(smd, balanced)),
      n = Inf)


###############################################################################
# 7. SAVE
###############################################################################

saveRDS(match_glp1r_overall$matched_data,
        file.path(PATHS$output_dir, "matched_glp1r_overall.rds"))
saveRDS(match_glp1r_tirz$matched_data,
        file.path(PATHS$output_dir, "matched_glp1r_tirzepatide.rds"))
saveRDS(match_glp1r_benign$matched_data,
        file.path(PATHS$output_dir, "matched_glp1r_benign.rds"))
saveRDS(match_gipr_overall$matched_data,
        file.path(PATHS$output_dir, "matched_gipr_overall.rds"))
saveRDS(match_gipr_tirz$matched_data,
        file.path(PATHS$output_dir, "matched_gipr_tirzepatide.rds"))
saveRDS(match_gipr_benign$matched_data,
        file.path(PATHS$output_dir, "matched_gipr_benign.rds"))
saveRDS(balance_summary,
        file.path(PATHS$output_dir, "balance_summary.rds"))

cat("\nMatched datasets saved to", PATHS$output_dir, "\n")
