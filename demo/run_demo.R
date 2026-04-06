###############################################################################
# demo/run_demo.R
#
# End-to-end demonstration of the analysis pipeline using synthetic data.
#
# Steps:
#   0. Generate synthetic demo data (calls demo/generate_demo_data.R)
#   1. Variant classification + outcome definition  (scripts/03_outcome_definition.R)
#   2. Matching                                     (scripts/04_matching.R)
#   3. Primary analysis                             (scripts/05_primary_analysis.R)
#   4. Subgroup analyses                            (scripts/06_subgroup_analyses.R)
#   5. Figures + Table 1                            (scripts/09_figures.R)
#
# Scripts 01 (Hail/HRN), 02 (OMOP CDM/HRN), 07 (regenie/HRN), and
# 08 (All of Us) require infrastructure access and are not run here.
#
# Usage:
#   Rscript demo/run_demo.R
#   — or —
#   source("demo/run_demo.R")   # from within the glp1r-pharmacogenetics/ directory
#
# Requirements:
#   R packages: tidyverse, MatchIt, sandwich, lmtest, RNOmni, patchwork
#   Install with: install.packages(c("tidyverse","MatchIt","sandwich",
#                                    "lmtest","RNOmni","patchwork"))
#
# Paper: "Rare GLP1R variants selectively attenuate weight loss with
#         semaglutide versus tirzepatide"
###############################################################################

# Ensure working directory is glp1r-pharmacogenetics/
if (!file.exists("config.R")) {
  stop("Run this script from the glp1r-pharmacogenetics/ directory:\n",
       "  cd glp1r-pharmacogenetics && Rscript demo/run_demo.R")
}

run_step <- function(label, script_path) {
  cat(sprintf("\n%s\n%s\n%s\n",
              strrep("=", 70), label, strrep("=", 70)))
  tryCatch(
    source(script_path, echo = FALSE),
    error = function(e) {
      cat(sprintf("\nERROR in %s:\n  %s\n", script_path, e$message))
      cat("Check that all required packages are installed.\n")
      stop(e)
    }
  )
  cat(sprintf("\n[DONE] %s\n", label))
}


###############################################################################
# STEP 0: GENERATE DEMO DATA
###############################################################################

run_step("Step 0: Generate synthetic demo data",
         "demo/generate_demo_data.R")


###############################################################################
# STEP 1: VARIANT CLASSIFICATION + OUTCOME DEFINITION
###############################################################################

run_step("Step 1: Variant classification and outcome definition (03)",
         "scripts/03_outcome_definition.R")


###############################################################################
# STEP 2: MATCHING
###############################################################################

run_step("Step 2: 1:10 nearest-neighbour matching (04)",
         "scripts/04_matching.R")


###############################################################################
# STEP 3: PRIMARY ANALYSIS
###############################################################################

run_step("Step 3: Primary weighted regression (05)",
         "scripts/05_primary_analysis.R")


###############################################################################
# STEP 4: SUBGROUP ANALYSES
###############################################################################

run_step("Step 4: Subgroup and sensitivity analyses (06)",
         "scripts/06_subgroup_analyses.R")


###############################################################################
# STEP 5: FIGURES + TABLE 1
###############################################################################

run_step("Step 5: Figures and Table 1 (09)",
         "scripts/09_figures.R")


###############################################################################
# SUMMARY
###############################################################################

source("config.R")

cat(sprintf("\n%s\n", strrep("=", 70)))
cat("Demo pipeline complete.\n\n")

# Print primary results
primary_file <- file.path(PATHS$output_dir, "results_primary.rds")
if (file.exists(primary_file)) {
  res <- readRDS(primary_file)
  cat("Primary results:\n")
  print(
    res %>%
      dplyr::mutate(
        p_fmt = formatC(p_value, digits = 3, format = "g"),
        ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
      ) %>%
      dplyr::select(analysis, n_carriers, n_controls,
                    mean_carrier, mean_control, estimate, ci, p_fmt),
    n = Inf, width = 120
  )
}

# Print interaction results
int_file <- file.path(PATHS$output_dir, "results_interaction.rds")
if (file.exists(int_file)) {
  int_res <- readRDS(int_file)
  cat("\nDrug x carrier interaction:\n")
  print(
    int_res %>%
      dplyr::mutate(
        p_fmt = formatC(p_value, digits = 3, format = "g"),
        ci    = sprintf("[%.2f, %.2f]", ci_lo, ci_hi)
      ) %>%
      dplyr::select(analysis, interaction_term, estimate, ci, p_fmt),
    n = Inf
  )
}

# List outputs
cat("\nOutput files:\n")
list.files(PATHS$output_dir, pattern = "\\.(rds|csv)$", recursive = TRUE) %>%
  sort() %>%
  { cat(paste(" ", ., collapse = "\n"), "\n") }

cat("\nFigures:\n")
list.files(file.path(PATHS$output_dir, "figures"),
           pattern = "\\.(pdf|png|csv)$") %>%
  sort() %>%
  { cat(paste(" ", ., collapse = "\n"), "\n") }

cat(sprintf("\nAll outputs in: %s\n", normalizePath(PATHS$output_dir)))
