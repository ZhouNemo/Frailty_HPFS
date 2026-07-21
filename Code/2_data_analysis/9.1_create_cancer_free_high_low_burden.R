# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  9.1_create_cancer_free_high_low_burden.R
# Author:  Nemo Zhou
# Date started:      2026-07-19
# Date last updated: 2026-07-20 (canonical earliest-cancer index date)
#
# Purpose:
#   Creates a separate full-HPFS-endpoint cancer-free-through-follow-up
#   incidence-density risk-set matched cohort for canonical high-burden and
#   low/moderate-burden cancer cases. Controls remain subject to the existing
#   alive-at-index, active-follow-up, age-caliper, cycle-caliper, and pre-index
#   FI-support requirements, with the additional restriction that their
#   cancer_index_dateca is missing or later than the maximum cancer ascertainment
#   month in the complete input dataset.
#
#   This workflow is intentionally separate from 2.1 and does not overwrite
#   primary 2.x matching datasets or any 3.xx report inputs.
#
# Output:
#   Data/riskset_matched_analysis_cancer_free_full_endpoint_long.rds
#   Data/riskset_matched_analysis_cancer_free_full_endpoint_exact_cycle_long.rds
#   Matching diagnostics are written under
#   Results/cancer/data/matching_diagnostics with unique output prefixes.
# =============================================================================

library(dplyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
input_path  <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
output_path <- file.path(
  data_dir,
  "riskset_matched_analysis_cancer_free_full_endpoint_long.rds"
)
exact_cycle_output_path <- file.path(
  data_dir,
  "riskset_matched_analysis_cancer_free_full_endpoint_exact_cycle_long.rds"
)

source(file.path(
  project_dir,
  "Code",
  "2_data_analysis",
  "9.0_cancer_free_matching_functions.R"
))

MATCH_RATIO   <- 5
AGE_CALIPER   <- 2
CYCLE_CALIPER <- 1
MIN_VISITS    <- 1
SEED          <- 20260703
target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")

cohort_levels <- c("Low/Moderate Burden Cohort", "High Burden Cohort")
classify_fn <- function(pl) {
  ifelse(
    coalesce0(pl$is_high_burden) == 1,
    "High Burden Cohort",
    "Low/Moderate Burden Cohort"
  )
}

matched_result <- build_riskset_matched_long(
  input_path = input_path,
  classification_vars = "is_high_burden",
  classify_fn = classify_fn,
  cohort_levels = cohort_levels,
  target_cycles = target_cycles,
  match_ratio = MATCH_RATIO,
  age_caliper = AGE_CALIPER,
  cycle_caliper = CYCLE_CALIPER,
  min_visits = MIN_VISITS,
  seed = SEED
)

save_riskset_match(
  matched_result,
  output_path,
  "High/low burden cancer-free-through-full-endpoint cohort"
)

# S7: rematch with an exact active-cycle requirement, reusing the primary
# full-endpoint cancer-free case-based age standardization constants.
exact_cycle_result <- build_riskset_matched_long(
  input_path = input_path,
  classification_vars = "is_high_burden",
  classify_fn = classify_fn,
  cohort_levels = cohort_levels,
  target_cycles = target_cycles,
  match_ratio = MATCH_RATIO,
  age_caliper = AGE_CALIPER,
  cycle_caliper = 0,
  min_visits = MIN_VISITS,
  seed = SEED,
  index_age_scaling = matched_result$scaling_metadata
)

save_riskset_match(
  exact_cycle_result,
  exact_cycle_output_path,
  "High/low burden cancer-free exact-cycle cohort (S7)"
)
