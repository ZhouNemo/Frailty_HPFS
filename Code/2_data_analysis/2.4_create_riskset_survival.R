# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  2.4_create_riskset_survival.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-18 (symmetric pre-index eligibility and S7 exact-cycle rematch)
#
# Purpose:
#   Creates and overwrites the incidence-density risk-set matched cohort dataset
#   for survival-stratified cancer trajectory analyses. Incident cancer cases are
#   classified as died <=5 years after diagnosis versus survived >5 years after
#   diagnosis. Alive cases with <5 years of post-diagnosis observed follow-up are
#   excluded because their >5-year survival status is not determined.
#
# Output:
#   Data/riskset_matched_survival_long.rds
#   Data/riskset_matched_survival_exact_cycle_long.rds (S7)
# =============================================================================

library(dplyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
input_path  <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
output_path <- file.path(data_dir, "riskset_matched_survival_long.rds")
exact_cycle_output_path <- file.path(data_dir, "riskset_matched_survival_exact_cycle_long.rds")

source(file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R"))

MATCH_RATIO  <- 5
AGE_CALIPER  <- 2
CYCLE_CALIPER <- 1
MIN_VISITS   <- 1
SEED         <- 20260703
SURV_CUT     <- 5
target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")

cohort_levels <- c("Died <=5y Cohort", "Survived >5y Cohort")
classify_fn <- function(pl) {
  surv <- pl$age_at_death - pl$true_age_at_cancer
  fu_after_dx <- pl$max_age - pl$true_age_at_cancer
  died <- pl$death_status == 1

  out <- rep(NA_character_, nrow(pl))
  out[died & !is.na(surv) & surv <= SURV_CUT] <- "Died <=5y Cohort"
  out[(died & !is.na(surv) & surv > SURV_CUT) |
        (!died & !is.na(fu_after_dx) & fu_after_dx >= SURV_CUT)] <- "Survived >5y Cohort"
  out
}

matched_result <- build_riskset_matched_long(
  input_path = input_path,
  classification_vars = character(0),
  classify_fn = classify_fn,
  cohort_levels = cohort_levels,
  target_cycles = target_cycles,
  match_ratio = MATCH_RATIO,
  age_caliper = AGE_CALIPER,
  cycle_caliper = CYCLE_CALIPER,
  min_visits = MIN_VISITS,
  seed = SEED
)

save_riskset_match(matched_result, output_path, "Survival-stratified risk-set cohort")

# S7 exact-cycle rematch; retain the primary case-based age scale.
exact_cycle_result <- build_riskset_matched_long(
  input_path = input_path,
  classification_vars = character(0),
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
  "Survival-stratified exact-cycle risk-set cohort (S7)"
)
