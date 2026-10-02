# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  2.5_create_riskset_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
#
# Purpose:
#   Creates and overwrites the incidence-density risk-set matched cohort dataset
#   for the OVERALL primary comparison: all incident (first) cancer cases versus
#   risk-set matched cancer-free controls, pooled into a single cohort with no
#   subgroup split. This is the headline "any incident cancer vs control"
#   analysis; the subgroup builders 2.1-2.4 (burden, smoking, obesity, survival)
#   are stratified versions of the same design.
#
#   Each incident case defines a real index date; controls alive, under
#   observation, and cancer-free at that index are matched within the age caliper
#   and inherit the case's index date (see
#   Documents/Methods/Risk_Set_Index_Time_Assignment_for_Controls.md).
#
# Eligibility: analytic-cohort entry by index; no FI availability or
# questionnaire-recency requirement. A complete assignment ledger is also saved.
#
# Output:
#   Data/riskset_matched_overall_long.rds
#   Data/riskset_matched_overall_assignments.rds
# =============================================================================

library(dplyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
input_path  <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
output_path <- file.path(data_dir, "riskset_matched_overall_long.rds")

source(file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R"))

MATCH_RATIO  <- 5
AGE_CALIPER  <- 2
SEED         <- 20260703
target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")

# Single pooled cohort: every incident case is assigned the same cohort label.
# Non-cases are dropped downstream by the is_case filter inside the builder.
cohort_levels <- c("All Cancer Cohort")
classify_fn <- function(pl) rep("All Cancer Cohort", nrow(pl))

matched_result <- build_riskset_matched_long(
  input_path = input_path,
  classification_vars = character(0),
  classify_fn = classify_fn,
  cohort_levels = cohort_levels,
  target_cycles = target_cycles,
  match_ratio = MATCH_RATIO,
  age_caliper = AGE_CALIPER,
  seed = SEED
)

save_riskset_match(matched_result, output_path, "Overall incident-cancer risk-set cohort")
