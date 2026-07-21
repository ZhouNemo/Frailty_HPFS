# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  6.5_riskset_survival_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-06-28
# Date last updated: 2026-06-29
#
# Purpose:
#   Fits quadratic frailty trajectories for survival-stratified risk-set matched
#   cancer cohorts: died <=5 years after diagnosis versus survived >5 years.
#   Risk-set matching is no longer created in this script; run 2.4 first to
#   overwrite the matched cohort dataset.
#
# Important caveat:
#   Cohort membership is defined by a post-index outcome, so the cross-cohort
#   contrast is descriptive and conditioned on survival.
#
# Input:
#   Data/riskset_matched_survival_long.rds
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer")

source(file.path(project_dir, "Code", "2_data_analysis", "6.0_quadratic_trajectory_functions.R"))

run_quadratic_trajectory_analysis(
  matched_path = file.path(data_dir, "riskset_matched_survival_long.rds"),
  results_dir = results_dir,
  analysis_title = "Frailty Trajectories by Survival After Diagnosis (Risk-Set Matched)",
  out_fig = "6.5_survival_quadratic_trajectories.png",
  out_coef = "6.5_survival_fixed_effects_CI.csv",
  out_pred = "6.5_survival_predicted_trajectories.csv",
  builder_script = "Code/2_data_analysis/2.4_create_riskset_survival.R"
)
