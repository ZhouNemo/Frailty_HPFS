# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  6.3_riskset_obesity_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-06-29
#
# Purpose:
#   Fits quadratic frailty trajectories for obesity-related versus non-obesity-
#   related risk-set matched cancer cohorts. Risk-set matching is no longer
#   created in this script; run 2.3 first to overwrite the matched cohort
#   dataset.
#
# Input:
#   Data/riskset_matched_obesity_long.rds
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer")

source(file.path(project_dir, "Code", "2_data_analysis", "6.0_quadratic_trajectory_functions.R"))

run_quadratic_trajectory_analysis(
  matched_path = file.path(data_dir, "riskset_matched_obesity_long.rds"),
  results_dir = results_dir,
  analysis_title = "Frailty Trajectories: Obesity-Related vs Non-Obesity Cancer (Risk-Set Matched)",
  out_fig = "6.3_obesity_quadratic_trajectories.png",
  out_coef = "6.3_obesity_fixed_effects_CI.csv",
  out_pred = "6.3_obesity_predicted_trajectories.csv",
  builder_script = "Code/2_data_analysis/2.3_create_riskset_obesity_related.R"
)
