# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  6.4_riskset_smoking_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-06-28
# Date last updated: 2026-06-29
#
# Purpose:
#   Fits quadratic frailty trajectories for smoking-related versus non-smoking-
#   related risk-set matched cancer cohorts. Risk-set matching is no longer
#   created in this script; run 2.2 first to overwrite the matched cohort
#   dataset.
#
# Input:
#   Data/riskset_matched_smoking_long.rds
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer")

source(file.path(project_dir, "Code", "2_data_analysis", "6.0_quadratic_trajectory_functions.R"))

run_quadratic_trajectory_analysis(
  matched_path = file.path(data_dir, "riskset_matched_smoking_long.rds"),
  results_dir = results_dir,
  analysis_title = "Frailty Trajectories: Smoking-Related vs Non-Smoking Cancer (Risk-Set Matched)",
  out_fig = "6.4_smoking_quadratic_trajectories.png",
  out_coef = "6.4_smoking_fixed_effects_CI.csv",
  out_pred = "6.4_smoking_predicted_trajectories.csv",
  builder_script = "Code/2_data_analysis/2.2_create_riskset_smoking_related.R"
)
