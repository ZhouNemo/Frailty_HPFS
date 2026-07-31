# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  5.1_ipcw_gee_overall_m2.R
# Author:  Nemo Zhou
# Date started:      2026-07-27
# Date last updated: 2026-07-28
#
# Purpose:
#   Runs the M2-equivalent full-covariate IPCW-GEE sensitivity analysis for the
#   overall incident-cancer versus risk-set matched control cohort.  The script
#   uses the Gate-G4-validated 2.5 matched cohort and the stored M2-compatible
#   df-3 spline basis from the 4.5 GLME analysis. It first fits the observation
#   model in a full-cohort person-cycle response ledger, transports predicted
#   probabilities to matched assignments, excludes death from censoring, and
#   fits only paired unweighted and IPCW-weighted GEE models with identical M2
#   fixed effects and analytic rows.
#
# Prerequisites:
#   1. Run 2.5_create_riskset_overall.R and verify Gate G4/provenance.
#   2. Run 4.5_GLME_overall.R to create current 4.5 spline metadata.
#
# Inputs:
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Data/riskset_matched_overall_long.rds
#   Results/cancer/data/4.5_spline_metadata.rds
#
# Outputs:
#   Data/riskset_matched_overall_long_ipcw_m2_gee.rds
#   Results/cancer/data/5.1_overall_m2_* CSV/RDS summaries
#
# This script does not run matching or GLME models.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

source(file.path(project_dir, "Code", "2_data_analysis", "5.0_ipcw_gee_functions.R"))

run_ipcw_m2_gee(
  panel_path = file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"),
  matched_path = file.path(data_dir, "riskset_matched_overall_long.rds"),
  spline_metadata_path = file.path(results_dir, "4.5_spline_metadata.rds"),
  cohort = "All Cancer Cohort",
  data_out_path = file.path(data_dir, "riskset_matched_overall_long_ipcw_m2_gee.rds"),
  results_dir = results_dir,
  out_prefix = "5.1_overall_m2"
)
