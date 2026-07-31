# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  5.2_ipcw_gee_subtype_m2.R
# Author:  Nemo Zhou
# Date started:      2026-07-27
# Date last updated: 2026-07-28
#
# Purpose:
#   Runs the M2-equivalent full-covariate IPCW-GEE sensitivity analysis for the
#   active non-survival cancer cohort families: low/moderate- and high-burden,
#   smoking-related, and obesity-related cancers.  Each cohort uses its own
#   Gate-G4-validated matched input and its own stored 4.x df-3 spline metadata.
#   The shared engine fits one full-cohort response model, transports its
#   probabilities to assignment-specific four-year IPCW ledgers, treats death
#   as a terminal event rather than censoring, and fits paired unweighted and
#   IPCW-weighted M2 GEE models only.
#
# Prerequisites:
#   1. Run 2.1, 2.2, and 2.3 and verify Gate G4/provenance for each input.
#   2. Run 4.1, 4.2, and 4.3 to create current spline metadata.
#
# Inputs:
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Data/riskset_matched_analysis_long.rds
#   Data/riskset_matched_smoking_long.rds
#   Data/riskset_matched_obesity_long.rds
#   Results/cancer/data/4.1_spline_metadata.rds
#   Results/cancer/data/4.2_spline_metadata.rds
#   Results/cancer/data/4.3_spline_metadata.rds
#
# Outputs:
#   Data/riskset_matched_*_long_ipcw_m2_gee.rds
#   Results/cancer/data/5.2_*_m2_* CSV/RDS summaries
#
# This script does not run matching or GLME models.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

source(file.path(project_dir, "Code", "2_data_analysis", "5.0_ipcw_gee_functions.R"))

# Fit the transported full-cohort response model once, then predict its
# observation probabilities for each matched subtype assignment ledger.
full_response_model <- ipcw_gee_fit_full_response_model(
  file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
)

cohort_specs <- list(
  list(
    cohort = "Low/Moderate Burden Cohort",
    matched_path = file.path(data_dir, "riskset_matched_analysis_long.rds"),
    spline_metadata_path = file.path(results_dir, "4.1_spline_metadata.rds"),
    data_out_path = file.path(data_dir, "riskset_matched_low_moderate_burden_long_ipcw_m2_gee.rds"),
    out_prefix = "5.2_low_moderate_burden_m2"
  ),
  list(
    cohort = "High Burden Cohort",
    matched_path = file.path(data_dir, "riskset_matched_analysis_long.rds"),
    spline_metadata_path = file.path(results_dir, "4.1_spline_metadata.rds"),
    data_out_path = file.path(data_dir, "riskset_matched_high_burden_long_ipcw_m2_gee.rds"),
    out_prefix = "5.2_high_burden_m2"
  ),
  list(
    cohort = "Smoking-Related Cancer Cohort",
    matched_path = file.path(data_dir, "riskset_matched_smoking_long.rds"),
    spline_metadata_path = file.path(results_dir, "4.2_spline_metadata.rds"),
    data_out_path = file.path(data_dir, "riskset_matched_smoking_long_ipcw_m2_gee.rds"),
    out_prefix = "5.2_smoking_m2"
  ),
  list(
    cohort = "Obesity-Related Cancer Cohort",
    matched_path = file.path(data_dir, "riskset_matched_obesity_long.rds"),
    spline_metadata_path = file.path(results_dir, "4.3_spline_metadata.rds"),
    data_out_path = file.path(data_dir, "riskset_matched_obesity_long_ipcw_m2_gee.rds"),
    out_prefix = "5.2_obesity_m2"
  )
)

for (spec in cohort_specs) {
  run_ipcw_m2_gee(
    panel_path = file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"),
    matched_path = spec$matched_path,
    spline_metadata_path = spec$spline_metadata_path,
    cohort = spec$cohort,
    data_out_path = spec$data_out_path,
    results_dir = results_dir,
    out_prefix = spec$out_prefix,
    response_model = full_response_model
  )
}
