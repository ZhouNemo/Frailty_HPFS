# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.3_ipcw_subtype_cohorts.R
# Author:  Nemo Zhou
# Date started:      2026-06-30
# Date last updated: 2026-07-17 (subtype cohorts now related cancer vs controls)
#
# Purpose:
#   Builds one shared set of stabilized inverse-probability-of-censoring weights
#   (IPCW) from the HPFS longitudinal panel and attaches post-index conditional
#   weights (`sw_ipcw`) to the active cancer-subtype risk-set matched cohorts:
#     * high-burden versus low/moderate-burden cancer,
#     * smoking-related cancer versus cancer-free controls, and
#     * obesity-related cancer versus cancer-free controls.
#
#   Prostate cancer subtype cohorts are deliberately not included here. Survival-
#   stratified cohorts are also not included because they classify cases by a
#   post-index prognosis/survival outcome rather than cancer subtype.
#
#   The censoring model and weight construction are documented in
#   7.0_ipcw_functions.R and Documents/Methods/IPCW_Censoring_Weights.md.
#   Death is not weighted away; the weights address loss to follow-up /
#   questionnaire nonresponse among survivors.
#
# Inputs:
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Data/riskset_matched_analysis_long.rds
#   Data/riskset_matched_smoking_long.rds
#   Data/riskset_matched_obesity_long.rds
# Outputs:
#   Data/riskset_matched_analysis_long_ipcw.rds
#   Data/riskset_matched_smoking_long_ipcw.rds
#   Data/riskset_matched_obesity_long_ipcw.rds
#   Results/cancer/data/7.3_*_ipcw_weight_summary.csv
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

source(file.path(project_dir, "Code", "2_data_analysis", "7.0_ipcw_functions.R"))

# Same 4-year analytic cycles used by the matching builders (2.x) and models.
target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")

ipcw_fit <- build_ipcw_weights(
  panel_path = file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"),
  target_cycles = target_cycles,
  age_df = 4
)

subtype_specs <- list(
  burden = list(
    matched_path = file.path(data_dir, "riskset_matched_analysis_long.rds"),
    out_path = file.path(data_dir, "riskset_matched_analysis_long_ipcw.rds"),
    out_prefix = "7.3_burden"
  ),
  smoking = list(
    matched_path = file.path(data_dir, "riskset_matched_smoking_long.rds"),
    out_path = file.path(data_dir, "riskset_matched_smoking_long_ipcw.rds"),
    out_prefix = "7.3_smoking"
  ),
  obesity = list(
    matched_path = file.path(data_dir, "riskset_matched_obesity_long.rds"),
    out_path = file.path(data_dir, "riskset_matched_obesity_long_ipcw.rds"),
    out_prefix = "7.3_obesity"
  )
)

for (nm in names(subtype_specs)) {
  spec <- subtype_specs[[nm]]
  cat("\n======== Attaching subtype IPCW:", nm, "========\n")
  attach_ipcw_to_matched(
    matched_path = spec$matched_path,
    idcycle = ipcw_fit$idcycle,
    out_path = spec$out_path,
    results_dir = results_dir,
    out_prefix = spec$out_prefix,
    trunc = c(0.01, 0.99)
  )
}
