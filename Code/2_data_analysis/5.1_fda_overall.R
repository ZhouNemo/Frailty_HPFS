# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  5.1_fda_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-17 (visual-data outputs moved; PNG writes removed)
#
# Purpose:
#   Converts the OVERALL risk-set matched dataset (any incident cancer vs control)
#   into the sparse functional-data format required by fdapace, and runs sparse
#   FPCA with a Cancer-vs-Control mean-function comparison. Functional unit = one
#   curve per trajectory assignment (Cohort__match_set__id__role).
#
#   Risk-set matching is NOT created here; run 2.5 first.
#   The same prep/FPCA functions work on any 2.x matched dataset -- just point
#   `matched_path` at riskset_matched_{smoking,obesity,survival,analysis}_long.rds.
#
# Input:
#   Data/riskset_matched_overall_long.rds      (from 2.5_create_riskset_overall.R)
# Outputs:
#   Data/fda_input_overall.rds                 (Ly/Lt/meta list for fdapace)
#   Results/cancer/data/5.1_fpca_scores.csv
#   Results/cancer/data/5.1_fpca_eigenfunctions.csv
#   Results/cancer/data/5.1_fpca_mean_by_group.csv
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

source(file.path(project_dir, "Code", "2_data_analysis", "5.0_fda_functions.R"))

# 1) Convert to fdapace Lt/Ly/meta and cache the functional-data object.
fda_input <- prepare_fda_input(
  matched_path   = file.path(data_dir, "riskset_matched_overall_long.rds"),
  out_path       = file.path(data_dir, "fda_input_overall.rds"),
  min_measures   = 3,            # curves need >= 3 points for sparse FPCA
  window_yrs     = 16,           # keep within the supported relative-time window
  builder_script = "Code/2_data_analysis/2.5_create_riskset_overall.R"
)

# 2) Run sparse FPCA + Cancer-vs-Control mean comparison (needs fdapace).
#    Controls dominate; cap curves per group to keep PACE tractable. Raise/remove
#    max_curves_per_group once you confirm runtime/memory are acceptable.
run_fpca(
  fda_input,
  results_dir          = results_dir,
  out_prefix           = "5.1",
  fve_threshold        = 0.95,
  max_curves_per_group = 6000
)
