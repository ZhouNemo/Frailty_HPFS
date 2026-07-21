# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.4_GLME_survival.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-20 (removed quadratic model; renumbered spline models)
#
# Purpose:
#   Time-bin GLME event-study analysis with the pre-specified natural-spline
#   Gaussian LME model set (M0 raw spline, M1 primary spline, M2 full spline,
#   M3 matching-set spline) for the
#   survival-stratified risk-set matched cohorts (cases who died <=5 years after
#   diagnosis vs cases who survived >5 years), each Control vs Cancer Case.
#   Methods: Documents/Methods/TimeBin_GLME_EventStudy_Analysis.md and
#   Documents/Methods/GLME_Natural_Spline_Trajectory_Analysis.md.
#
#   NOTE: cohort membership is defined by a post-index survival outcome, so the
#   cross-cohort contrast is descriptive and conditioned on survival (the same
#   immortal-time / selective-survival caveat documented for 6.5). Post-index bins
#   are especially affected by informative death; read with that caveat.
#
#   Risk-set matching is NOT created here; run 2.4 first to build the matched data.
#
# Input:
#   Data/riskset_matched_survival_long.rds   (from 2.4_create_riskset_survival.R;
#   must pass Gate G4 and match its provenance hash)
# Outputs: 4.4_eventstudy_* and 4.4_spline_* CSV summaries in Results/cancer/data,
#   plus Results/cancer/visuals/4.4_GLME_figures.pdf containing the three figures.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")
visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")

source(file.path(project_dir, "Code", "2_data_analysis", "4.0_GLME_spline_functions.R"))

run_eventstudy_spline_analysis(
  matched_path   = file.path(data_dir, "riskset_matched_survival_long.rds"),
  results_dir    = results_dir,
  visuals_dir    = visuals_dir,
  out_prefix     = "4.4",
  analysis_title = "Frailty Event-Study by Post-Diagnosis Survival (Risk-Set Matched)",
  builder_script = "Code/2_data_analysis/2.4_create_riskset_survival.R"
)
