# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.3_GLME_obesity_related.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-20 (removed quadratic model; renumbered spline models)
#
# Purpose:
#   Time-bin GLME event-study analysis with the pre-specified natural-spline
#   Gaussian LME model set (M0 raw spline, M1 primary spline, M2 full spline,
#   M3 matching-set spline) for the
#   obesity-related cancer versus cancer-free control risk-set matched cohort.
#   Methods: Documents/Methods/TimeBin_GLME_EventStudy_Analysis.md and
#   Documents/Methods/GLME_Natural_Spline_Trajectory_Analysis.md.
#
#   Risk-set matching is NOT created here; run 2.3 first to build the matched data.
#
# Input:
#   Data/riskset_matched_obesity_long.rds   (from 2.3_create_riskset_obesity_related.R;
#   must pass Gate G4 and match its provenance hash)
# Outputs: 4.3_eventstudy_* and 4.3_spline_* CSV summaries in Results/cancer/data,
#   plus Results/cancer/visuals/4.3_GLME_figures.pdf containing the three figures.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")
visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")

source(file.path(project_dir, "Code", "2_data_analysis", "4.0_GLME_spline_functions.R"))

run_eventstudy_spline_analysis(
  matched_path   = file.path(data_dir, "riskset_matched_obesity_long.rds"),
  results_dir    = results_dir,
  visuals_dir    = visuals_dir,
  out_prefix     = "4.3",
  analysis_title = "Frailty Event-Study by Obesity-Related Cancer (Risk-Set Matched)",
  builder_script = "Code/2_data_analysis/2.3_create_riskset_obesity_related.R"
)
