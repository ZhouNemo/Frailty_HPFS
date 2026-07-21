# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.1_ipcw_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-17 (visual-data outputs moved; PNG writes removed)
#
# Purpose:
#   Builds stabilized inverse-probability-of-censoring weights (IPCW) for the
#   OVERALL any-cancer-vs-control comparison and attaches post-index conditional
#   weights to the overall matched dataset, so the event-study (4.5) and quadratic
#   (6.x) models can be re-fit with weights = sw_ipcw to address informative loss
#   to follow-up / questionnaire nonresponse after the index.
#
#   The censoring model and weight construction are documented in
#   7.0_ipcw_functions.R and Documents/Methods/IPCW_Censoring_Weights.md.
#   Death is NOT weighted away (handled by the estimand, not as recoverable
#   censoring). Covariates (lagged): race, marital_status, living_arr, smoke,
#   pckyr, bmi, act, alco (section 10), plus prior FI, cancer status/timing, age,
#   cycle, and response history. No time-varying comorbidity items.
#
#   Run 2.5 first (builds the overall matched dataset). The same functions apply
#   to any 2.x dataset by changing matched_path / out_path.
#
# Inputs:
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds   (longitudinal panel)
#   Data/riskset_matched_overall_long.rds               (from 2.5)
# Outputs:
#   Data/riskset_matched_overall_long_ipcw.rds          (matched long + sw_ipcw)
#   Results/cancer/data/7.1_ipcw_weight_summary.csv
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

source(file.path(project_dir, "Code", "2_data_analysis", "7.0_ipcw_functions.R"))

# Same 4-year analytic cycles used by the matching builders (2.x) and models.
target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")

run_ipcw(
  panel_path    = file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"),
  matched_path  = file.path(data_dir, "riskset_matched_overall_long.rds"),
  out_path      = file.path(data_dir, "riskset_matched_overall_long_ipcw.rds"),
  target_cycles = target_cycles,
  results_dir   = results_dir,
  out_prefix    = "7.1",
  age_df        = 4,
  trunc         = c(0.01, 0.99)
)

# -----------------------------------------------------------------------------
# Using the weights downstream (example): re-fit the event-study weighted.
#   ml <- readRDS(file.path(data_dir, "riskset_matched_overall_long_ipcw.rds"))
#   library(lme4)
#   m <- lmer(fi_score_nocancer ~ Group * rel_time_bin + index_age_z +
#             base_race + base_marital + base_pckgr + (1 | id),
#             data = ml, weights = sw_ipcw, REML = TRUE,
#             control = lmerControl(optimizer = "bobyqa"))
# A weighted variant of 4.0's engine can read this dataset and pass weights = sw_ipcw.
# -----------------------------------------------------------------------------
