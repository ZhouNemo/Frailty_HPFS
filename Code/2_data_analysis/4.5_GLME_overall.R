# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.5_GLME_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-20 (removed quadratic model; renumbered spline models)
#
# Purpose:
#   OVERALL primary time-bin GLME event-study analysis with the pre-specified
#   natural-spline Gaussian LME model set (M0 raw, M1 primary spline,
#   M2 full spline, M3 matching-set spline): all incident
#   cancer cases versus risk-set matched cancer-free
#   controls, pooled in a single cohort. Estimates when frailty trajectories
#   diverge between cases and controls (before / at / after the index). Methods:
#   Documents/Methods/TimeBin_GLME_EventStudy_Analysis.md and
#   Documents/Methods/GLME_Natural_Spline_Trajectory_Analysis.md.
#
#   Risk-set matching is NOT created here; run 2.5 first to build the matched data.
#
# Input:
#   Data/riskset_matched_overall_long.rds   (from 2.5_create_riskset_overall.R)
# Outputs:
#   Results/cancer/data/4.5_eventstudy_* and 4.5_spline_* primary summaries;
#   Results/cancer/data/4.5_sensitivities/* S1-S12 results;
#   Results/cancer/visuals/4.5_GLME_overall.html and
#   Results/cancer/visuals/4.5_GLME_figures.pdf. No persistent PNG is created.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")
visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")
sensitivity_dir <- file.path(results_dir, "4.5_sensitivities")
matching_diagnostics_dir <- file.path(results_dir, "matching_diagnostics")
matched_path <- file.path(data_dir, "riskset_matched_overall_long.rds")
exact_cycle_path <- file.path(data_dir, "riskset_matched_overall_exact_cycle_long.rds")
gate_path <- file.path(matching_diagnostics_dir, "riskset_matched_overall_long_gate_g4.csv")
matching_run_path <- file.path(matching_diagnostics_dir, "riskset_matched_overall_long_run_metadata.rds")

source(file.path(project_dir, "Code", "2_data_analysis", "4.0_GLME_spline_functions.R"))

if (!file.exists(gate_path) || !file.exists(matching_run_path)) {
  stop("Gate G4 or matching provenance is missing. Rerun 2.5 before the overall GLME.")
}
gate_g4 <- read.csv(gate_path, stringsAsFactors = FALSE)
if (!("gate_pass" %in% names(gate_g4)) || !all(as.logical(gate_g4$gate_pass))) {
  stop("Gate G4 did not pass for the overall matched cohort.")
}
matching_run <- readRDS(matching_run_path)
current_matching_md5 <- unname(tools::md5sum(matched_path))
if (is.null(matching_run$output_md5) || !identical(unname(matching_run$output_md5), current_matching_md5)) {
  stop("The overall matched RDS does not match its Gate G4 provenance. Rerun 2.5.")
}

primary_result <- run_eventstudy_spline_analysis(
  matched_path   = matched_path,
  results_dir    = results_dir,
  visuals_dir    = visuals_dir,
  out_prefix     = "4.5",
  analysis_title = "Frailty Event-Study: Overall Incident Cancer vs Control (Risk-Set Matched)",
  builder_script = "Code/2_data_analysis/2.5_create_riskset_overall.R"
)

required_primary_models <- c("M1_primary_spline", "M3_primary_matching_spline")
primary_status <- primary_result$sp_status
primary_ok <- !is.null(primary_status) && all(required_primary_models %in% primary_status$model_id) &&
  all(primary_status$status[match(required_primary_models, primary_status$model_id)] %in%
        c("fit", "fit_boundary_matching_variance_zero")) &&
  all(primary_status$convergence[match(required_primary_models, primary_status$model_id)] %in% TRUE)
if (!primary_ok) {
  stop("Required M1 primary spline and M3 matching-set spline models did not both pass strict convergence. See 4.5_spline_model_status.csv.")
}

sensitivity_result <- run_overall_glme_sensitivities(
  matched_path = matched_path,
  exact_cycle_path = exact_cycle_path,
  sensitivity_dir = sensitivity_dir
)
required_sensitivities <- paste0("S", 1:12)
sensitivity_status <- bind_rows(lapply(sensitivity_result, `[[`, "status"))
sensitivity_ok <- all(required_sensitivities %in% sensitivity_status$spec) &&
  all(sensitivity_status$status[match(required_sensitivities, sensitivity_status$spec)] == "fit") &&
  all(sensitivity_status$convergence[match(required_sensitivities, sensitivity_status$spec)] %in% TRUE)
if (!sensitivity_ok) {
  stop("One or more required S1-S12 analyses failed. See 4.5_sensitivity_status.csv.")
}

run_configuration <- list(
  schema_version = 1L,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis = "overall cancer risk-set GLME",
  primary_input = normalizePath(matched_path, mustWork = TRUE),
  primary_input_md5 = unname(tools::md5sum(matched_path)),
  exact_cycle_input = if (file.exists(exact_cycle_path)) normalizePath(exact_cycle_path) else NA_character_,
  exact_cycle_input_md5 = if (file.exists(exact_cycle_path)) unname(tools::md5sum(exact_cycle_path)) else NA_character_,
  fitting_window = "all eligible relative-time observations",
  prediction_window = c(-20, 20),
  theta_windows = list(primary = c(-8, 0, 8), S5 = c(-4, 0, 4), S6 = c(-12, 0, 12)),
  covariance_fallback = c("CR2 clustered on id", "CR0 clustered on id", "model-based"),
  sensitivities = paste0("S", 1:12),
  persistent_png = FALSE
)
saveRDS(run_configuration, file.path(results_dir, "4.5_run_configuration.rds"))

report_ok <- render_overall_glme_html_report(
  results_dir = results_dir,
  sensitivity_dir = sensitivity_dir,
  visuals_dir = visuals_dir,
  report_rmd = file.path(project_dir, "Code", "2_data_analysis", "4.5_GLME_overall_report.Rmd")
)
if (!isTRUE(report_ok)) {
  stop("The required overall GLME HTML report was not rendered successfully.")
}
