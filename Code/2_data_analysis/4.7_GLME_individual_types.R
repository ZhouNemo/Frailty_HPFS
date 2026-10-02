# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.7_GLME_individual_types.R
# Author:  Nemo Zhou
# Date started:      2026-09-29
# Date last updated: 2026-09-29
#
# Purpose:
#   Fit the existing M0--M3 natural-spline Gaussian mixed-effects model set
#   separately for the six incident first-cancer site cohorts built by 2.6.
#   Each site uses its own risk-set controls, assignment ledger, Gate G4, and
#   result prefix. --site=<slug> runs one site; the default runs all six.
#   M1 must fit and converge strictly. A failed site is recorded, never
#   presented as fitted, and does not prevent attempts for remaining sites.
#
# Inputs:  Data/riskset_matched_site_<slug>_long.rds and companion provenance.
# Outputs: Results/cancer/data/4.7_<slug>_* and
#          Results/cancer/data/4.7_individual_type_run_status.csv.
#          No persistent PNG files are produced.
# Run only after Nemo has run 2.6; Nemo runs this long modeling script.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
source(file.path(project_dir, "Code", "2_data_analysis", "2.6_create_riskset_individual_types.R"))
source(file.path(project_dir, "Code", "2_data_analysis", "4.0_GLME_spline_functions.R"))

verify_site_matching <- function(site, matched_path, status_path, input_md5,
                                 builder_md5, matching_engine_md5) {
  if (!file.exists(status_path)) stop("2.6 matching status is missing.", call. = FALSE)
  status <- read.csv(status_path, stringsAsFactors = FALSE)
  row <- status[status$site == site, , drop = FALSE]
  if (nrow(row) != 1L || !identical(row$status[[1]], "matched"))
    stop("Site ", site, " has no successful current 2.6 run.", call. = FALSE)
  if (!file.exists(matched_path)) stop("Matched site input is missing.", call. = FALSE)
  matched_md5 <- unname(tools::md5sum(matched_path))
  if (!identical(row$input_md5[[1]], input_md5) ||
      !identical(row$builder_md5[[1]], builder_md5) ||
      !identical(row$matching_engine_md5[[1]], matching_engine_md5) ||
      !identical(row$output_md5[[1]], matched_md5))
    stop("Site matching status is stale or its matched RDS changed. Rerun 2.6 for ", site,
         ".", call. = FALSE)
  provenance <- validate_matching_provenance(matched_path)
  meta <- provenance$run_metadata
  expected_codes <- as.integer(site_icd_codes[[site]])
  valid <- identical(meta$site, site) &&
    identical(as.integer(meta$site_icds), expected_codes) &&
    identical(meta$site_label, unname(site_labels[[site]])) &&
    identical(meta$case_definition, "canonical earliest dated cancer ICDs only") &&
    identical(meta$tie_policy, "include in each applicable site cohort") &&
    identical(meta$input_md5, input_md5) &&
    identical(meta$builder_md5, builder_md5) &&
    identical(meta$matching_engine_md5, matching_engine_md5) &&
    identical(provenance$input_md5, matched_md5)
  if (!valid) stop("Site definition or source provenance mismatch for ", site,
                   ". Rerun 2.6 for this site.", call. = FALSE)
  if (nrow(provenance$gate) != 1L || !isTRUE(as.logical(provenance$gate$gate_pass[[1]])))
    stop("Gate G4 did not pass for ", site, ".", call. = FALSE)
  provenance
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  sites <- requested_sites(args)
  data_dir <- file.path(project_dir, "Data")
  results_dir <- file.path(project_dir, "Results", "cancer", "data")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")
  matching_status <- file.path(results_dir, "matching_diagnostics", "2.6_individual_type_run_status.csv")
  status_path <- file.path(results_dir, "4.7_individual_type_run_status.csv")
  input_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
  builder_path <- file.path(project_dir, "Code", "2_data_analysis", "2.6_create_riskset_individual_types.R")
  matching_engine_path <- file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R")
  glme_engine_path <- file.path(project_dir, "Code", "2_data_analysis", "4.0_GLME_spline_functions.R")
  input_md5 <- unname(tools::md5sum(input_path))
  builder_md5 <- unname(tools::md5sum(builder_path))
  matching_engine_md5 <- unname(tools::md5sum(matching_engine_path))
  glme_engine_md5 <- unname(tools::md5sum(glme_engine_path))

  write_status <- function(site, status, matched_md5 = NA_character_, m1_status = NA_character_, message = "") {
    previous <- if (file.exists(status_path)) read.csv(status_path, stringsAsFactors = FALSE) else data.frame()
    if (nrow(previous)) previous <- previous[previous$site != site, , drop = FALSE]
    row <- data.frame(
      site = site, status = status, input_md5 = input_md5,
      matched_md5 = matched_md5, builder_md5 = builder_md5,
      matching_engine_md5 = matching_engine_md5, glme_engine_md5 = glme_engine_md5,
      m1_status = m1_status, message = message,
      updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      stringsAsFactors = FALSE
    )
    write.csv(rbind(previous, row), status_path, row.names = FALSE)
  }
  for (site in sites) write_status(site, "pending")
  failures <- character()
  for (site in sites) {
    write_status(site, "running")
    tryCatch({
      matched_path <- file.path(data_dir, paste0("riskset_matched_site_", site, "_long.rds"))
      provenance <- verify_site_matching(site, matched_path, matching_status,
                                         input_md5, builder_md5, matching_engine_md5)
      prefix <- paste0("4.7_", site)
      result <- run_spline_analysis(
        matched_path = matched_path,
        results_dir = results_dir, visuals_dir = visuals_dir,
        out_prefix = prefix,
        analysis_title = paste("Frailty Trajectories:", unname(site_labels[[site]]),
                               "vs Risk-Set Controls"),
        builder_script = "Code/2_data_analysis/2.6_create_riskset_individual_types.R"
      )
      m1 <- result$sp_status[result$sp_status$model_id == "M1_primary_spline" &
                               as.character(result$sp_status$Cohort) == unname(site_labels[[site]]),,
                             drop = FALSE]
      if (nrow(m1) != 1L ||
          !(as.character(m1$status[[1]]) %in% c("fit", "fit_boundary_matching_variance_zero")) ||
          !isTRUE(as.logical(m1$convergence[[1]])))
        stop("Required M1 did not fit and converge for ", site,
             "; inspect ", prefix, "_spline_model_status.csv.", call. = FALSE)
      write_status(site, "fit", provenance$input_md5, as.character(m1$status[[1]]))
    }, error = function(e) {
      failures <<- c(failures, site)
      write_status(site, "failed", message = conditionMessage(e))
      message("Site ", site, " failed: ", conditionMessage(e))
    })
  }
  if (length(failures)) stop("GLME failed for: ", paste(failures, collapse = ", "),
                             ". See ", status_path, call. = FALSE)
  invisible(status_path)
}

if (sys.nframe() == 0L) main()
