# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script: 2.0_matching_provenance.R
# Author: Nemo Zhou
# Date started: 2026-09-28
# Date last updated: 2026-09-29 (add individual cancer-site builder provenance)
# Purpose: Read-only validation of outcome-independent eligibility, Gate G4,
# and exact outcome/assignment-ledger hashes. Shared by downstream analyses and reports.
# Legacy and retired cohorts fail with the matching builder Nemo must rerun.
# No matching, model fitting, dataset writes, or report generation occurs here.
# =============================================================================

matching_builder_for <- function(matched_path) {
  builders <- c(
    riskset_matched_analysis_long = "2.1_create_riskset_high_low_burden.R",
    riskset_matched_smoking_long = "2.2_create_riskset_smoking_related.R",
    riskset_matched_obesity_long = "2.3_create_riskset_obesity_related.R",
    riskset_matched_survival_long = "2.4_create_riskset_survival.R",
    riskset_matched_overall_long = "2.5_create_riskset_overall.R",
    riskset_matched_site_lung_long = "2.6_create_riskset_individual_types.R",
    riskset_matched_site_colorectal_long = "2.6_create_riskset_individual_types.R",
    riskset_matched_site_prostate_long = "2.6_create_riskset_individual_types.R",
    riskset_matched_site_bladder_long = "2.6_create_riskset_individual_types.R",
    riskset_matched_site_pancreas_long = "2.6_create_riskset_individual_types.R",
    riskset_matched_site_kidney_long = "2.6_create_riskset_individual_types.R")
  stem <- tools::file_path_sans_ext(basename(matched_path))
  if (stem %in% names(builders)) file.path("Code/2_data_analysis", builders[[stem]])
  else "the corresponding Code/2_data_analysis/2.x matching builder"
}

matching_assignment_path <- function(matched_path) {
  if (!grepl("_long[.]rds$", matched_path))
    stop("Expected a primary *_long.rds path for assignment provenance.", call. = FALSE)
  sub("_long[.]rds$", "_assignments.rds", matched_path)
}

assert_matching_metadata <- function(metadata, matched_path) {
  valid <- identical(metadata$eligibility_version, "no_fi_requirement_v1") &&
    identical(metadata$fi_requirement, "none") &&
    identical(metadata$cohort_entry_rule, "first_participated_analytic_cycle_return") &&
    identical(metadata$control_entry_on_or_before_index, TRUE)
  if (!valid) stop("Matched cohort uses legacy or incompatible eligibility; require no_fi_requirement_v1. ",
    "Nemo must run: Rscript ", matching_builder_for(matched_path), call. = FALSE)
  invisible(TRUE)
}

validate_assignment_provenance <- function(metadata, matched_path) {
  path <- matching_assignment_path(matched_path)
  if (!file.exists(path)) stop("Assignment ledger missing. Nemo must run: Rscript ",
                              matching_builder_for(matched_path), call. = FALSE)
  expected <- unname(as.character(metadata$assignment_md5))
  observed <- unname(tools::md5sum(path))
  if (length(expected) != 1L || is.na(expected) || !identical(expected, observed))
    stop("Assignment ledger hash does not match matching provenance. Nemo must run: Rscript ",
         matching_builder_for(matched_path), call. = FALSE)
  list(assignment_path = path, assignment_md5 = observed)
}

assert_assignment_rows <- function(data, label = "Matched input") {
  required <- c("eligibility_version", "trajectory_id", "first_return", "index_date")
  if (!all(required %in% names(data)))
    stop(label, " predates no-FI eligibility. Rebuild matching and regenerate this artifact.", call. = FALSE)
  valid <- data$eligibility_version == "no_fi_requirement_v1" &
    is.finite(data$first_return) & is.finite(data$index_date) &
    data$first_return <= data$index_date
  if (!all(valid %in% TRUE)) stop(label, " fails cohort-entry/version checks.", call. = FALSE)
  invisible(TRUE)
}

read_matching_assignments <- function(matched_path) {
  provenance <- validate_matching_provenance(matched_path)
  ledger <- readRDS(provenance$assignment_path)
  assert_assignment_rows(ledger, "Assignment ledger")
  if (anyDuplicated(ledger$trajectory_id)) stop("Duplicate assignment ledger keys.")
  ledger
}

# Preserve every matched assignment when adding observed-outcome summaries.
# Baseline/index values come from the ledger; absent observation counts are zero.
complete_assignment_summary <- function(ledger, summary, keys, count_columns) {
  if (anyDuplicated(ledger[keys]) || anyDuplicated(summary[keys]))
    stop("Duplicate keys in assignment-level summaries.")
  ledger <- dplyr::mutate(ledger, dplyr::across(dplyr::all_of(keys), as.character))
  summary <- dplyr::mutate(summary, dplyr::across(dplyr::all_of(keys), as.character))
  extra <- setdiff(names(summary), names(ledger))
  out <- dplyr::left_join(ledger, summary[unique(c(keys, extra))], by = keys)
  for (v in intersect(count_columns, names(out))) out[[v]][is.na(out[[v]])] <- 0L
  out
}

validate_matching_provenance <- function(matched_path) {
  builder <- matching_builder_for(matched_path)
  if (grepl("exact_cycle|cancer_free_full_endpoint", basename(matched_path)))
    stop("This matched design is retired; historical outputs are not active analysis inputs.", call. = FALSE)
  if (!file.exists(matched_path)) {
    stop("Matched dataset not found at ", matched_path, call. = FALSE)
  }
  project_dir <- dirname(dirname(normalizePath(matched_path, mustWork = FALSE)))
  diagnostics_dir <- file.path(project_dir, "Results", "cancer", "data",
                                "matching_diagnostics")
  output_stem <- tools::file_path_sans_ext(basename(matched_path))
  gate_path <- file.path(diagnostics_dir, paste0(output_stem, "_gate_g4.csv"))
  run_path <- file.path(diagnostics_dir, paste0(output_stem, "_run_metadata.rds"))
  if (!file.exists(gate_path) || !file.exists(run_path)) {
    stop("Gate G4 or matching provenance is missing for ", output_stem,
         ". Nemo must run: Rscript ", builder,
         call. = FALSE)
  }
  gate <- tryCatch(read.csv(gate_path, stringsAsFactors = FALSE),
                   error = function(e) {
                     stop("Could not read matching Gate G4 file ", gate_path,
                          ": ", conditionMessage(e), call. = FALSE)
                   })
  if (!nrow(gate) || !("gate_pass" %in% names(gate)) ||
      !all(as.logical(gate$gate_pass) %in% TRUE)) {
    stop("Gate G4 did not pass for ", output_stem,
         ". Refusing to fit the matched GLME input.", call. = FALSE)
  }
  run_metadata <- readRDS(run_path)
  assert_matching_metadata(run_metadata, matched_path)
  expected_md5 <- if (!is.null(run_metadata$output_md5)) {
    unname(as.character(run_metadata$output_md5))
  } else {
    character(0)
  }
  observed_md5 <- unname(tools::md5sum(matched_path))
  if (length(expected_md5) != 1L || !nzchar(expected_md5) ||
      !identical(expected_md5, observed_md5)) {
    stop("Matched RDS hash does not match its Gate G4 provenance for ",
         output_stem, ". Nemo must run: Rscript ", builder, call. = FALSE)
  }
  assignment <- validate_assignment_provenance(run_metadata, matched_path)
  list(output_stem = output_stem, gate = gate, run_metadata = run_metadata,
       gate_path = gate_path, run_path = run_path, input_md5 = observed_md5,
       assignment_path = assignment$assignment_path, assignment_md5 = assignment$assignment_md5)
}

# Derived weighted files must retain the exact primary inputs they came from.
validate_derived_assignment_provenance <- function(data, label = "Derived input") {
  assert_assignment_rows(data, label)
  recorded <- attr(data, "assignment_provenance")
  if (is.null(recorded$assignment_path)) stop(label, " lacks assignment provenance; rerun its 7.x preparation.")
  primary <- sub("_assignments[.]rds$", "_long.rds", recorded$assignment_path)
  current <- validate_matching_provenance(primary)
  if (!identical(recorded$assignment_md5, current$assignment_md5) ||
      !identical(recorded$input_md5, current$input_md5))
    stop(label, " has stale primary ledger/outcome hashes; rerun its 7.x preparation.")
  invisible(current)
}
