# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.0.2_test_descriptive_contracts.R
# Author:  Nemo Zhou
# Date started: 2026-09-21
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
# Purpose: Model-free contracts for the observed descriptive trajectory builder.
# Tests M1 filters, assignment reuse, bin summaries, support and provenance.
# No matching, GLME fitting, bootstrap, or report rendering is performed.
# Usage: Rscript Code/tests/4.0.2_test_descriptive_contracts.R
# =============================================================================

source("Code/2_data_analysis/4.0.2_GLME_descriptive_trajectories.R")

expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), fixed = TRUE))
}

fixture <- data.frame(
  Cohort = "Fixture Cohort",
  Group = c("Cancer Case", "Control", "Control", "Cancer Case", "Control", "Cancer Case"),
  id = c("case", "control", "control", "drop_post", "drop_covariate", "drop_missing_fi"),
  match_set = c("a", "a", "b", "a", "a", "a"),
  role = c("Case", "Control", "Control", "Case", "Control", "Case"),
  cycle = c("2000", "2000", "2000", "2000", "2000", "2000"),
  Age_Centered = c(-1, -1, -1, 1, 1, 1),
  fi_score_nocancer = c(.20, .10, .10, .30, .40, NA_real_),
  post_own_cancer = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  index_age_z = c(0, 0, 0, 0, NA, 0),
  base_race = "White", base_marital = "Married", base_living = "With others",
  stringsAsFactors = FALSE
)
prepared <- .prepare_descriptive_m1_data(fixture)
stopifnot(nrow(prepared) == 3L, sum(prepared$id == "control") == 2L,
          n_distinct(prepared$trajectory_id) == 3L)
provenance <- list(input_md5 = "fixture")
summary <- .summarize_descriptive_cohort(prepared, "fixture", provenance,
  min_case_bin = 1L, min_ctrl_bin = 1L)
means <- summary$observed
means <- subset(means, rel_time_bin == "-2 to 0")
difference <- subset(summary$difference, rel_time_bin == "-2 to 0")
support <- subset(summary$support, rel_time_bin == "-2 to 0")
stopifnot(means$n_observations[means$Group == "Control"] == 2L,
          means$n_assignments[means$Group == "Control"] == 2L,
          means$n_participants[means$Group == "Control"] == 1L,
          means$mean_fi[means$Group == "Control"] == .10)
stopifnot(difference$observed_difference == .10,
          support$n_case_participants == 1L,
          support$n_control_participants == 1L,
          support$display_supported)
stopifnot(identical(.descriptive_bin_midpoint(c("-2 to 0", "0 to +2")), c(-1, 1)),
          is.na(.descriptive_bin_midpoint("<= -20 years")))
bad <- fixture; bad$post_own_cancer <- NA
expect_error(.prepare_descriptive_m1_data(bad), "missing identifiers")
dup <- rbind(fixture[1:3, ], fixture[1, ])
expect_error(.prepare_descriptive_m1_data(dup), "Duplicated trajectory_id x cycle")

# A changed matched input must fail before descriptive data are written.
project <- tempfile("descriptive_provenance_")
dir.create(file.path(project, "Data"), recursive = TRUE)
diagnostics <- file.path(project, "Results", "cancer", "data", "matching_diagnostics")
dir.create(diagnostics, recursive = TRUE)
on.exit(unlink(project, recursive = TRUE), add = TRUE)
fixture_path <- file.path(project, "Data", "riskset_matched_analysis_long.rds")
saveRDS(fixture[1:3, ], fixture_path)
write.csv(data.frame(gate_pass = TRUE),
          file.path(diagnostics, "riskset_matched_analysis_long_gate_g4.csv"), row.names = FALSE)
saveRDS(list(eligibility_version = "no_fi_requirement_v1", fi_requirement = "none",
    cohort_entry_rule = "first_participated_analytic_cycle_return", control_entry_on_or_before_index = TRUE, output_md5 = "intentionally-wrong"),
        file.path(diagnostics, "riskset_matched_analysis_long_run_metadata.rds"))
expect_error(run_descriptive_trajectories(project), "hash does not match")

stopifnot(!any(grepl("lmer|bootstrap|vcov", deparse(run_descriptive_trajectories), ignore.case = TRUE)),
          !any(grepl("lmer|bootstrap|vcov", deparse(render_descriptive_trajectories), ignore.case = TRUE)))
report_source <- readLines("Code/2_data_analysis/4.0.2_GLME_descriptive_trajectories.Rmd", warn = FALSE)
stopifnot(!any(grepl("lmer\\(|bootstrap|vcovCR", report_source, ignore.case = TRUE)))
cat("PASS: descriptive contracts; no matching, GLME, bootstrap, or report rendering executed.\n")
