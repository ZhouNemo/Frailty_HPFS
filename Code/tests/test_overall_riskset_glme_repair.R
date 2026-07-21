# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  test_overall_riskset_glme_repair.R
# Author:  Nemo Zhou
# Date started:      2026-07-18
# Date last updated: 2026-07-20 (canonical earliest-cancer index date support)
#
# Purpose:
#   Runs deterministic synthetic regression tests for the repaired overall-
#   cancer risk-set matching contract. The tests verify pre-index-only support,
#   symmetric cycle eligibility, one-visit inclusion, unrestricted relative-
#   time retention, case-reference age scaling, exact-cycle behavior, and
#   post-expansion matched-set/duplicate-key safeguards. No project dataset or
#   result is created or modified; all synthetic inputs live in tempdir().
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
source(file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R"))

expect_error <- function(expr, pattern) {
  observed <- tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
  if (is.na(observed) || !grepl(pattern, observed, fixed = TRUE)) {
    stop("Expected error containing '", pattern, "'; observed: ", observed)
  }
  invisible(observed)
}

make_row <- function(id, dob, cancer_date, cycle, return_date, fi, participated = 1) {
  data.frame(
    id = id,
    cycle = cycle,
    cancer_dateca = cancer_date,
    cancer_index_dateca = cancer_date,
    dtdth = NA_real_,
    dbmy09 = dob,
    worked_rtmnyr = return_date,
    date_of_return = return_date,
    fi_score = fi,
    fi_score_nocancer = fi,
    fi_score_nocancer_nocarry = fi,
    participated = participated,
    race = "White",
    marital_status = "Married",
    living_arr = "With others",
    pckgr = "Never",
    base_living = "With others",
    base_calor = 2200,
    base_sat = 22,
    base_diet_chol = 280,
    base_nblnk = 0,
    base_alco = 8,
    stringsAsFactors = FALSE
  )
}

target_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")
cycle_date <- c(`86` = 1032, `88` = 1056, `92` = 1104, `04` = 1248,
                `08` = 1296, `12` = 1344, `20` = 1440)

# Two cases with different diagnosis ages ensure finite case-reference scaling.
# Each person has one nonmissing pre-index FI, an index-cycle participation row,
# and a later FI. Case/control A also has a +28-year row to test no ±20 filter.
spec <- list(
  case_a = list(dob = 240, cancer = 1104, rows = c("86", "88", "92", "20")),
  ctrl_a = list(dob = 246, cancer = NA_real_, rows = c("86", "88", "92", "20")),
  case_b = list(dob = 360, cancer = 1296, rows = c("86", "04", "08", "12")),
  ctrl_b = list(dob = 366, cancer = NA_real_, rows = c("86", "04", "08", "12"))
)

synthetic <- do.call(rbind, lapply(names(spec), function(pid) {
  x <- spec[[pid]]
  do.call(rbind, lapply(x$rows, function(cyc) {
    pre_fi_cycle <- if (pid %in% c("case_a", "ctrl_a")) "88" else "04"
    fi <- if (cyc == pre_fi_cycle) 0.10 else if (cyc %in% c("12", "20")) 0.18 else NA_real_
    make_row(pid, x$dob, x$cancer, cyc, cycle_date[[cyc]], fi)
  }))
}))
rownames(synthetic) <- NULL

tmp_input <- tempfile(fileext = ".rds")
saveRDS(synthetic, tmp_input)
classify_all <- function(person_level) rep("All Cancer Cohort", nrow(person_level))

primary <- build_riskset_matched_long(
  input_path = tmp_input,
  classification_vars = character(0),
  classify_fn = classify_all,
  cohort_levels = "All Cancer Cohort",
  target_cycles = target_cycles,
  match_ratio = 1,
  age_caliper = 2,
  cycle_caliper = 1,
  min_visits = 1,
  seed = 20260703,
  max_unmatched_fraction = 1,
  max_absolute_smd = Inf
)

# One pre-index FI is sufficient even when the only second FI is post-index.
stopifnot(nrow(primary$inactive_cases) == 0)
stopifnot(all(primary$caliper_integrity$all_preindex_fi_ok))
stopifnot(all(primary$matched_long$n_preindex_fi == 1L))

# Post-index FI cannot rescue eligibility.
future_only <- data.frame(
  id = "future_only", cycle = c("08", "12"), worked_rtmnyr = c(1296, 1344),
  fi_score_nocancer = c(NA_real_, 0.2), stringsAsFactors = FALSE
)
future_support <- preindex_support(future_only, 1296, "08", 1, min_visits = 1)
stopifnot(future_support$has_active_return, future_support$cycle_ok)
stopifnot(future_support$n_preindex_fi == 0L, !future_support$eligible)

# The same cycle boundary applies to cases and controls; exact-cycle S7 rejects
# an adjacent return that the primary ≤4-year caliper accepts.
adjacent <- data.frame(
  id = "adjacent", cycle = "04", worked_rtmnyr = 1248,
  fi_score_nocancer = 0.1, stringsAsFactors = FALSE
)
stopifnot(preindex_support(adjacent, 1296, "08", 1, 1)$eligible)
stopifnot(!preindex_support(adjacent, 1296, "08", 0, 1)$eligible)

# All eligible rows are retained, including the +28-year observation.
stopifnot(any(primary$matched_long$Age_Centered > 20))
stopifnot(any(!primary$matched_long$in_win_20))

# Scaling uses exactly one retained case per set, not expanded rows or controls.
case_assignments <- primary$matched_long |>
  dplyr::filter(role == "Case") |>
  dplyr::distinct(Cohort, match_set, index_age)
stopifnot(primary$scaling_metadata$n_reference_cases == nrow(case_assignments))
stopifnot(isTRUE(all.equal(
  primary$scaling_metadata$index_age_mean,
  mean(case_assignments$index_age), tolerance = 1e-12
)))
stopifnot(isTRUE(all.equal(
  primary$scaling_metadata$index_age_sd,
  stats::sd(case_assignments$index_age), tolerance = 1e-12
)))

# Every returned set has one case and at least one control; exact keys are unique.
stopifnot(all(primary$set_integrity$valid_set))
stopifnot(!anyDuplicated(primary$matched_long[c("trajectory_id", "cycle")]))

# Duplicate analytic rows are rejected before an invalid matched object can be saved.
duplicate_input <- rbind(synthetic, synthetic[synthetic$id == "case_a" & synthetic$cycle == "88", ])
tmp_duplicate <- tempfile(fileext = ".rds")
saveRDS(duplicate_input, tmp_duplicate)
expect_error(
  build_riskset_matched_long(
    input_path = tmp_duplicate,
    classification_vars = character(0),
    classify_fn = classify_all,
    cohort_levels = "All Cancer Cohort",
    target_cycles = target_cycles,
    match_ratio = 1,
    age_caliper = 2,
    cycle_caliper = 1,
    min_visits = 1,
    seed = 20260703,
    max_unmatched_fraction = 1,
    max_absolute_smd = Inf
  ),
  "Duplicated trajectory_id x cycle rows"
)

message("All overall risk-set/GLME repair synthetic tests passed.")
