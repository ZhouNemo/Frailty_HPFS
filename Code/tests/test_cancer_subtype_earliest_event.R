# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  test_cancer_subtype_earliest_event.R
# Author:  Nemo Zhou
# Date started: 2026-07-20
# Date last updated: 2026-07-20
#
# Purpose:
#   Runs deterministic, synthetic regression tests for the earliest-date cancer
#   event reconstruction in 7.4_cancer_subtypes.R. The test uses temporary RDS
#   input/output paths and does not modify the project dataset or result files.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
script_path <- file.path(project_dir, "Code", "1_data_cleaning", "7.4_cancer_subtypes.R")

make_person <- function(id, raw_date,
                        first_icd, first_yy, first_mm,
                        second_icd = NA_character_, second_yy = NA,
                        second_mm = NA) {
  data.frame(
    id = rep(id, 2),
    cycle = c("86", "88"),
    cancer_dateca = rep(raw_date, 2),
    cancer_icdx = rep(first_icd, 2),
    cancer_newag = rep("0", 2),
    cancer_prosnoa1 = rep("0", 2),
    cancer_status = rep("0", 2),
    cancer_firsticd = rep(first_icd, 2),
    cancer_firstyy = rep(first_yy, 2),
    cancer_firstmm = rep(first_mm, 2),
    cancer_secondicd = rep(second_icd, 2),
    cancer_secondyy = rep(second_yy, 2),
    cancer_secondmm = rep(second_mm, 2),
    cancer_thirdicd = NA_character_, cancer_thirdyy = NA, cancer_thirdmm = NA,
    cancer_fourthicd = NA_character_, cancer_fourthyy = NA, cancer_fourthmm = NA,
    stringsAsFactors = FALSE
  )
}

synthetic <- rbind(
  # Later high-burden cancer must not reclassify the earlier prostate cancer.
  make_person("low_then_high", 1081, "185", 90, 1, "162", 90, 6),
  # A later ICD slot can have an earlier diagnosis date than firsticd.
  make_person("later_slot_earlier", 1086, "162", 90, 6, "185", 90, 1),
  # Same-date high and low cancers: high burden wins.
  make_person("same_date_tie", 1086, "185", 90, 6, "162", 90, 6),
  # Later smoking cancer must not create a smoking-related case.
  make_person("unrelated_then_smoking", 1081, "172", 90, 1, "162", 90, 6),
  # Missing first-slot date uses the endpoint-level date fallback.
  make_person("missing_first_date", 1081, "185", NA, NA, "162", 90, 6),
  # Valid but unmapped ICDs are low/moderate by default.
  make_person("unknown_icd", 1081, "199", 90, 1),
  # Duplicate same-date ICD records collapse to one distinct index ICD.
  make_person("duplicate_event", 1081, "185", 90, 1, "185", 90, 1)
)

input_path <- tempfile(fileext = ".rds")
output_path <- tempfile(fileext = ".rds")
saveRDS(synthetic, input_path)

old_input <- Sys.getenv("CANCER_SUBTYPE_INPUT_PATH", unset = NA_character_)
old_output <- Sys.getenv("CANCER_SUBTYPE_OUTPUT_PATH", unset = NA_character_)
Sys.setenv(
  CANCER_SUBTYPE_INPUT_PATH = input_path,
  CANCER_SUBTYPE_OUTPUT_PATH = output_path
)
run_output <- system2("Rscript", shQuote(script_path), stdout = TRUE, stderr = TRUE)
if (is.na(old_input)) Sys.unsetenv("CANCER_SUBTYPE_INPUT_PATH") else Sys.setenv(CANCER_SUBTYPE_INPUT_PATH = old_input)
if (is.na(old_output)) Sys.unsetenv("CANCER_SUBTYPE_OUTPUT_PATH") else Sys.setenv(CANCER_SUBTYPE_OUTPUT_PATH = old_output)
status <- attr(run_output, "status")
if (!is.null(status) && status != 0) {
  stop(paste(run_output, collapse = "\n"))
}

result <- readRDS(output_path)
person <- result[!duplicated(result$id), ]
row_for <- function(id) person[person$id == id, , drop = FALSE]

stopifnot(row_for("low_then_high")$cancer_index_dateca == 1081)
stopifnot(row_for("low_then_high")$is_high_burden == 0L)
stopifnot(row_for("low_then_high")$is_smoking_cancer == 0L)

stopifnot(row_for("later_slot_earlier")$cancer_index_dateca == 1081)
stopifnot(row_for("later_slot_earlier")$cancer_index_icds == "185")
stopifnot(row_for("later_slot_earlier")$is_high_burden == 0L)

stopifnot(row_for("same_date_tie")$cancer_index_dateca == 1086)
stopifnot(row_for("same_date_tie")$cancer_index_n_distinct_icds == 2L)
stopifnot(row_for("same_date_tie")$is_high_burden == 1L)

stopifnot(row_for("unrelated_then_smoking")$is_smoking_cancer == 0L)
stopifnot(row_for("unrelated_then_smoking")$is_obesity_cancer == 0L)

stopifnot(row_for("missing_first_date")$cancer_index_dateca == 1081)
stopifnot(row_for("missing_first_date")$is_high_burden == 0L)

stopifnot(row_for("unknown_icd")$is_high_burden == 0L)
stopifnot(row_for("duplicate_event")$cancer_index_n_distinct_icds == 1L)
stopifnot(row_for("duplicate_event")$cancer_index_icds == "185")

constant_check <- result |>
  dplyr::group_by(id) |>
  dplyr::summarize(
    n_dates = dplyr::n_distinct(cancer_index_dateca, na.rm = TRUE),
    n_flags = dplyr::n_distinct(
      paste(is_high_burden, is_smoking_cancer, is_obesity_cancer, sep = "|")
    ),
    .groups = "drop"
  )
stopifnot(all(constant_check$n_dates <= 1L))
stopifnot(all(constant_check$n_flags <= 1L))

unlink(c(input_path, output_path), force = TRUE)
message("Earliest cancer subtype synthetic tests passed.")
