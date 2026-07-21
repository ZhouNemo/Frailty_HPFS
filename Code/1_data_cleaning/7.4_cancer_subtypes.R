# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health
#          Professionals Follow-up Study
# Script: 7.4_cancer_subtypes.R
# Author: Nemo Zhou
# Date started: 2026-06-29
# Date last updated: 2026-07-20
# Purpose:
#   Reconstructs each participant's earliest cancer diagnosis from the HPFS
#   endpoint ICD records and classifies cancer subgroups using only cancers
#   sharing that earliest diagnosis date. Later cancers cannot reclassify a
#   participant. If simultaneous cancers include a high-burden ICD, the index
#   cancer is high burden.
#
#   The raw endpoint date (cancer_dateca) is preserved. The canonical analytic
#   date written by this script is cancer_index_dateca, measured in months since
#   1900. The event metadata fields provide an audit trail for date selection,
#   simultaneous ICDs, duplicate collapse, and undated later ICD records.
#
#   Smoking-related, obesity-related, high-burden, sensitivity, and prostate
#   flags are derived from the earliest dated cancer event only. Valid ICDs that
#   do not map to an explicitly high-burden site are low/moderate burden by
#   default.
#
# Input : Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
# Output: Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#
# Run order: after 8_compute FI without cancer.R and before active matching,
#            IPCW, or trajectory analyses.
# ===============================================================================

library(dplyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
default_io_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
io_path <- Sys.getenv("CANCER_SUBTYPE_INPUT_PATH", unset = default_io_path)
output_path <- Sys.getenv("CANCER_SUBTYPE_OUTPUT_PATH", unset = io_path)

fi <- readRDS(io_path)

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

to01 <- function(x) {
  x <- to_num(x)
  as.integer(!is.na(x) & x == 1)
}

# HPFS yy values are years since 1900 (for example, 90 = 1990 and
# 108 = 2008); mm is the diagnosis month. The result is the same month-scale
# used by dateca. Invalid source values remain NA and are reported below.
diagnosis_month <- function(yy, mm) {
  yy <- to_num(yy)
  mm <- to_num(mm)
  valid <- is.finite(yy) & yy >= 0 & yy <= 150 &
    is.finite(mm) & mm >= 1 & mm <= 12
  out <- rep(NA_real_, length(yy))
  out[valid] <- 12 * yy[valid] + mm[valid]
  out
}

event_site_flags <- function(icd) {
  icd <- to_num(icd)
  list(
    smoking = (!is.na(icd) &
      ((icd >= 140 & icd <= 149) | icd %in% c(150, 151, 153, 154, 155, 157, 162, 188, 189))),
    smoking_core = (!is.na(icd) &
      ((icd >= 140 & icd <= 149) | icd %in% c(150, 162, 188, 189))),
    obesity = (!is.na(icd) & icd %in% c(150, 151, 153, 154, 155, 157, 189, 203)),
    obesity_sens = (!is.na(icd) & icd %in% c(150, 151, 153, 154, 155, 157, 189, 191, 192, 203)),
    high_burden = (!is.na(icd) & icd %in% c(150, 151, 155, 157, 162, 191, 192)),
    high_burden_sens = (!is.na(icd) &
      icd %in% c(150, 151, 155, 157, 162, 191, 192, 204, 205, 206, 207, 208)),
    prostate = (!is.na(icd) & icd == 185)
  )
}

# -----------------------------------------------------------------------------
# 1. Build and de-duplicate the endpoint cancer event table
# -----------------------------------------------------------------------------
slot_names <- c("first", "second", "third", "fourth")
person_source <- fi %>%
  distinct(id, .keep_all = TRUE) %>%
  select(
    id,
    cancer_dateca,
    cancer_icdx,
    cancer_newag,
    cancer_prosnoa1,
    cancer_status,
    all_of(unlist(lapply(slot_names, function(s) {
      paste0("cancer_", s, c("icd", "mm", "yy"))
    })))
  )

slot_events <- bind_rows(lapply(slot_names, function(s) {
  icd_name <- paste0("cancer_", s, "icd")
  mm_name <- paste0("cancer_", s, "mm")
  yy_name <- paste0("cancer_", s, "yy")
  icd <- if (s == "first") {
    coalesce(to_num(person_source[[icd_name]]), to_num(person_source$cancer_icdx))
  } else {
    to_num(person_source[[icd_name]])
  }
  date <- diagnosis_month(person_source[[yy_name]], person_source[[mm_name]])
  # The endpoint-level date is the documented fallback for the first ICD when
  # its month/year fields are missing. Later undated ICDs remain undated.
  if (s == "first") {
    date <- ifelse(is.na(date), to_num(person_source$cancer_dateca), date)
  }
  data.frame(
    id = as.character(person_source$id),
    record_slot = s,
    icd = icd,
    event_dateca = date,
    stringsAsFactors = FALSE
  )
}))

raw_endpoint_events <- person_source %>%
  transmute(
    id = as.character(id),
    record_slot = "endpoint",
    icd = NA_real_,
    event_dateca = to_num(cancer_dateca)
  ) %>%
  filter(!is.na(event_dateca))

event_records <- bind_rows(slot_events, raw_endpoint_events)
raw_event_count <- nrow(event_records)
undated_icd_count <- sum(!is.na(event_records$icd) & is.na(event_records$event_dateca))

# Repeated endpoint records with the same participant, ICD, and diagnosis
# month are one event for earliest-date and tie classification.
event_records <- event_records %>%
  distinct(id, icd, event_dateca, .keep_all = TRUE)
duplicate_event_count <- raw_event_count - nrow(event_records)

event_flags <- event_site_flags(event_records$icd)
event_records <- event_records %>%
  mutate(
    is_smoking_event = as.integer(event_flags$smoking),
    is_smoking_core_event = as.integer(event_flags$smoking_core),
    is_obesity_event = as.integer(event_flags$obesity),
    is_obesity_sens_event = as.integer(event_flags$obesity_sens),
    is_high_burden_event = as.integer(event_flags$high_burden),
    is_high_burden_sens_event = as.integer(event_flags$high_burden_sens),
    is_prostate_event = as.integer(event_flags$prostate)
  )

# -----------------------------------------------------------------------------
# 2. Select the earliest dated cancer and classify only tied index events
# -----------------------------------------------------------------------------
earliest_events <- event_records %>%
  filter(!is.na(event_dateca)) %>%
  group_by(id) %>%
  filter(event_dateca == min(event_dateca)) %>%
  ungroup()

earliest_summary <- earliest_events %>%
  group_by(id) %>%
  summarize(
    cancer_index_dateca = first(event_dateca),
    cancer_index_icds = {
      codes <- sort(unique(icd[!is.na(icd)]))
      if (length(codes)) paste(codes, collapse = ";") else NA_character_
    },
    cancer_index_n_events = n(),
    cancer_index_n_distinct_icds = n_distinct(icd[!is.na(icd)]),
    cancer_index_tie_high_burden = as.integer(any(is_high_burden_event == 1L)),
    cancer_index_tie_high_burden_sens = as.integer(any(is_high_burden_sens_event == 1L)),
    is_smoking_cancer = as.integer(any(is_smoking_event == 1L)),
    is_smoking_cancer_core = as.integer(any(is_smoking_core_event == 1L)),
    is_obesity_cancer = as.integer(any(is_obesity_event == 1L)),
    is_obesity_cancer_sens = as.integer(any(is_obesity_sens_event == 1L)),
    is_high_burden = as.integer(any(is_high_burden_event == 1L)),
    is_high_burden_sens = as.integer(any(is_high_burden_sens_event == 1L)),
    index_has_prostate = as.integer(any(is_prostate_event == 1L)),
    .groups = "drop"
  )

# Prostate stage flags are participant-level endpoint metadata. They are gated
# by an ICD-185 cancer at the earliest date so later prostate records cannot
# create an index prostate subgroup.
earliest_summary <- earliest_summary %>%
  left_join(
    person_source %>%
      transmute(
        id,
        endpoint_prostate_agg = to01(cancer_newag),
        endpoint_prostate_indolent = as.integer(
          (to01(cancer_prosnoa1) == 1L | to01(cancer_status) == 1L) &
            to01(cancer_newag) != 1L
        )
      ),
    by = "id"
  ) %>%
  mutate(
    is_prostate_agg = as.integer(index_has_prostate == 1L & endpoint_prostate_agg == 1L),
    is_prostate_indolent = as.integer(
      index_has_prostate == 1L & endpoint_prostate_indolent == 1L
    )
  ) %>%
  select(-endpoint_prostate_agg, -endpoint_prostate_indolent, -index_has_prostate)

classification_cols <- c(
  "cancer_index_dateca", "cancer_index_icds", "cancer_index_n_events",
  "cancer_index_n_distinct_icds", "cancer_index_tie_high_burden",
  "cancer_index_tie_high_burden_sens", "is_smoking_cancer",
  "is_smoking_cancer_core", "is_obesity_cancer", "is_obesity_cancer_sens",
  "is_high_burden", "is_high_burden_sens", "is_prostate_agg",
  "is_prostate_indolent"
)

fi <- fi %>%
  select(-any_of(classification_cols)) %>%
  left_join(earliest_summary, by = "id") %>%
  mutate(
    across(
      c(
        is_smoking_cancer, is_smoking_cancer_core, is_obesity_cancer,
        is_obesity_cancer_sens, is_high_burden, is_high_burden_sens,
        is_prostate_agg, is_prostate_indolent
      ),
      ~ coalesce(as.integer(.x), 0L)
    )
  )

# -----------------------------------------------------------------------------
# 3. Quality checks and save
# -----------------------------------------------------------------------------
new_flags <- c(
  "is_smoking_cancer", "is_smoking_cancer_core", "is_obesity_cancer",
  "is_obesity_cancer_sens", "is_high_burden", "is_high_burden_sens",
  "is_prostate_agg", "is_prostate_indolent"
)

person_output <- fi %>% distinct(id, .keep_all = TRUE)
stopifnot(all(c(classification_cols, new_flags) %in% names(fi)))
stopifnot(all(!is.na(person_output$cancer_index_dateca) |
                is.na(person_output$cancer_index_icds)))

constant_check <- fi %>%
  group_by(id) %>%
  summarize(
    n_index_dates = n_distinct(cancer_index_dateca, na.rm = TRUE),
    n_flag_patterns = n_distinct(do.call(
      paste,
      c(across(all_of(new_flags)), sep = "|")
    )),
    .groups = "drop"
  )
stopifnot(all(constant_check$n_index_dates <= 1L))
stopifnot(all(constant_check$n_flag_patterns <= 1L))

cat("Endpoint event records before duplicate collapse:", raw_event_count, "\n")
cat("Duplicate event records collapsed:", duplicate_event_count, "\n")
cat("Undated non-first ICD records retained for diagnostics:", undated_icd_count, "\n")
cat("Participants with an earliest dated cancer:", nrow(earliest_summary), "\n\n")
cat("Earliest-date person-level classification counts:\n")
print(colSums(person_output[, new_flags, drop = FALSE], na.rm = TRUE))

output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
saveRDS(fi, output_path)
cat("\nEarliest-date cancer classifications written to:", output_path, "\n")
