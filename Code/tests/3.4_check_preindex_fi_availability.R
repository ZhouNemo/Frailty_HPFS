# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-29
# Date last updated: 2026-09-29
# Purpose: Audit FI availability in the final rebuilt matched cohorts.
# Reads all five complete assignment ledgers and observed outcome files; checks
# Gate G4, no_fi_requirement_v1, file hashes, assignment integrity and FI counts.
# Uses actual month-coded return/index dates, not cycle labels or rounded ages.
# Reports strict pre-index FI (<), index-month FI (=), and post-index FI (>).
# Primary summaries exclude control rows at/after their later own cancer;
# a separate saved-outcome view retains them for reconciliation. No matching,
# model fitting, synthetic data or changes to source datasets occur here.
# Outputs: timestamped aggregate CSV/RDS audit under Results/cancer/data/
# 3.4_preindex_fi_availability/. No participant-level output or PNG is saved.
# Run from project root: Rscript Code/tests/3.4_check_preindex_fi_availability.R
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
source("Code/2_data_analysis/2.0_matching_provenance.R")

check_fi <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

summarize_fi_assignments <- function(a, rows, view) {
  counts <- rows %>% group_by(trajectory_id) %>% summarize(
    n_before = n_distinct(worked_rtmnyr[worked_rtmnyr < index_date]),
    n_index_month = n_distinct(worked_rtmnyr[worked_rtmnyr == index_date]),
    n_after = n_distinct(worked_rtmnyr[worked_rtmnyr > index_date]),
    n_all = n_distinct(worked_rtmnyr), .groups = "drop")
  out <- a %>% left_join(counts, by = "trajectory_id") %>%
    mutate(across(c(n_before, n_index_month, n_after, n_all), ~coalesce(.x, 0L)),
      view = view,
      no_strict_preindex_fi = n_before == 0,
      no_fi_on_or_before_index = n_before + n_index_month == 0,
      category = case_when(
        n_before > 0 ~ "FI strictly before index",
        n_index_month > 0 ~ "No earlier FI; FI in index month (with or without later FI)",
        n_after > 0 ~ "Post-index FI only",
        TRUE ~ "No usable FI anywhere"))
  check_fi(all(out$n_all == out$n_before + out$n_index_month + out$n_after),
           "FI date partitions do not reconcile.")
  out
}

run_preindex_fi_audit <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  stems <- c("analysis", "smoking", "obesity", "survival", "overall")
  summaries <- categories <- people <- gates <- provenance <- list()
  for (stem in stems) {
    path <- file.path(project_dir, "Data", paste0("riskset_matched_", stem, "_long.rds"))
    pv <- validate_matching_provenance(path)
    a <- readRDS(pv$assignment_path)
    d <- readRDS(path)
    assert_assignment_rows(a, paste(stem, "ledger"))
    assert_assignment_rows(d, paste(stem, "outcomes"))
    keys <- c("trajectory_id", "id", "Cohort", "role", "Group", "match_set")
    a <- a %>% mutate(across(all_of(keys), as.character))
    d <- d %>% mutate(across(all_of(keys), as.character), cycle = as.character(cycle))
    check_fi(!anyNA(a[keys]) && !anyDuplicated(a$trajectory_id), "Invalid ledger keys.")
    check_fi(!anyDuplicated(d[c("trajectory_id", "cycle")]), "Duplicate outcome keys.")
    check_fi(!anyNA(d[keys]) && all(is.finite(d$worked_rtmnyr)) &&
      all(is.finite(d$fi_score_nocancer)), "Missing outcome identifiers, FI or timing.")
    link <- c(keys, "index_date")
    check_fi(nrow(anti_join(distinct(d, across(all_of(link))), a, by = link)) == 0,
             "Observed rows disagree with the assignment ledger.")
    sets <- a %>% group_by(Cohort, match_set) %>% summarize(
      cases = sum(role == "Case"), controls = sum(role == "Control"), .groups = "drop")
    check_fi(all(sets$cases == 1 & sets$controls >= 1 & sets$controls <= 5), "Invalid ledger sets.")
    controls <- a %>% filter(role == "Control")
    check_fi(all(abs(a$age_gap) <= 2) &&
      all(is.na(controls$dtdth) | controls$dtdth > controls$index_date) &&
      all(is.na(controls$cancer_dateca) | controls$cancer_dateca > controls$index_date) &&
      all(controls$id != controls$donor_case_id), "Control eligibility integrity failed.")
    expected <- d$role == "Control" & !is.na(d$cancer_dateca) &
      d$cancer_dateca > d$index_date & d$worked_rtmnyr >= d$cancer_dateca
    check_fi(!anyNA(d$post_own_cancer) && all(d$post_own_cancer == expected),
             "Own-cancer censoring flags disagree with dates.")
    saved <- summarize_fi_assignments(a, d, "Saved FI rows before own-cancer exclusion")
    primary <- summarize_fi_assignments(a, filter(d, !post_own_cancer), "Primary uncensored FI")
    check_fi(all(saved$n_all == a$n_analytic_fi_visits) &&
      all(saved$n_before + saved$n_index_month == a$n_preindex_fi) &&
      all(saved$n_after == a$n_postindex_fi) &&
      all(primary$n_all == a$n_uncensored_fi_visits), "Ledger FI counts do not reconcile with observed rows.")
    z <- bind_rows(saved, primary)
    summaries[[stem]] <- z %>% group_by(Cohort, Group, view) %>% summarize(
      n_assignments = n(), n_unique_participants = n_distinct(id),
      n_no_strict_preindex_fi = sum(no_strict_preindex_fi),
      pct_no_strict_preindex_fi = 100 * mean(no_strict_preindex_fi),
      n_no_fi_on_or_before_index = sum(no_fi_on_or_before_index),
      pct_no_fi_on_or_before_index = 100 * mean(no_fi_on_or_before_index),
      n_index_month_without_earlier_fi = sum(n_before == 0 & n_index_month > 0),
      n_postindex_only = sum(n_before == 0 & n_index_month == 0 & n_after > 0),
      n_no_fi_anywhere = sum(n_all == 0),
      pct_no_fi_anywhere = 100 * mean(n_all == 0), .groups = "drop") %>% mutate(dataset = stem, .before = 1)
    check_fi(all(with(summaries[[stem]], n_no_strict_preindex_fi ==
      n_index_month_without_earlier_fi + n_postindex_only + n_no_fi_anywhere)), "Availability categories do not reconcile.")
    categories[[stem]] <- z %>% count(Cohort, Group, view, category, name = "n_assignments") %>%
      group_by(Cohort, Group, view) %>% mutate(denominator = sum(n_assignments),
        percent = 100 * n_assignments / denominator) %>% ungroup() %>% mutate(dataset = stem, .before = 1)
    people[[stem]] <- z %>% group_by(Cohort, Group, view, id) %>% summarize(
      any_assignment_without_prior_fi = any(no_strict_preindex_fi),
      all_assignments_without_prior_fi = all(no_strict_preindex_fi),
      any_assignment_without_any_fi = any(n_all == 0),
      all_assignments_without_any_fi = all(n_all == 0), .groups = "drop") %>%
      group_by(Cohort, Group, view) %>% summarize(n_unique_participants = n(),
        across(starts_with(c("any_assignment", "all_assignments")),
          list(n = ~sum(.x), pct = ~100 * mean(.x))), .groups = "drop") %>% mutate(dataset = stem, .before = 1)
    gates[[stem]] <- pv$gate %>% mutate(dataset = stem, .before = 1)
    provenance[[stem]] <- list(outcome_path = path, assignment_path = pv$assignment_path,
      output_md5 = pv$input_md5, assignment_md5 = pv$assignment_md5,
      eligibility_version = pv$run_metadata$eligibility_version, matching_created_at = pv$run_metadata$created_at,
      n_assignments = nrow(a), n_outcome_rows = nrow(d))
    check_fi(identical(unname(tools::md5sum(path)), pv$input_md5) &&
      identical(unname(tools::md5sum(pv$assignment_path)), pv$assignment_md5), "Source files changed during audit.")
    message(stem, ": Gate G4, provenance, assignment integrity and outcome counts passed.")
  }
  out <- file.path(project_dir, "Results/cancer/data/3.4_preindex_fi_availability",
                   paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()))
  check_fi(!dir.exists(out), "Audit output directory already exists.")
  dir.create(out, recursive = TRUE)
  write.csv(bind_rows(gates), file.path(out, "gate_g4.csv"), row.names = FALSE)
  write.csv(bind_rows(summaries), file.path(out, "assignment_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(categories), file.path(out, "availability_categories.csv"), row.names = FALSE)
  write.csv(bind_rows(people), file.path(out, "unique_participant_summary.csv"), row.names = FALSE)
  saveRDS(list(created_at = Sys.time(), inputs = provenance,
    script_md5 = unname(tools::md5sum(file.path(project_dir, "Code/tests/3.4_check_preindex_fi_availability.R"))),
    definition = "Strict pre-index: worked_rtmnyr < index_date; index month: equality. Counts are distinct dates.",
    interpretation = "Primary denominators are selected assignments, not model-complete rows. Unique participants are within cohort and role; reuse and cross-role overlap prohibit summing across groups.",
    session = sessionInfo()), file.path(out, "run_metadata.rds"))
  print(bind_rows(summaries) %>% filter(dataset == "overall", view == "Primary uncensored FI"), width = Inf)
  message("Saved audit: ", out)
  invisible(out)
}
if (sys.nframe() == 0L) run_preindex_fi_audit()
