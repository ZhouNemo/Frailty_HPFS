# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  9.0_cancer_free_matching_functions.R
# Author:  Nemo Zhou
# Date started:      2026-07-19
# Date last updated: 2026-07-20 (canonical earliest-cancer index date)
#
# Purpose:
#   Self-contained cancer-free-through-full-HPFS-endpoint risk-set matching
#   utilities for the high/low-burden cancer sensitivity analysis. This file
#   preserves the active 2.0 matching design for age, active follow-up, cycle,
#   and pre-index FI eligibility, but requires selected controls to have no
#   cancer diagnosis on or before the maximum observed cancer ascertainment
#   month in the input dataset.
#
#   The helper is intentionally separate from 2.0 so the existing primary
#   matching datasets and all downstream 3.xx reports remain unchanged.
#
# Units:
#   Date variables are months since 1900. Age is calculated as
#   (date - dbmy09) / 12. Age_Centered is years relative to the assigned index.
# =============================================================================

library(dplyr)

coalesce0 <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  dplyr::coalesce(x, 0)
}

cancer_free_through_endpoint <- function(cancer_dates, endpoint_month) {
  endpoint_month <- suppressWarnings(as.numeric(endpoint_month))
  if (length(endpoint_month) != 1L || !is.finite(endpoint_month)) {
    stop("endpoint_month must be one finite numeric month value.", call. = FALSE)
  }
  cancer_dates <- suppressWarnings(as.numeric(as.character(cancer_dates)))
  is.na(cancer_dates) | cancer_dates > endpoint_month
}

cycle_to_number <- function(x) {
  x_chr <- as.character(x)
  out <- suppressWarnings(as.integer(x_chr))
  # HPFS cycle labels cross the century at "00".  The analytic calendar is
  # 1988, 1992, 1996, 2000, ..., 2020; "00" must therefore map to 2000.
  out <- ifelse(
    !is.na(out) & out >= 0 & out <= 20,
    2000L + out,
    ifelse(!is.na(out) & out >= 80 & out < 100, 1900L + out, out)
  )
  out
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[[1]]
}

optional_cols <- function(df, cols) {
  intersect(cols, names(df))
}

ensure_optional_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df
}

assign_index_cycle <- function(target_cycles, index_date) {
  target_cycles <- as.character(target_cycles)
  target_years <- cycle_to_number(target_cycles)

  if (is.na(index_date) || length(target_cycles) == 0 || anyNA(target_years)) {
    return(list(cycle = NA_character_, rule = "no_cycle"))
  }

  diagnosis_year <- 1900L + floor(as.numeric(index_date) / 12)
  on_or_after <- which(target_years >= diagnosis_year)
  if (length(on_or_after) == 0) {
    return(list(cycle = NA_character_, rule = "no_future_analytic_cycle"))
  }

  idx <- on_or_after[[1]]
  rule <- if (target_years[[idx]] == diagnosis_year) "containing" else "first_after"
  list(cycle = target_cycles[[idx]], rule = rule)
}

latest_observed_cycle_before <- function(obs_dates, obs_cycles, index_date) {
  ok <- !is.na(obs_dates) & obs_dates <= index_date
  if (!any(ok)) {
    return(list(cycle = NA_character_, date = NA_real_))
  }

  pos <- which(ok)
  idx <- pos[which.max(obs_dates[pos])]
  list(cycle = as.character(obs_cycles[[idx]]), date = obs_dates[[idx]])
}

preindex_support <- function(obs, index_date, index_cycle, cycle_caliper,
                             min_visits = 1L) {
  empty <- list(
    active_cycle = NA_character_, active_date = NA_real_, cycle_gap = NA_real_,
    n_preindex_fi = 0L, latest_preindex_fi_date = NA_real_,
    has_active_return = FALSE, cycle_ok = FALSE, fi_ok = FALSE, eligible = FALSE
  )
  if (is.null(obs) || is.na(index_date) || is.na(index_cycle)) return(empty)

  valid_dates <- !is.na(obs$worked_rtmnyr)
  if (!all(valid_dates)) obs <- obs[valid_dates, , drop = FALSE]
  if (nrow(obs) == 0) return(empty)
  if (is.unsorted(obs$worked_rtmnyr)) {
    obs <- obs[order(obs$worked_rtmnyr), , drop = FALSE]
  }
  if (!("cum_distinct_fi" %in% names(obs))) {
    fi_seen <- !is.na(obs$fi_score_nocancer)
    distinct_fi <- logical(nrow(obs))
    fi_pos <- which(fi_seen)
    distinct_fi[fi_pos] <- !duplicated(obs$worked_rtmnyr[fi_pos])
    obs$cum_distinct_fi <- cumsum(distinct_fi)
    latest_marker <- ifelse(fi_seen, obs$worked_rtmnyr, -Inf)
    obs$cum_latest_fi_date <- cummax(latest_marker)
    obs$cum_latest_fi_date[!is.finite(obs$cum_latest_fi_date)] <- NA_real_
  }
  active_pos <- findInterval(index_date, obs$worked_rtmnyr)
  if (active_pos == 0L) return(empty)
  active_cycle <- as.character(obs$cycle[[active_pos]])
  active_date <- as.numeric(obs$worked_rtmnyr[[active_pos]])
  cycle_gap <- abs(cycle_to_number(active_cycle) - cycle_to_number(index_cycle))
  n_preindex_fi <- as.integer(obs$cum_distinct_fi[[active_pos]])
  latest_fi <- as.numeric(obs$cum_latest_fi_date[[active_pos]])
  cycle_ok <- !is.na(cycle_gap) && cycle_gap <= (4L * cycle_caliper)
  fi_ok <- n_preindex_fi >= min_visits

  list(
    active_cycle = active_cycle,
    active_date = active_date,
    cycle_gap = cycle_gap,
    n_preindex_fi = as.integer(n_preindex_fi),
    latest_preindex_fi_date = latest_fi,
    has_active_return = TRUE,
    cycle_ok = cycle_ok,
    fi_ok = fi_ok,
    eligible = cycle_ok && fi_ok
  )
}

prepare_riskset_inputs <- function(input_path,
                                   target_cycles,
                                   min_visits,
                                   classification_vars = character(0)) {
  fi_long <- readRDS(input_path)

  if (!"cancer_index_dateca" %in% names(fi_long)) {
    stop(
      "Missing cancer_index_dateca in input dataset. Run ",
      "Code/1_data_cleaning/7.4_cancer_subtypes.R first.",
      call. = FALSE
    )
  }

  fi_time <- fi_long %>%
    mutate(
      id                = as.character(id),
      cycle             = as.character(cycle),
      cancer_dateca_endpoint = as.numeric(as.character(cancer_dateca)),
      cancer_dateca     = as.numeric(as.character(cancer_index_dateca)),
      dtdth             = as.numeric(as.character(dtdth)),
      dbmy09            = as.numeric(as.character(dbmy09)),
      worked_rtmnyr     = as.numeric(as.character(worked_rtmnyr)),
      date_of_return    = if ("date_of_return" %in% names(.)) date_of_return else NA,
      fi_score_nocancer = as.numeric(as.character(fi_score_nocancer)),
      participated      = as.numeric(as.character(participated)),
      age_at_cancer     = (cancer_dateca - dbmy09) / 12
    )

  cancer_date_values <- fi_time$cancer_dateca[is.finite(fi_time$cancer_dateca)]
  if (length(cancer_date_values) == 0L) {
    stop("No finite cancer_dateca values found; cannot derive the HPFS endpoint.",
         call. = FALSE)
  }
  cancer_free_cutoff <- max(cancer_date_values)

  missing_class_vars <- setdiff(classification_vars, names(fi_time))
  if (length(missing_class_vars) > 0) {
    stop(
      "Missing cancer classification columns in input dataset: ",
      paste(missing_class_vars, collapse = ", "),
      ". Run Code/1_data_cleaning/7.4_cancer_subtypes.R first."
    )
  }

  if (length(classification_vars) > 0) {
    fi_time <- fi_time %>%
      mutate(across(all_of(classification_vars), coalesce0))

    classification_consistency <- fi_time %>%
      group_by(id) %>%
      summarize(across(all_of(classification_vars), n_distinct), .groups = "drop")
    if (any(as.matrix(classification_consistency[, classification_vars, drop = FALSE]) > 1L)) {
      stop(
        "Cancer classification fields are not constant within participant;",
        " refusing to reconstruct an ever-status from later rows.",
        call. = FALSE
      )
    }
  }

  fi_time <- ensure_optional_cols(
    fi_time,
    c("base_living", "base_calor", "base_sat", "base_diet_chol", "base_nblnk", "base_alco")
  )

  fi_trajectory <- fi_time %>%
    filter(
      participated == 1,
      cycle %in% target_cycles,
      !is.na(fi_score_nocancer),
      !is.na(worked_rtmnyr),
      !is.na(dbmy09)
    ) %>%
    mutate(age_at_cycle = (worked_rtmnyr - dbmy09) / 12) %>%
    arrange(id, age_at_cycle)

  observation_panel <- fi_time %>%
    filter(
      participated == 1,
      cycle %in% target_cycles,
      !is.na(worked_rtmnyr),
      !is.na(dbmy09)
    ) %>%
    mutate(
      observed_raw = TRUE,
      age_at_cycle = (worked_rtmnyr - dbmy09) / 12
    ) %>%
    arrange(id, worked_rtmnyr)

  baseline_covs <- fi_time %>%
    filter(cycle == "86") %>%
    distinct(id, .keep_all = TRUE) %>%
    transmute(
      id,
      base_race    = factor(if_else(is.na(race), "Missing", as.character(race))),
      base_marital = factor(if_else(is.na(marital_status), "Missing", as.character(marital_status))),
      base_living  = factor(if_else(is.na(base_living),
                                    if_else(is.na(living_arr), "Missing", as.character(living_arr)),
                                    as.character(base_living))),
      base_pckgr   = factor(if_else(is.na(pckgr), "Missing", as.character(pckgr))),
      base_calor = as.numeric(as.character(base_calor)),
      base_sat = as.numeric(as.character(base_sat)),
      base_diet_chol = as.numeric(as.character(base_diet_chol)),
      base_nblnk = as.numeric(as.character(base_nblnk)),
      base_alco = as.numeric(as.character(base_alco))
    )

  # Participation, rather than the number or span of FI outcomes over the full
  # record, defines the population that can enter an index-time risk set.
  # Index-specific FI support is evaluated later using only rows on/before index.
  person_level <- observation_panel %>%
    group_by(id) %>%
    arrange(age_at_cycle, .by_group = TRUE) %>%
    summarize(
      dbmy09             = first(dbmy09),
      dtdth              = first(dtdth),
      cancer_dateca      = first(cancer_dateca),
      true_age_at_cancer = first(age_at_cancer),
      base_age           = first(age_at_cycle),
      first_return       = first(worked_rtmnyr),
      last_return        = last(worked_rtmnyr),
      n_visits           = n_distinct(age_at_cycle[!is.na(fi_score_nocancer)]),
      min_age            = min(age_at_cycle),
      max_age            = max(age_at_cycle),
      across(all_of(classification_vars), ~ first(coalesce0(.x))),
      .groups = "drop"
    ) %>%
    mutate(
      age_at_death     = (dtdth - dbmy09) / 12,
      death_status     = if_else(!is.na(dtdth), 1L, 0L),
      is_prevalent     = if_else(!is.na(true_age_at_cancer) & true_age_at_cancer <= base_age, 1L, 0L),
      is_case          = if_else(!is.na(true_age_at_cancer) & true_age_at_cancer >  base_age, 1L, 0L)
    ) %>%
    filter(is_prevalent == 0)

  index_info <- observation_panel %>%
    semi_join(person_level, by = "id") %>%
    group_by(id) %>%
    summarize(
      index_cycle_info = list(assign_index_cycle(target_cycles, first(cancer_dateca))),
      .groups = "drop"
    ) %>%
    mutate(
      index_cycle = vapply(index_cycle_info, `[[`, character(1), "cycle"),
      index_cycle_rule = vapply(index_cycle_info, `[[`, character(1), "rule")
    ) %>%
    select(-index_cycle_info)

  person_level <- person_level %>%
    left_join(index_info, by = "id")

  list(
    fi_trajectory = fi_trajectory,
    observation_panel = observation_panel,
    baseline_covs = baseline_covs,
    person_level = person_level,
    cancer_free_cutoff = cancer_free_cutoff
  )
}

match_one_cohort <- function(case_df,
                             pool,
                             observation_by_id,
                             ratio,
                             age_caliper,
                             cycle_caliper,
                             min_visits,
                             cohort_label,
                             seed = 20260703,
                             control_cancer_free_cutoff) {
  if (length(control_cancer_free_cutoff) != 1L ||
      !is.finite(control_cancer_free_cutoff)) {
    stop("control_cancer_free_cutoff must be one finite numeric month value.",
         call. = FALSE)
  }
  match_case <- function(i) {
    idate <- case_df$cancer_dateca[[i]]
    iage  <- case_df$true_age_at_cancer[[i]]
    cid   <- case_df$id[[i]]
    icyc  <- case_df$index_cycle[[i]]

    age_at_index <- (idate - pool$dob) / 12
    not_self <- pool$id != cid
    alive <- is.na(pool$dth) | pool$dth > idate
    cancer_free <- cancer_free_through_endpoint(
      pool$canc,
      control_cancer_free_cutoff
    )
    age_ok <- abs(age_at_index - iage) <= age_caliper
    base_elig <- not_self & alive & cancer_free & age_ok
    base_idx <- which(base_elig)
    base_ids <- pool$id[base_idx]
    base_age_at_index <- age_at_index[base_idx]

    support <- lapply(base_ids, function(pid) {
      obs <- observation_by_id[[as.character(pid)]]
      preindex_support(obs, idate, icyc, cycle_caliper, min_visits)
    })
    active_cycle <- vapply(support, `[[`, character(1), "active_cycle")
    active_date <- vapply(support, `[[`, numeric(1), "active_date")
    cycle_gap <- vapply(support, `[[`, numeric(1), "cycle_gap")
    n_preindex_fi <- vapply(support, `[[`, integer(1), "n_preindex_fi")
    latest_preindex_fi_date <- vapply(support, `[[`, numeric(1), "latest_preindex_fi_date")
    active_followup <- vapply(support, `[[`, logical(1), "has_active_return")
    cycle_ok <- vapply(support, `[[`, logical(1), "cycle_ok")
    fi_ok <- vapply(support, `[[`, logical(1), "fi_ok")
    elig_post_base <- active_followup & cycle_ok & fi_ok
    cand <- data.frame(
      id = base_ids[elig_post_base],
      control_index_age = base_age_at_index[elig_post_base],
      age_gap = base_age_at_index[elig_post_base] - iage,
      control_active_cycle = active_cycle[elig_post_base],
      control_active_date = active_date[elig_post_base],
      cycle_gap = cycle_gap[elig_post_base],
      n_preindex_fi = n_preindex_fi[elig_post_base],
      latest_preindex_fi_date = latest_preindex_fi_date[elig_post_base],
      stringsAsFactors = FALSE
    )

    fail <- data.frame(
      match_set = cid,
      Cohort = cohort_label,
      failed_not_self = sum(!not_self, na.rm = TRUE),
      failed_alive = sum(not_self & !alive, na.rm = TRUE),
      failed_cancer_free = sum(not_self & alive & !cancer_free, na.rm = TRUE),
      failed_cancer_free_through_endpoint =
        sum(not_self & alive & !cancer_free, na.rm = TRUE),
      failed_age_caliper = sum(not_self & alive & cancer_free & !age_ok, na.rm = TRUE),
      failed_active_followup = sum(!active_followup, na.rm = TRUE),
      failed_cycle_caliper = sum(active_followup & !cycle_ok, na.rm = TRUE),
      failed_preindex_fi = sum(active_followup & cycle_ok & !fi_ok, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (nrow(cand) == 0) {
      return(list(assignment = NULL, n_eligible = 0L, fail = fail))
    }

    # Generate tie-breaking draws from a case-specific seed so parallel
    # execution remains reproducible and independent of worker scheduling.
    set.seed(as.integer((as.double(seed) + i) %% 2147483647))
    picked <- cand %>%
      mutate(tie_break = runif(n())) %>%
      arrange(abs(age_gap), tie_break) %>%
      slice_head(n = min(ratio, nrow(cand))) %>%
      select(-tie_break)

    assignment <- data.frame(
      id = picked$id,
      match_set = cid,
      index_date = idate,
      index_cycle = icyc,
      index_cycle_rule = case_df$index_cycle_rule[[i]],
      donor_case_id = cid,
      riskset_size = nrow(cand),
      age_gap = picked$age_gap,
      cycle_gap = picked$cycle_gap,
      control_active_cycle = picked$control_active_cycle,
      control_active_date = picked$control_active_date,
      n_preindex_fi = picked$n_preindex_fi,
      latest_preindex_fi_date = picked$latest_preindex_fi_date,
      role = "Control",
      Cohort = cohort_label,
      stringsAsFactors = FALSE
    )
    list(assignment = assignment, n_eligible = nrow(cand), fail = fail)
  }

  if (nrow(case_df) == 0) {
    case_results <- list()
  } else {
    detected_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(detected_cores)) detected_cores <- 4L
    mc_cores <- min(4L, max(1L, detected_cores - 1L))
    case_results <- parallel::mclapply(
      seq_len(nrow(case_df)),
      match_case,
      mc.preschedule = TRUE,
      mc.cores = mc_cores
    )
  }

  assign_list <- lapply(case_results, `[[`, "assignment")
  n_eligible <- vapply(case_results, `[[`, integer(1), "n_eligible")
  fail_counts <- lapply(case_results, `[[`, "fail")
  controls <- bind_rows(assign_list)
  cases_a <- data.frame(
    id = case_df$id,
    match_set = case_df$id,
    index_date = case_df$cancer_dateca,
    index_cycle = case_df$index_cycle,
    index_cycle_rule = case_df$index_cycle_rule,
    donor_case_id = case_df$id,
    riskset_size = n_eligible,
    age_gap = 0,
    cycle_gap = case_df$case_cycle_gap,
    control_active_cycle = case_df$case_active_cycle,
    control_active_date = case_df$case_active_date,
    n_preindex_fi = case_df$n_preindex_fi,
    latest_preindex_fi_date = case_df$latest_preindex_fi_date,
    role = "Case",
    Cohort = cohort_label,
    stringsAsFactors = FALSE
  )

  unmatched_cases <- cases_a %>%
    filter(riskset_size == 0)

  list(
    assignments = bind_rows(
      cases_a %>% filter(riskset_size > 0),
      controls
    ),
    unmatched_cases = unmatched_cases,
    n_eligible = n_eligible,
    fail_counts = bind_rows(fail_counts)
  )
}

build_riskset_matched_long <- function(input_path,
                                       classification_vars,
                                       classify_fn,
                                       cohort_levels,
                                       target_cycles,
                                       match_ratio,
                                       age_caliper,
                                       cycle_caliper = 1L,
                                       min_visits,
                                       seed,
                                       index_age_scaling = NULL,
                                       max_unmatched_fraction = 0.10,
                                       max_absolute_smd = 0.05) {
  prepped <- prepare_riskset_inputs(
    input_path = input_path,
    target_cycles = target_cycles,
    min_visits = min_visits,
    classification_vars = classification_vars
  )

  person_level <- prepped$person_level
  control_cancer_free_cutoff <- prepped$cancer_free_cutoff
  person_level$case_cohort <- classify_fn(person_level)

  observation_by_id <- split(
    prepped$observation_panel[, c("id", "cycle", "worked_rtmnyr", "fi_score_nocancer")],
    prepped$observation_panel$id
  )
  observation_by_id <- lapply(observation_by_id, function(obs) {
    obs <- obs[order(obs$worked_rtmnyr), , drop = FALSE]
    fi_seen <- !is.na(obs$fi_score_nocancer)
    distinct_fi <- logical(nrow(obs))
    fi_pos <- which(fi_seen)
    distinct_fi[fi_pos] <- !duplicated(obs$worked_rtmnyr[fi_pos])
    obs$cum_distinct_fi <- cumsum(distinct_fi)
    latest_marker <- ifelse(fi_seen, obs$worked_rtmnyr, -Inf)
    obs$cum_latest_fi_date <- cummax(latest_marker)
    obs$cum_latest_fi_date[!is.finite(obs$cum_latest_fi_date)] <- NA_real_
    obs
  })

  case_candidates <- person_level %>%
    filter(
      is_case == 1,
      !is.na(true_age_at_cancer),
      !is.na(case_cohort),
      !is.na(index_cycle)
    ) %>%
    arrange(cancer_dateca)

  case_support <- lapply(seq_len(nrow(case_candidates)), function(i) {
    preindex_support(
      observation_by_id[[case_candidates$id[[i]]]],
      case_candidates$cancer_dateca[[i]],
      case_candidates$index_cycle[[i]],
      cycle_caliper,
      min_visits
    )
  })
  if (nrow(case_candidates) > 0) {
    case_candidates <- case_candidates %>%
      mutate(
        case_active_cycle = vapply(case_support, `[[`, character(1), "active_cycle"),
        case_active_date = vapply(case_support, `[[`, numeric(1), "active_date"),
        case_cycle_gap = vapply(case_support, `[[`, numeric(1), "cycle_gap"),
        n_preindex_fi = vapply(case_support, `[[`, integer(1), "n_preindex_fi"),
        latest_preindex_fi_date = vapply(case_support, `[[`, numeric(1), "latest_preindex_fi_date"),
        has_active_return = vapply(case_support, `[[`, logical(1), "has_active_return"),
        active_cycle_ok = vapply(case_support, `[[`, logical(1), "cycle_ok"),
        preindex_fi_ok = vapply(case_support, `[[`, logical(1), "fi_ok"),
        active_at_index = vapply(case_support, `[[`, logical(1), "eligible"),
        exclusion_reason = case_when(
          !has_active_return ~ "no_participated_return_on_or_before_index",
          !active_cycle_ok ~ "latest_participated_cycle_outside_caliper",
          !preindex_fi_ok ~ "insufficient_preindex_fi",
          TRUE ~ NA_character_
        )
      )
  } else {
    case_candidates <- case_candidates %>%
      mutate(
        case_active_cycle = character(),
        case_active_date = numeric(),
        case_cycle_gap = numeric(),
        n_preindex_fi = integer(),
        latest_preindex_fi_date = numeric(),
        has_active_return = logical(),
        active_cycle_ok = logical(),
        preindex_fi_ok = logical(),
        active_at_index = logical(),
        exclusion_reason = character()
      )
  }
  inactive_cases <- case_candidates %>% filter(!active_at_index)
  cases <- case_candidates %>% filter(active_at_index)

  pool <- list(
    id = person_level$id,
    dob = person_level$dbmy09,
    dth = person_level$dtdth,
    canc = person_level$cancer_dateca,
    first = person_level$first_return
  )

  set.seed(seed)
  matched_by_cohort <- lapply(cohort_levels, function(cl) {
    cdf <- cases %>% filter(case_cohort == cl)
    match_one_cohort(
      case_df = cdf,
      pool = pool,
      observation_by_id = observation_by_id,
      ratio = match_ratio,
      age_caliper = age_caliper,
      cycle_caliper = cycle_caliper,
      min_visits = min_visits,
      cohort_label = cl,
      seed = seed,
      control_cancer_free_cutoff = control_cancer_free_cutoff
    )
  })
  names(matched_by_cohort) <- cohort_levels

  assignments <- bind_rows(lapply(matched_by_cohort, `[[`, "assignments"))
  unmatched_cases <- bind_rows(lapply(matched_by_cohort, `[[`, "unmatched_cases"))

  keep_cols <- c(
    "id", "cycle", "worked_rtmnyr", "age_at_cycle",
    "fi_score", "fi_score_nocancer", "fi_score_nocancer_nocarry",
    "frailty_cat", "frailty_cat_nocancer",
    "n_items_observed", "n_answered", "n_answered_nocancer",
    "items_asked", "items_asked_nocancer"
  )

  fi_keep_cols <- optional_cols(prepped$fi_trajectory, keep_cols)

  fi_rows <- prepped$fi_trajectory %>%
    semi_join(person_level, by = "id") %>%
    select(all_of(fi_keep_cols)) %>%
    left_join(prepped$baseline_covs, by = "id") %>%
    left_join(select(person_level, id, dbmy09, cancer_dateca, true_age_at_cancer), by = "id")

  assignment_ages <- assignments %>%
    left_join(select(person_level, id, dbmy09), by = "id") %>%
    mutate(index_age = (index_date - dbmy09) / 12)
  age_reference <- assignment_ages %>%
    filter(role == "Case") %>%
    distinct(Cohort, match_set, index_age)
  if (is.null(index_age_scaling)) {
    index_age_mean <- mean(age_reference$index_age, na.rm = TRUE)
    index_age_sd <- stats::sd(age_reference$index_age, na.rm = TRUE)
    scaling_source <- "estimated from this matched cohort"
  } else {
    required_scaling <- c("index_age_mean", "index_age_sd")
    missing_scaling <- setdiff(required_scaling, names(index_age_scaling))
    if (length(missing_scaling) > 0) {
      stop("Supplied index_age_scaling is missing: ",
           paste(missing_scaling, collapse = ", "), call. = FALSE)
    }
    index_age_mean <- as.numeric(index_age_scaling$index_age_mean[[1]])
    index_age_sd <- as.numeric(index_age_scaling$index_age_sd[[1]])
    scaling_source <-
      "reused from primary full-endpoint cancer-free matched cohort"
  }
  if (!is.finite(index_age_mean) || !is.finite(index_age_sd) || index_age_sd <= 0) {
    stop("Case-based index-age scaling constants are not finite and positive.")
  }
  scaling_metadata <- data.frame(
    reference_population = "one retained cancer case per matched set",
    n_reference_cases = nrow(age_reference),
    index_age_mean = index_age_mean,
    index_age_sd = index_age_sd,
    scaling_source = scaling_source,
    stringsAsFactors = FALSE
  )

  matched_long <- assignments %>%
    inner_join(fi_rows, by = "id", relationship = "many-to-many") %>%
    mutate(
      index_age    = (index_date - dbmy09) / 12,
      Age_Centered = age_at_cycle - index_age,
      Post         = if_else(Age_Centered > 0, 1L, 0L),
      Age_Post     = if_else(Post == 1L, Age_Centered, 0),
      Group        = if_else(role == "Case", "Cancer Case", "Control"),
      in_win_8     = abs(Age_Centered) <= 8,
      in_win_12    = abs(Age_Centered) <= 12,
      in_win_16    = abs(Age_Centered) <= 16,
      in_win_20    = abs(Age_Centered) <= 20,
      own_cancer_after_index = !is.na(cancer_dateca) & cancer_dateca > index_date,
      post_own_cancer = role == "Control" &
        own_cancer_after_index &
        !is.na(worked_rtmnyr) &
        worked_rtmnyr >= cancer_dateca,
      control_cancer_free_through_endpoint =
        role != "Control" |
        is.na(cancer_dateca) |
        cancer_dateca > control_cancer_free_cutoff,
      trajectory_id = paste(Cohort, match_set, id, role, sep = "__")
    ) %>%
    mutate(index_age_z = (index_age - index_age_mean) / index_age_sd) %>%
    mutate(
      Cohort    = factor(Cohort, levels = cohort_levels),
      Group     = factor(Group, levels = c("Control", "Cancer Case")),
      match_set = factor(match_set),
      id        = factor(id),
      role      = factor(role, levels = c("Control", "Case")),
      index_cycle = factor(index_cycle, levels = target_cycles)
    )

  bad_endpoint_controls <- matched_long %>%
    filter(
      role == "Control",
      !control_cancer_free_through_endpoint
    )
  if (nrow(bad_endpoint_controls) > 0L) {
    stop(
      "Cancer-free-through-endpoint assertion failed: retained controls have ",
      "a cancer date on or before the derived HPFS endpoint.",
      call. = FALSE
    )
  }

  duplicate_rows <- matched_long %>%
    count(trajectory_id, cycle) %>%
    filter(n > 1)

  if (nrow(duplicate_rows) > 0) {
    stop("Duplicated trajectory_id x cycle rows detected in matched output.")
  }

  set_integrity <- matched_long %>%
    distinct(Cohort, match_set, id, role, trajectory_id) %>%
    group_by(Cohort, match_set) %>%
    summarize(
      n_cases = sum(role == "Case"),
      n_controls = sum(role == "Control"),
      valid_set = n_cases == 1L & n_controls >= 1L,
      .groups = "drop"
    )
  if (any(!set_integrity$valid_set)) {
    stop("Post-expansion matched-set integrity failed: every set must contain exactly one case and at least one control.")
  }
  caliper_integrity <- assignments %>%
    summarize(
      max_abs_age_gap = max(abs(age_gap), na.rm = TRUE),
      max_cycle_gap_years = max(cycle_gap, na.rm = TRUE),
      all_age_ok = all(abs(age_gap) <= age_caliper),
      all_cycle_ok = all(cycle_gap <= 4L * cycle_caliper),
      all_preindex_fi_ok = all(n_preindex_fi >= min_visits),
      .groups = "drop"
    )
  if (!all(caliper_integrity$all_age_ok, caliper_integrity$all_cycle_ok,
           caliper_integrity$all_preindex_fi_ok)) {
    stop("Post-matching age, cycle, or pre-index FI integrity check failed.")
  }

  diagnostics <- bind_rows(lapply(cohort_levels, function(cl) {
    res <- matched_by_cohort[[cl]]
    data.frame(
      Cohort = cl,
      cases_identified = sum(case_candidates$case_cohort == cl),
      cases_excluded_inactive = sum(inactive_cases$case_cohort == cl),
      cases_attempted = length(res$n_eligible),
      cases_with_controls = sum(res$n_eligible > 0),
      cases_without_controls = sum(res$n_eligible == 0),
      control_records = sum(res$assignments$role == "Control"),
      distinct_controls = n_distinct(res$assignments$id[res$assignments$role == "Control"]),
      min_visits = min_visits,
      match_ratio = match_ratio,
      age_caliper = age_caliper,
      cycle_caliper_adjacent_cycles = cycle_caliper,
      seed = seed,
      stringsAsFactors = FALSE
    )
  }))

  fail_diagnostics <- bind_rows(lapply(matched_by_cohort, `[[`, "fail_counts"))
  reuse_diagnostics <- matched_long %>%
    filter(role == "Control") %>%
    distinct(Cohort, match_set, id, trajectory_id) %>%
    count(id, name = "control_assignment_count") %>%
    arrange(desc(control_assignment_count))

  post_own_cancer_diagnostics <- matched_long %>%
    filter(role == "Control") %>%
    summarize(
      control_assignments_with_later_own_cancer = n_distinct(trajectory_id[own_cancer_after_index]),
      rows_flagged_post_own_cancer = sum(post_own_cancer, na.rm = TRUE),
      .groups = "drop"
    )

  cancer_free_endpoint_diagnostics <- matched_long %>%
    filter(role == "Control") %>%
    summarize(
      control_assignments_with_cancer_on_or_before_endpoint =
        n_distinct(trajectory_id[!control_cancer_free_through_endpoint]),
      rows_flagged_cancer_on_or_before_endpoint =
        sum(!control_cancer_free_through_endpoint, na.rm = TRUE),
      cancer_free_cutoff_month = control_cancer_free_cutoff,
      .groups = "drop"
    )

  visit_support_diagnostics <- matched_long %>%
    distinct(Cohort, Group, trajectory_id, cycle, worked_rtmnyr) %>%
    count(Cohort, Group, trajectory_id, name = "n_analytic_fi_visits") %>%
    count(Cohort, Group, n_analytic_fi_visits, name = "n_trajectories")

  standardized_difference <- function(x_case, x_control) {
    denom <- sqrt((stats::var(x_case, na.rm = TRUE) + stats::var(x_control, na.rm = TRUE)) / 2)
    if (!is.finite(denom) || denom == 0) return(NA_real_)
    (mean(x_case, na.rm = TRUE) - mean(x_control, na.rm = TRUE)) / denom
  }
  balance_source <- assignment_ages %>%
    mutate(index_cycle_year = cycle_to_number(index_cycle))
  balance_diagnostics <- data.frame(
    variable = c("index_age", "index_cycle_year"),
    absolute_smd = abs(c(
      standardized_difference(balance_source$index_age[balance_source$role == "Case"],
                              balance_source$index_age[balance_source$role == "Control"]),
      standardized_difference(balance_source$index_cycle_year[balance_source$role == "Case"],
                              balance_source$index_cycle_year[balance_source$role == "Control"])
    )),
    stringsAsFactors = FALSE
  )
  gate_status <- diagnostics %>%
    transmute(
      Cohort,
      unmatched_fraction = if_else(cases_attempted > 0,
                                   cases_without_controls / cases_attempted,
                                   0),
      unmatched_ok = unmatched_fraction < max_unmatched_fraction
    ) %>%
    mutate(
      max_absolute_smd_observed = suppressWarnings(max(balance_diagnostics$absolute_smd,
                                                       na.rm = TRUE)),
      balance_ok = all(is.finite(balance_diagnostics$absolute_smd)) &&
        all(balance_diagnostics$absolute_smd < max_absolute_smd),
      gate_pass = unmatched_ok & balance_ok
    )
  if (any(!gate_status$gate_pass)) {
    stop(
      "Gate G4 failed: require unmatched fraction < ", max_unmatched_fraction,
      " and finite absolute SMDs < ", max_absolute_smd, ".",
      call. = FALSE
    )
  }
  run_metadata <- list(
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_path = normalizePath(input_path, mustWork = FALSE),
    input_md5 = unname(tools::md5sum(input_path)),
    target_cycles = target_cycles,
    match_ratio = match_ratio,
    age_caliper_years = age_caliper,
    cycle_caliper_adjacent_cycles = cycle_caliper,
    min_preindex_fi_visits = min_visits,
    max_unmatched_fraction = max_unmatched_fraction,
    max_absolute_smd = max_absolute_smd,
    seed = seed,
    unrestricted_relative_time = TRUE,
    control_cancer_free_rule = "through_full_hpfs_endpoint",
    control_cancer_free_cutoff = control_cancer_free_cutoff,
    control_cancer_free_cutoff_units = "months_since_1900",
    control_cancer_free_cutoff_source =
      "max(cancer_dateca) from the complete input dataset"
  )

  list(
    matched_long = matched_long,
    diagnostics = diagnostics,
    fail_diagnostics = fail_diagnostics,
    inactive_cases = inactive_cases,
    unmatched_cases = unmatched_cases,
    reuse_diagnostics = reuse_diagnostics,
    post_own_cancer_diagnostics = post_own_cancer_diagnostics,
    cancer_free_endpoint_diagnostics = cancer_free_endpoint_diagnostics,
    visit_support_diagnostics = visit_support_diagnostics,
    set_integrity = set_integrity,
    caliper_integrity = caliper_integrity,
    balance_diagnostics = balance_diagnostics,
    gate_status = gate_status,
    scaling_metadata = scaling_metadata,
    run_metadata = run_metadata
  )
}

save_riskset_match <- function(matched_result, output_path, label) {
  saveRDS(matched_result$matched_long, output_path)
  diagnostics_dir <- file.path(dirname(dirname(output_path)), "Results", "cancer", "data", "matching_diagnostics")
  if (!dir.exists(diagnostics_dir)) dir.create(diagnostics_dir, recursive = TRUE)
  output_stem <- tools::file_path_sans_ext(basename(output_path))
  write.csv(matched_result$diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_yield.csv")),
            row.names = FALSE)
  write.csv(matched_result$fail_diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_failed_criteria.csv")),
            row.names = FALSE)
  write.csv(matched_result$inactive_cases,
            file.path(diagnostics_dir, paste0(output_stem, "_inactive_cases.csv")),
            row.names = FALSE)
  write.csv(matched_result$unmatched_cases,
            file.path(diagnostics_dir, paste0(output_stem, "_unmatched_cases.csv")),
            row.names = FALSE)
  write.csv(matched_result$reuse_diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_control_reuse.csv")),
            row.names = FALSE)
  write.csv(matched_result$post_own_cancer_diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_post_own_cancer.csv")),
            row.names = FALSE)
  write.csv(matched_result$cancer_free_endpoint_diagnostics,
            file.path(diagnostics_dir,
                      paste0(output_stem, "_cancer_free_endpoint.csv")),
            row.names = FALSE)
  write.csv(matched_result$visit_support_diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_visit_support.csv")),
            row.names = FALSE)
  write.csv(matched_result$set_integrity,
            file.path(diagnostics_dir, paste0(output_stem, "_set_integrity.csv")),
            row.names = FALSE)
  write.csv(matched_result$caliper_integrity,
            file.path(diagnostics_dir, paste0(output_stem, "_caliper_integrity.csv")),
            row.names = FALSE)
  write.csv(matched_result$balance_diagnostics,
            file.path(diagnostics_dir, paste0(output_stem, "_balance.csv")),
            row.names = FALSE)
  write.csv(matched_result$gate_status,
            file.path(diagnostics_dir, paste0(output_stem, "_gate_g4.csv")),
            row.names = FALSE)
  saveRDS(matched_result$scaling_metadata,
          file.path(diagnostics_dir, paste0(output_stem, "_scaling_metadata.rds")))
  matched_result$run_metadata$output_path <- normalizePath(output_path, mustWork = FALSE)
  matched_result$run_metadata$output_md5 <- unname(tools::md5sum(output_path))
  saveRDS(matched_result$run_metadata,
          file.path(diagnostics_dir, paste0(output_stem, "_run_metadata.rds")))
  saveRDS(
    list(
      scaling = matched_result$scaling_metadata,
      run = matched_result$run_metadata,
      balance = matched_result$balance_diagnostics,
      set_integrity = matched_result$set_integrity,
      caliper_integrity = matched_result$caliper_integrity
    ),
    file.path(diagnostics_dir, paste0(output_stem, "_matching_metadata.rds"))
  )

  cat("\n", label, "\n", sep = "")
  cat("Saved matched dataset to: ", output_path, "\n", sep = "")
  cat("Matched analytic rows: ", nrow(matched_result$matched_long), "\n", sep = "")
  print(table(Cohort = matched_result$matched_long$Cohort,
              Group = matched_result$matched_long$Group))
  cat("\nMatching diagnostics:\n")
  print(matched_result$diagnostics)
  cat("\nDiagnostic CSVs written to: ", diagnostics_dir, "\n", sep = "")
}
