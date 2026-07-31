# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  5.0_ipcw_gee_functions.R
# Author:  Nemo Zhou
# Date started:      2026-07-27
# Date last updated: 2026-07-28
#
# Purpose:
#   Shared implementation of the inverse-probability-of-censoring-weighted
#   (IPCW) generalized estimating equation (GEE) sensitivity analysis for the
#   active risk-set matched cancer cohorts.  The engine implements only the
#   user-selected M2-equivalent full baseline-adjustment model from the active
#   4.x GLME workflow; it does not fit the M0, M1, or M3 model sets.
#
#   For each supplied matched cohort, this script:
#     1. verifies Gate G4 and the matched-data MD5 provenance;
#     2. builds a complete full-cohort response ledger from the canonical FI
#        panel, then transports its observation probabilities to a complete
#        matched-assignment x 4-year-cycle ledger, preserving reused controls;
#     3. treats death as a terminal event, never a censoring event;
#     4. estimates primary stabilized IPCW in the full-cohort response ledger
#        using calendar cycle, attained age, and cancer history in the numerator
#        and M2 baseline covariates, lagged FI, and response history in the
#        denominator;
#     5. applies monotone post-index nonresponse censoring and attaches
#        assignment-specific weights; and
#     6. fits unweighted and IPCW-weighted M2 GEE models with identical rows,
#        fixed effects, stored df-3 spline basis, and participant-id clusters.
#
#   The canonical spline knots are imported from the corresponding 4.x metadata
#   artifact.  This engine never re-estimates spline knots, never runs matching
#   or GLME models, and writes only data/RDS/CSV outputs (no persistent PNGs).
#
# Inputs supplied by a runner:
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Data/riskset_matched_*_long.rds
#   Results/cancer/data/4.x_spline_metadata.rds
#
# Outputs supplied by a runner:
#   Data/*_ipcw_m2_gee.rds
#   Results/cancer/data/5.x_* CSV/RDS diagnostic and model summaries
#
# Methods:
#   Documents/Methods/IPCW_Censoring_Weights.md
#   Documents/Methods/GLME_Natural_Spline_Trajectory_Analysis.md
# =============================================================================

library(dplyr)
library(geepack)
library(splines)
library(tidyr)

.ipcw_gee_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")
.ipcw_gee_factor_covars <- c("base_race", "base_marital", "base_living", "base_pckgr")
.ipcw_gee_numeric_covars <- c(
  "index_age_z", "base_calor", "base_sat", "base_diet_chol", "base_alco"
)
.ipcw_gee_m2_covars <- c(.ipcw_gee_factor_covars, .ipcw_gee_numeric_covars)

ipcw_gee_first_nonmissing <- function(x) {
  keep <- which(!is.na(x))
  if (length(keep)) x[[keep[[1]]]] else x[[1]]
}

ipcw_gee_locf <- function(x) {
  out <- x
  last <- NA_real_
  for (i in seq_along(out)) {
    if (!is.na(out[[i]])) last <- out[[i]]
    out[[i]] <- last
  }
  out
}

ipcw_gee_cycle_dates <- function(cycles = .ipcw_gee_cycles) {
  cycle_num <- as.integer(cycles)
  calendar_year <- ifelse(cycle_num >= 50L, 1900L + cycle_num, 2000L + cycle_num)
  tibble(
    cycle = as.character(cycles),
    cycle_order = seq_along(cycles),
    cycle_date = (calendar_year - 1900) * 12 + 6
  )
}

ipcw_gee_validate_matching_provenance <- function(matched_path) {
  if (!file.exists(matched_path)) {
    stop("Matched dataset not found: ", matched_path, call. = FALSE)
  }
  project_dir <- dirname(dirname(normalizePath(matched_path, mustWork = FALSE)))
  diagnostics_dir <- file.path(project_dir, "Results", "cancer", "data", "matching_diagnostics")
  output_stem <- tools::file_path_sans_ext(basename(matched_path))
  gate_path <- file.path(diagnostics_dir, paste0(output_stem, "_gate_g4.csv"))
  run_path <- file.path(diagnostics_dir, paste0(output_stem, "_run_metadata.rds"))

  if (!file.exists(gate_path) || !file.exists(run_path)) {
    stop("Gate G4 or matching provenance is missing for ", output_stem,
         ". Rerun the corresponding 2.x builder before IPCW-GEE.", call. = FALSE)
  }

  gate <- read.csv(gate_path, stringsAsFactors = FALSE)
  if (!"gate_pass" %in% names(gate) || !all(as.logical(gate$gate_pass) %in% TRUE)) {
    stop("Gate G4 did not pass for ", output_stem,
         ". IPCW-GEE will not use an unaccepted matched cohort.", call. = FALSE)
  }

  run_metadata <- readRDS(run_path)
  expected_md5 <- unname(as.character(run_metadata$output_md5))
  observed_md5 <- unname(tools::md5sum(matched_path))
  if (length(expected_md5) != 1L || !nzchar(expected_md5) ||
      !identical(expected_md5, observed_md5)) {
    stop("Matched RDS MD5 does not match Gate G4 provenance for ", output_stem,
         ". Rerun the corresponding 2.x builder.", call. = FALSE)
  }

  list(
    output_stem = output_stem,
    gate_path = gate_path,
    run_path = run_path,
    input_md5 = observed_md5
  )
}

ipcw_gee_read_spline_metadata <- function(metadata_path, cohort) {
  if (!file.exists(metadata_path)) {
    stop("Required 4.x spline metadata is missing: ", metadata_path,
         ". Run the corresponding 4.x GLME wrapper before IPCW-GEE.", call. = FALSE)
  }
  metadata <- readRDS(metadata_path)
  cohort_metadata <- if (!is.null(metadata$cohort)) {
    metadata
  } else if (cohort %in% names(metadata)) {
    metadata[[cohort]]
  } else {
    stop("Spline metadata does not contain cohort '", cohort, "': ", metadata_path,
         call. = FALSE)
  }

  basis <- cohort_metadata$spline_bases$adjusted_spline
  required <- c("terms", "knots", "boundary_knots", "spline_df")
  missing <- setdiff(required, names(basis))
  if (length(missing)) {
    stop("Malformed adjusted-spline metadata; missing ",
         paste(missing, collapse = ", "), " in ", metadata_path, call. = FALSE)
  }
  if (!identical(as.integer(basis$spline_df), 3L)) {
    stop("IPCW-GEE M2 requires the active df-3 adjusted spline; metadata has df = ",
         basis$spline_df, call. = FALSE)
  }

  list(
    cohort_metadata = cohort_metadata,
    terms = as.character(basis$terms),
    knots = as.numeric(basis$knots),
    boundary_knots = as.numeric(basis$boundary_knots),
    spline_df = as.integer(basis$spline_df),
    metadata_path = normalizePath(metadata_path, mustWork = TRUE)
  )
}

ipcw_gee_add_spline_columns <- function(data, spline_info) {
  spline_basis <- splines::ns(
    data$Age_Centered,
    knots = spline_info$knots,
    Boundary.knots = spline_info$boundary_knots
  )
  if (ncol(spline_basis) != length(spline_info$terms)) {
    stop("Stored spline term count does not match the reconstructed basis.", call. = FALSE)
  }
  spline_data <- as.data.frame(unclass(spline_basis))
  names(spline_data) <- spline_info$terms
  bind_cols(data, spline_data)
}

ipcw_gee_most_common_level <- function(x) {
  x <- as.character(x[!is.na(x)])
  if (!length(x)) stop("No observed values available for a prediction factor.", call. = FALSE)
  tab <- sort(table(x), decreasing = TRUE)
  candidates <- names(tab)[tab == tab[[1]]]
  sort(candidates)[[1]]
}

ipcw_gee_fit_binomial <- function(formula, data, label) {
  separation_warning <- FALSE
  fit <- withCallingHandlers(
    stats::glm(formula, data = data, family = stats::binomial()),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1|did not converge", conditionMessage(w),
                ignore.case = TRUE)) {
        separation_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )

  method <- "glm"
  if (separation_warning) {
    if (!requireNamespace("brglm2", quietly = TRUE)) {
      stop(label, " showed separation/convergence warnings, but brglm2 is unavailable.",
           call. = FALSE)
    }
    fit <- stats::glm(
      formula,
      data = data,
      family = stats::binomial(),
      method = brglm2::brglmFit
    )
    method <- "brglm2::brglmFit"
  }
  list(fit = fit, method = method, separation_warning = separation_warning)
}

ipcw_gee_robust_vcov <- function(fit) {
  summary_fit <- summary(fit)
  V <- summary_fit$cov.scaled
  if (is.null(V)) V <- fit$geese$vbeta
  V <- as.matrix(V)
  beta_names <- names(stats::coef(fit))
  if (!identical(dim(V), c(length(beta_names), length(beta_names)))) {
    stop("Could not align the GEE empirical covariance with its coefficients.", call. = FALSE)
  }
  dimnames(V) <- list(beta_names, beta_names)
  V
}

ipcw_gee_make_reference_grid <- function(gee_data, times, spline_info) {
  factor_covars <- .ipcw_gee_factor_covars
  numeric_covars <- .ipcw_gee_numeric_covars
  grid <- expand.grid(
    Age_Centered = times,
    Group = factor(c("Control", "Cancer Case"), levels = c("Control", "Cancer Case")),
    stringsAsFactors = FALSE
  )

  for (covar in factor_covars) {
    selected <- ipcw_gee_most_common_level(gee_data[[covar]])
    grid[[covar]] <- factor(selected, levels = levels(gee_data[[covar]]))
  }
  for (covar in numeric_covars) {
    grid[[covar]] <- if (identical(covar, "index_age_z")) {
      0
    } else {
      mean(gee_data[[covar]], na.rm = TRUE)
    }
  }
  ipcw_gee_add_spline_columns(grid, spline_info)
}

ipcw_gee_full_baseline <- function(raw) {
  required <- c("race", "marital_status", "living_arr", "pckgr", "base_living",
                "base_calor", "base_sat", "base_diet_chol", "base_alco")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("Canonical FI panel is missing baseline variables needed for the full-cohort IPCW model: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  raw %>%
    mutate(id = as.character(id), cycle = as.character(cycle)) %>%
    filter(cycle == "86") %>%
    distinct(id, .keep_all = TRUE) %>%
    transmute(
      id,
      base_race = if_else(is.na(race), "Missing", as.character(race)),
      base_marital = if_else(is.na(marital_status), "Missing", as.character(marital_status)),
      base_living = if_else(
        is.na(base_living),
        if_else(is.na(living_arr), "Missing", as.character(living_arr)),
        as.character(base_living)
      ),
      base_pckgr = if_else(is.na(pckgr), "Missing", as.character(pckgr)),
      base_calor = as.numeric(as.character(base_calor)),
      base_sat = as.numeric(as.character(base_sat)),
      base_diet_chol = as.numeric(as.character(base_diet_chol)),
      base_alco = as.numeric(as.character(base_alco))
    )
}

ipcw_gee_prepare_full_model_frame <- function(data) {
  factor_covars <- .ipcw_gee_factor_covars
  numeric_covars <- setdiff(.ipcw_gee_numeric_covars, "index_age_z")
  data <- data %>% mutate(cycle_f = factor(cycle, levels = .ipcw_gee_cycles))

  for (covar in factor_covars) {
    data[[covar]] <- factor(if_else(is.na(data[[covar]]), "Missing",
                                    as.character(data[[covar]])))
  }

  numeric_to_impute <- c(numeric_covars, "lag_fi", "prior_response_count")
  medians <- vapply(numeric_to_impute, function(covar) {
    value <- stats::median(data[[covar]], na.rm = TRUE)
    if (is.finite(value)) value else 0
  }, numeric(1))
  for (covar in numeric_to_impute) {
    data[[covar]] <- ifelse(is.na(data[[covar]]), medians[[covar]], data[[covar]])
  }

  list(
    data = data,
    factor_levels = lapply(factor_covars, function(covar) levels(data[[covar]])) %>%
      stats::setNames(factor_covars),
    numeric_medians = medians
  )
}

ipcw_gee_fit_full_response_model <- function(panel_path) {
  if (!file.exists(panel_path)) {
    stop("Canonical FI panel not found: ", panel_path, call. = FALSE)
  }
  raw <- readRDS(panel_path)
  required <- c("id", "cycle", "participated", "fi_score_nocancer", "dbmy09", "dtdth",
                "cancer_index_dateca", "worked_rtmnyr", "age_at_cycle")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("Canonical FI panel is missing full-cohort IPCW columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  baseline <- ipcw_gee_full_baseline(raw)
  raw_4yr <- raw %>%
    transmute(
      id = as.character(id),
      cycle = as.character(cycle),
      participated = as.numeric(as.character(participated)),
      fi_score_nocancer = as.numeric(as.character(fi_score_nocancer)),
      dbmy09 = as.numeric(as.character(dbmy09)),
      dtdth = as.numeric(as.character(dtdth)),
      cancer_index_dateca = as.numeric(as.character(cancer_index_dateca)),
      age_at_cycle = as.numeric(as.character(age_at_cycle))
    ) %>%
    filter(cycle %in% .ipcw_gee_cycles) %>%
    distinct(id, cycle, .keep_all = TRUE)

  people <- raw_4yr %>%
    group_by(id) %>%
    summarize(
      dbmy09 = ipcw_gee_first_nonmissing(dbmy09),
      dtdth = suppressWarnings(max(dtdth, na.rm = TRUE)),
      cancer_index_dateca = suppressWarnings(min(cancer_index_dateca, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      dtdth = if_else(is.infinite(dtdth), NA_real_, dtdth),
      cancer_index_dateca = if_else(is.infinite(cancer_index_dateca), NA_real_, cancer_index_dateca)
    ) %>%
    filter(!is.na(dbmy09))

  observed_cycles <- raw_4yr %>%
    group_by(id, cycle) %>%
    summarize(
      observed = as.integer(any(participated == 1 & !is.na(fi_score_nocancer))),
      fi_score_nocancer = ipcw_gee_first_nonmissing(fi_score_nocancer),
      age_at_cycle = ipcw_gee_first_nonmissing(age_at_cycle),
      .groups = "drop"
    )

  ledger <- tidyr::crossing(id = people$id, ipcw_gee_cycle_dates()) %>%
    left_join(people, by = "id") %>%
    left_join(observed_cycles, by = c("id", "cycle")) %>%
    left_join(baseline, by = "id") %>%
    mutate(
      observed = coalesce(observed, 0L),
      age_nom = (cycle_date - dbmy09) / 12,
      alive_at_cycle = is.na(dtdth) | dtdth >= cycle_date,
      cancer_now = as.integer(!is.na(cancer_index_dateca) & cancer_index_dateca <= cycle_date),
      years_since_cancer = if_else(
        cancer_now == 1L,
        (cycle_date - cancer_index_dateca) / 12,
        0
      )
    ) %>%
    group_by(id) %>%
    arrange(cycle_order, .by_group = TRUE) %>%
    mutate(
      observed_fi_for_lag = if_else(observed == 1L, fi_score_nocancer, NA_real_),
      last_observed_fi = ipcw_gee_locf(observed_fi_for_lag),
      lag_fi = dplyr::lag(last_observed_fi),
      entered_prior = dplyr::lag(cumsum(observed) > 0L, default = FALSE),
      prior_observed = dplyr::lag(observed, default = 0L),
      prior_response_count = dplyr::lag(cumsum(observed), default = 0L),
      candidate_at_risk = alive_at_cycle & entered_prior & prior_observed == 1L,
      raw_nonresponse = candidate_at_risk & observed == 0L,
      prior_nonresponse = dplyr::lag(cummax(as.integer(raw_nonresponse)) > 0L,
                                      default = FALSE),
      at_risk = candidate_at_risk & !prior_nonresponse
    ) %>%
    ungroup()

  model_rows <- ledger %>% filter(at_risk)
  if (length(unique(model_rows$observed)) < 2L) {
    stop("The full-cohort response ledger has no observed/nonobserved variation among alive at-risk cycles.",
         call. = FALSE)
  }

  prepared <- ipcw_gee_prepare_full_model_frame(model_rows)
  model_data <- prepared$data
  factor_keep <- c("cycle_f", .ipcw_gee_factor_covars)
  factor_keep <- factor_keep[vapply(factor_keep, function(covar) {
    nlevels(model_data[[covar]]) >= 2L
  }, logical(1))]
  numeric_candidates <- c("cancer_now", "years_since_cancer",
                          setdiff(.ipcw_gee_numeric_covars, "index_age_z"),
                          "lag_fi", "prior_response_count")
  numeric_keep <- numeric_candidates[vapply(numeric_candidates, function(covar) {
    stats::sd(model_data[[covar]], na.rm = TRUE) > 0
  }, logical(1))]

  age_term <- "splines::ns(age_nom, df = 3)"
  numerator_terms <- c(age_term, intersect(c("cycle_f", "cancer_now", "years_since_cancer"),
                                            c(factor_keep, numeric_keep)))
  denominator_terms <- c(age_term, factor_keep, numeric_keep)
  numerator_formula <- stats::as.formula(paste("observed ~", paste(numerator_terms, collapse = " + ")))
  denominator_formula <- stats::as.formula(paste("observed ~", paste(denominator_terms, collapse = " + ")))

  numerator_fit <- ipcw_gee_fit_binomial(numerator_formula, model_data,
                                          "Full-cohort IPCW numerator")
  denominator_fit <- ipcw_gee_fit_binomial(denominator_formula, model_data,
                                            "Full-cohort IPCW denominator")

  list(
    ledger = ledger,
    model_data = model_data,
    numerator_fit = numerator_fit,
    denominator_fit = denominator_fit,
    numerator_formula = numerator_formula,
    denominator_formula = denominator_formula,
    factor_levels = prepared$factor_levels,
    numeric_medians = prepared$numeric_medians,
    population = "full-cohort person-cycle response ledger"
  )
}

ipcw_gee_build_assignment_panel <- function(matched_path, cohort, response_model) {
  matched <- readRDS(matched_path)
  matched_required <- c("id", "Cohort", "Group", "match_set", "role", "index_date", "index_age",
                        .ipcw_gee_m2_covars)
  missing_matched <- setdiff(matched_required, names(matched))
  if (length(missing_matched)) {
    stop("Matched input is missing M2 IPCW-GEE columns: ",
         paste(missing_matched, collapse = ", "), call. = FALSE)
  }

  assignment_columns <- c("Cohort", "Group", "match_set", "id", "role", "index_date", "index_age",
                          .ipcw_gee_m2_covars)
  assignment_meta <- matched %>%
    mutate(
      id = as.character(id), Cohort = as.character(Cohort), Group = as.character(Group),
      match_set = as.character(match_set), role = as.character(role),
      index_date = as.numeric(index_date), index_age = as.numeric(index_age)
    ) %>%
    filter(Cohort == cohort) %>%
    group_by(Cohort, Group, match_set, id, role) %>%
    summarize(
      across(all_of(setdiff(assignment_columns, c("Cohort", "Group", "match_set", "id", "role"))),
             ipcw_gee_first_nonmissing),
      .groups = "drop"
    ) %>%
    mutate(
      trajectory_id = paste(Cohort, match_set, id, role, sep = "__"),
      Group = factor(Group, levels = c("Control", "Cancer Case"))
    )
  if (!nrow(assignment_meta)) {
    stop("Matched input contains no rows for requested cohort '", cohort, "'.", call. = FALSE)
  }
  if (anyDuplicated(assignment_meta$trajectory_id)) {
    stop("Matched assignment metadata has duplicate trajectory_id values.", call. = FALSE)
  }

  response_columns <- response_model$ledger %>%
    select(id, cycle, cycle_order, cycle_date, fi_score_nocancer, age_at_cycle, age_nom,
           alive_at_cycle, cancer_index_dateca, cancer_now, years_since_cancer,
           lag_fi, prior_response_count, observed)
  assignment_meta %>%
    tidyr::crossing(ipcw_gee_cycle_dates()) %>%
    left_join(response_columns, by = c("id", "cycle", "cycle_order", "cycle_date")) %>%
    mutate(
      Age_Centered = coalesce(age_at_cycle - index_age, age_nom - index_age),
      post_own_cancer = as.character(Group) == "Control" &
        !is.na(cancer_index_dateca) & cancer_index_dateca > index_date &
        cycle_date >= cancer_index_dateca,
      assignment_eligible = !post_own_cancer,
      fi_observed = if_else(alive_at_cycle & assignment_eligible, observed == 1L, NA)
    ) %>%
    arrange(trajectory_id, cycle_order)
}

ipcw_gee_prepare_assignment_prediction_frame <- function(data, response_model) {
  out <- data
  for (covar in .ipcw_gee_factor_covars) {
    values <- if_else(is.na(out[[covar]]), "Missing", as.character(out[[covar]]))
    levels_expected <- response_model$factor_levels[[covar]]
    unknown <- setdiff(unique(values), levels_expected)
    if (length(unknown)) {
      stop("Matched assignment contains factor levels absent from the full-cohort response model for ",
           covar, ": ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    out[[covar]] <- factor(values, levels = levels_expected)
  }
  for (covar in c(setdiff(.ipcw_gee_numeric_covars, "index_age_z"), "lag_fi", "prior_response_count")) {
    out[[covar]] <- ifelse(is.na(out[[covar]]), response_model$numeric_medians[[covar]], out[[covar]])
  }
  out$cycle_f <- factor(out$cycle, levels = .ipcw_gee_cycles)
  out
}

ipcw_gee_prepare_weights <- function(panel, spline_info, response_model,
                                     truncation = c(0.01, 0.99)) {
  if (length(truncation) != 2L || any(!is.finite(truncation)) ||
      truncation[[1]] <= 0 || truncation[[2]] >= 1 || truncation[[1]] >= truncation[[2]]) {
    stop("truncation must contain two increasing probabilities strictly between 0 and 1.", call. = FALSE)
  }

  panel <- panel %>%
    mutate(
      across(all_of(.ipcw_gee_factor_covars), ~ factor(.x)),
      across(all_of(.ipcw_gee_numeric_covars), as.numeric)
    ) %>%
    group_by(trajectory_id) %>%
    arrange(cycle_order, .by_group = TRUE) %>%
    mutate(
      post_index = !is.na(Age_Centered) & Age_Centered >= 0,
      candidate_at_risk = post_index & alive_at_cycle & assignment_eligible & !is.na(lag_fi),
      raw_nonresponse = candidate_at_risk & !dplyr::coalesce(fi_observed, FALSE),
      prior_nonresponse = dplyr::lag(cummax(as.integer(raw_nonresponse)) > 0, default = FALSE),
      at_risk = candidate_at_risk & !prior_nonresponse,
      censor_event = at_risk & !dplyr::coalesce(fi_observed, FALSE),
      later_observed = rev(cummax(rev(as.integer(dplyr::coalesce(fi_observed, FALSE))))),
      returner_after_censor_event = censor_event & later_observed == 1L,
      retained_post = at_risk & dplyr::coalesce(fi_observed, FALSE),
      retained_pre = !post_index & alive_at_cycle & assignment_eligible &
        dplyr::coalesce(fi_observed, FALSE),
      retained_gee_row = retained_pre | retained_post
    ) %>%
    ungroup()

  panel <- ipcw_gee_add_spline_columns(panel, spline_info)
  weight_data <- panel %>% filter(at_risk)
  if (!nrow(weight_data) || length(unique(weight_data$Group)) < 2L) {
    stop("Insufficient two-group assignment-cycles for attached IPCW probabilities.", call. = FALSE)
  }

  prediction_data <- ipcw_gee_prepare_assignment_prediction_frame(weight_data, response_model)
  weight_data <- weight_data %>%
    mutate(
      p_num = as.numeric(stats::predict(response_model$numerator_fit$fit,
                                        newdata = prediction_data, type = "response")),
      p_den = as.numeric(stats::predict(response_model$denominator_fit$fit,
                                        newdata = prediction_data, type = "response"))
    )
  if (any(!is.finite(weight_data$p_num)) || any(!is.finite(weight_data$p_den)) ||
      any(weight_data$p_num <= 0) || any(weight_data$p_den <= 0)) {
    stop("Transported IPCW fitted probabilities are non-finite or non-positive.", call. = FALSE)
  }

  panel <- panel %>%
    left_join(weight_data %>% select(trajectory_id, cycle, p_num, p_den),
              by = c("trajectory_id", "cycle")) %>%
    group_by(trajectory_id) %>%
    arrange(cycle_order, .by_group = TRUE) %>%
    mutate(
      step_ratio = if_else(retained_post, p_num / p_den, 1),
      sw_ipcw_raw = cumprod(dplyr::coalesce(step_ratio, 1))
    ) %>%
    ungroup()

  retained_post_weights <- panel %>%
    filter(retained_post, is.finite(sw_ipcw_raw), sw_ipcw_raw > 0) %>%
    pull(sw_ipcw_raw)
  if (!length(retained_post_weights)) {
    stop("No retained post-index rows received valid transported IPCW weights.", call. = FALSE)
  }
  truncation_bounds <- stats::quantile(retained_post_weights, probs = truncation,
                                        na.rm = TRUE, names = FALSE)
  panel <- panel %>%
    mutate(
      sw_ipcw = case_when(
        !post_index ~ 1,
        retained_post ~ pmin(pmax(sw_ipcw_raw, truncation_bounds[[1]]), truncation_bounds[[2]]),
        TRUE ~ NA_real_
      )
    )

  list(
    panel = panel,
    weight_data = weight_data,
    numerator_fit = response_model$numerator_fit,
    denominator_fit = response_model$denominator_fit,
    numerator_formula = response_model$numerator_formula,
    denominator_formula = response_model$denominator_formula,
    truncation_bounds = truncation_bounds,
    model_population = response_model$population,
    full_response_summary = tibble(
      model_population = response_model$population,
      n_at_risk_cycles = nrow(response_model$model_data),
      n_observed = sum(response_model$model_data$observed == 1L),
      n_nonobserved = sum(response_model$model_data$observed == 0L)
    )
  )
}

ipcw_gee_fit_pair <- function(panel, spline_info, prediction_window = c(-16, 16)) {
  required <- c("fi_score_nocancer", "Age_Centered", "sw_ipcw", .ipcw_gee_m2_covars,
                spline_info$terms)
  gee_data <- panel %>%
    filter(retained_gee_row) %>%
    filter(if_all(all_of(required), ~ !is.na(.x))) %>%
    mutate(
      id = factor(id),
      Group = factor(Group, levels = c("Control", "Cancer Case")),
      across(all_of(.ipcw_gee_factor_covars), ~ droplevels(factor(.x)))
    ) %>%
    arrange(id, cycle_order, trajectory_id) %>%
    droplevels()

  if (!nrow(gee_data) || length(unique(gee_data$Group)) < 2L) {
    stop("Insufficient two-group complete-case rows for the M2 GEE pair.", call. = FALSE)
  }
  if (any(!is.finite(gee_data$sw_ipcw)) || any(gee_data$sw_ipcw <= 0)) {
    stop("The retained M2 GEE data contain missing or non-positive IPCW weights.", call. = FALSE)
  }

  time_rhs <- paste0("Group * (", paste(spline_info$terms, collapse = " + "), ")")
  outcome_rhs <- stats::as.formula(paste("~", time_rhs, "+",
                                         paste(.ipcw_gee_m2_covars, collapse = " + ")))
  outcome_formula <- stats::update.formula(outcome_rhs, fi_score_nocancer ~ .)

  gee_unweighted <- geepack::geeglm(
    outcome_formula,
    id = id,
    data = gee_data,
    weights = rep(1, nrow(gee_data)),
    corstr = "independence",
    std.err = "san.se"
  )
  gee_ipcw <- geepack::geeglm(
    outcome_formula,
    id = id,
    data = gee_data,
    weights = sw_ipcw,
    corstr = "independence",
    std.err = "san.se"
  )

  make_output <- function(fit, model_label) {
    beta <- stats::coef(fit)
    V <- ipcw_gee_robust_vcov(fit)
    coefficient_table <- tibble(
      Model = model_label,
      Term = names(beta),
      Estimate = as.numeric(beta),
      SE = sqrt(diag(V)),
      CI_low = Estimate - 1.96 * SE,
      CI_high = Estimate + 1.96 * SE,
      vcov_type = "GEE empirical sandwich"
    )

    prediction_grid <- ipcw_gee_make_reference_grid(
      gee_data,
      seq(prediction_window[[1]], prediction_window[[2]], by = 0.25),
      spline_info
    )
    X <- stats::model.matrix(outcome_rhs, prediction_grid)
    missing_terms <- setdiff(names(beta), colnames(X))
    if (length(missing_terms)) {
      stop("Prediction design matrix is missing GEE coefficients: ",
           paste(missing_terms, collapse = ", "), call. = FALSE)
    }
    X <- X[, names(beta), drop = FALSE]
    prediction_grid$pred <- as.vector(X %*% beta)
    prediction_grid$se <- sqrt(rowSums((X %*% V) * X))
    prediction_grid$lwr <- prediction_grid$pred - 1.96 * prediction_grid$se
    prediction_grid$upr <- prediction_grid$pred + 1.96 * prediction_grid$se
    prediction_grid$Model <- model_label
    prediction_grid$vcov_type <- "GEE empirical sandwich"

    case_grid <- prediction_grid %>%
      filter(Group == "Cancer Case") %>%
      arrange(Age_Centered)
    control_grid <- prediction_grid %>%
      filter(Group == "Control") %>%
      arrange(Age_Centered)
    X_case <- stats::model.matrix(outcome_rhs, case_grid)[, names(beta), drop = FALSE]
    X_control <- stats::model.matrix(outcome_rhs, control_grid)[, names(beta), drop = FALSE]
    X_difference <- X_case - X_control
    contrast <- tibble(
      Model = model_label,
      Age_Centered = case_grid$Age_Centered,
      contrast = "Cancer Case - Control",
      estimate = as.vector(X_difference %*% beta),
      se = sqrt(rowSums((X_difference %*% V) * X_difference)),
      vcov_type = "GEE empirical sandwich"
    ) %>%
      mutate(CI_low = estimate - 1.96 * se, CI_high = estimate + 1.96 * se)

    time_mid <- (contrast$Age_Centered[-1] + contrast$Age_Centered[-nrow(contrast)]) / 2
    X_derivative <- (X_difference[-1, , drop = FALSE] -
      X_difference[-nrow(X_difference), , drop = FALSE]) /
      diff(contrast$Age_Centered)
    derivative <- tibble(
      Model = model_label,
      Age_Centered = time_mid,
      estimate = as.vector(X_derivative %*% beta),
      se = sqrt(rowSums((X_derivative %*% V) * X_derivative)),
      vcov_type = "GEE empirical sandwich"
    ) %>%
      mutate(CI_low = estimate - 1.96 * se, CI_high = estimate + 1.96 * se)

    pre_rows <- derivative$Age_Centered >= -8 & derivative$Age_Centered < 0
    post_rows <- derivative$Age_Centered > 0 & derivative$Age_Centered <= 8
    if (!any(pre_rows) || !any(post_rows)) {
      stop("The prediction grid does not contain the required -8 to 0 and 0 to 8 theta intervals.",
           call. = FALSE)
    }
    theta_l <- colMeans(X_derivative[post_rows, , drop = FALSE]) -
      colMeans(X_derivative[pre_rows, , drop = FALSE])
    theta_estimate <- as.numeric(theta_l %*% beta)
    theta_se <- sqrt(as.numeric(theta_l %*% V %*% theta_l))
    theta <- tibble(
      Model = model_label,
      theta = theta_estimate,
      SE = theta_se,
      CI_low = theta_estimate - 1.96 * theta_se,
      CI_high = theta_estimate + 1.96 * theta_se,
      z = theta_estimate / theta_se,
      p_value = 2 * stats::pnorm(abs(theta_estimate / theta_se), lower.tail = FALSE),
      vcov_type = "GEE empirical sandwich"
    )

    list(
      fit = fit,
      coefficients = coefficient_table,
      predictions = prediction_grid,
      contrast = contrast,
      derivative = derivative,
      theta = theta
    )
  }

  unweighted <- make_output(gee_unweighted, "Unweighted GEE")
  weighted <- make_output(gee_ipcw, "IPCW-weighted GEE")
  list(
    gee_data = gee_data,
    outcome_formula = outcome_formula,
    unweighted = unweighted,
    weighted = weighted
  )
}

ipcw_gee_write_outputs <- function(result, out_prefix, results_dir, data_out_path,
                                   cohort, matched_path, spline_info, provenance) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  data_dir <- dirname(data_out_path)
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  panel <- result$weights$panel
  gee_data <- result$gee$gee_data
  saveRDS(panel, data_out_path)

  coefficient_table <- bind_rows(
    result$gee$unweighted$coefficients,
    result$gee$weighted$coefficients
  ) %>% mutate(Cohort = cohort, .before = 1)
  predictions <- bind_rows(result$gee$unweighted$predictions, result$gee$weighted$predictions) %>%
    mutate(Cohort = cohort, .before = 1)
  contrast <- bind_rows(result$gee$unweighted$contrast, result$gee$weighted$contrast) %>%
    mutate(Cohort = cohort, .before = 1)
  derivative <- bind_rows(result$gee$unweighted$derivative, result$gee$weighted$derivative) %>%
    mutate(Cohort = cohort, .before = 1)
  theta <- bind_rows(result$gee$unweighted$theta, result$gee$weighted$theta) %>%
    mutate(Cohort = cohort, .before = 1)

  censoring_events <- panel %>%
    group_by(Cohort, Group) %>%
    summarize(
      n_assignment_cycles = n(),
      n_distinct_ids = n_distinct(id),
      n_alive_eligible_cycles = sum(alive_at_cycle & assignment_eligible, na.rm = TRUE),
      n_at_risk = sum(at_risk, na.rm = TRUE),
      n_nonresponse_censor_events = sum(censor_event, na.rm = TRUE),
      n_returner_assignments = sum(returner_after_censor_event, na.rm = TRUE),
      returner_percent = dplyr::if_else(
        n_nonresponse_censor_events > 0,
        100 * n_returner_assignments / n_nonresponse_censor_events,
        NA_real_
      ),
      n_post_death_cycles_excluded = sum(!alive_at_cycle, na.rm = TRUE),
      n_post_own_cancer_cycles_excluded = sum(post_own_cancer, na.rm = TRUE),
      n_retained_gee_rows = sum(retained_gee_row, na.rm = TRUE),
      .groups = "drop"
    )

  retained_weights <- panel %>%
    filter(retained_gee_row) %>%
    group_by(Cohort, Group) %>%
    summarize(
      n_rows = n(),
      n_distinct_ids = n_distinct(id),
      mean_weight = mean(sw_ipcw, na.rm = TRUE),
      sd_weight = sd(sw_ipcw, na.rm = TRUE),
      min_weight = min(sw_ipcw, na.rm = TRUE),
      q01_weight = stats::quantile(sw_ipcw, 0.01, na.rm = TRUE, names = FALSE),
      q25_weight = stats::quantile(sw_ipcw, 0.25, na.rm = TRUE, names = FALSE),
      q50_weight = stats::quantile(sw_ipcw, 0.50, na.rm = TRUE, names = FALSE),
      q75_weight = stats::quantile(sw_ipcw, 0.75, na.rm = TRUE, names = FALSE),
      q99_weight = stats::quantile(sw_ipcw, 0.99, na.rm = TRUE, names = FALSE),
      max_weight = max(sw_ipcw, na.rm = TRUE),
      effective_sample_size = sum(sw_ipcw, na.rm = TRUE)^2 / sum(sw_ipcw^2, na.rm = TRUE),
      .groups = "drop"
    )

  positivity <- result$weights$weight_data %>%
    mutate(
      rel_time_bin = cut(Age_Centered, breaks = c(-Inf, -8, -4, 0, 4, 8, Inf),
                         right = FALSE)
    ) %>%
    group_by(Group, rel_time_bin, .drop = FALSE) %>%
    summarize(
      n_assignment_cycles = n(),
      n_pden_lt_005 = sum(p_den < 0.05),
      percent_pden_lt_005 = 100 * n_pden_lt_005 / n_assignment_cycles,
      positivity_warning = percent_pden_lt_005 > 1,
      .groups = "drop"
    ) %>%
    mutate(Cohort = cohort, .before = 1)

  summarize_weight_model <- function(fit, model_label) {
    coefficient_matrix <- as.data.frame(summary(fit)$coefficients)
    tibble(
      model = model_label,
      term = rownames(coefficient_matrix),
      estimate = coefficient_matrix[[1]],
      std_error = coefficient_matrix[[2]],
      z_value = coefficient_matrix[[3]],
      p_value = coefficient_matrix[[4]]
    )
  }
  weight_coefficients <- bind_rows(
    summarize_weight_model(result$weights$numerator_fit$fit, "numerator"),
    summarize_weight_model(result$weights$denominator_fit$fit, "denominator")
  ) %>% mutate(Cohort = cohort, .before = 1)

  write.csv(censoring_events, file.path(results_dir, paste0(out_prefix, "_censoring_events.csv")),
            row.names = FALSE)
  write.csv(retained_weights, file.path(results_dir, paste0(out_prefix, "_weight_summary.csv")),
            row.names = FALSE)
  write.csv(positivity, file.path(results_dir, paste0(out_prefix, "_positivity.csv")),
            row.names = FALSE)
  write.csv(result$weights$full_response_summary,
            file.path(results_dir, paste0(out_prefix, "_full_response_model_summary.csv")),
            row.names = FALSE)
  write.csv(weight_coefficients,
            file.path(results_dir, paste0(out_prefix, "_weight_model_coefficients.csv")),
            row.names = FALSE)
  write.csv(coefficient_table, file.path(results_dir, paste0(out_prefix, "_gee_fixed_effects.csv")),
            row.names = FALSE)
  write.csv(predictions, file.path(results_dir, paste0(out_prefix, "_gee_predicted_trajectories.csv")),
            row.names = FALSE)
  write.csv(contrast, file.path(results_dir, paste0(out_prefix, "_gee_group_difference.csv")),
            row.names = FALSE)
  write.csv(derivative, file.path(results_dir, paste0(out_prefix, "_gee_derivative.csv")),
            row.names = FALSE)
  write.csv(theta, file.path(results_dir, paste0(out_prefix, "_gee_theta.csv")), row.names = FALSE)

  metadata <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    cohort = cohort,
    requested_4x_model = "M2_full_spline",
    matched_path = normalizePath(matched_path, mustWork = TRUE),
    matched_md5 = provenance$input_md5,
    spline_metadata_path = spline_info$metadata_path,
    spline_df = spline_info$spline_df,
    knots = spline_info$knots,
    boundary_knots = spline_info$boundary_knots,
    spline_terms = spline_info$terms,
    censoring_model_population = result$weights$model_population,
    numerator_formula = result$weights$numerator_formula,
    denominator_formula = result$weights$denominator_formula,
    outcome_formula = result$gee$outcome_formula,
    numerator_fit_method = result$weights$numerator_fit$method,
    denominator_fit_method = result$weights$denominator_fit$method,
    truncation_bounds = result$weights$truncation_bounds,
    prediction_window = c(-16, 16),
    death_treated_as_censoring = FALSE,
    weighted_comparator = "IPCW-weighted GEE minus unweighted GEE"
  )
  saveRDS(
    list(
      metadata = metadata,
      numerator_model = result$weights$numerator_fit$fit,
      denominator_model = result$weights$denominator_fit$fit,
      gee_unweighted = result$gee$unweighted$fit,
      gee_ipcw = result$gee$weighted$fit
    ),
    file.path(results_dir, paste0(out_prefix, "_models.rds"))
  )

  invisible(list(
    data_out_path = data_out_path,
    results_dir = results_dir,
    out_prefix = out_prefix,
    theta = theta,
    censoring_events = censoring_events,
    weight_summary = retained_weights,
    full_response_summary = result$weights$full_response_summary
  ))
}

run_ipcw_m2_gee <- function(panel_path, matched_path, spline_metadata_path,
                            cohort, data_out_path, results_dir, out_prefix,
                            truncation = c(0.01, 0.99),
                            response_model = NULL) {
  provenance <- ipcw_gee_validate_matching_provenance(matched_path)
  spline_info <- ipcw_gee_read_spline_metadata(spline_metadata_path, cohort)
  if (is.null(response_model)) {
    response_model <- ipcw_gee_fit_full_response_model(panel_path)
  }
  panel <- ipcw_gee_build_assignment_panel(matched_path, cohort, response_model)
  weights <- ipcw_gee_prepare_weights(
    panel, spline_info, response_model, truncation = truncation
  )
  gee <- ipcw_gee_fit_pair(weights$panel, spline_info, prediction_window = c(-16, 16))
  result <- list(weights = weights, gee = gee)
  ipcw_gee_write_outputs(
    result = result,
    out_prefix = out_prefix,
    results_dir = results_dir,
    data_out_path = data_out_path,
    cohort = cohort,
    matched_path = matched_path,
    spline_info = spline_info,
    provenance = provenance
  )
}
