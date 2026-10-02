# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.0.3_GLME_m1_df_aic_comparison.R
# Author:  Nemo Zhou
# Date started:      2026-09-22
# Date last updated: 2026-09-23 (added paired fitted-trajectory outputs)
#
# Purpose:
#   Compare M1 natural-spline specifications with df = 3 and df = 4 using AIC
#   within each active 4.x cohort: the two burden cohorts in 4.1, smoking in
#   4.2, obesity in 4.3, and overall cancer in 4.5. Each pair uses the same
#   engine-defined M1 complete-case rows and covariates, the same correlated
#   participant random intercept and attained-age random slope, and maximum
#   likelihood (ML). Only spline degrees of freedom change. The comparison is
#   separate from the formal df = 3 REML production analyses. In addition to
#   AIC/status/provenance, the runner saves model-based fixed-effect trajectory
#   predictions for visual comparison; these are not CR2 inference or theta.
#
#   Input Gate G4/MD5 provenance and the shared engine's M1 filtering, natural
#   spline basis, scaling, and strict convergence checks are reused. Optimizer
#   retries are allowed, but the random-effects structure is held fixed. The
#   script saves paired AIC results, support-gated predictions, per-fit status,
#   convergence attempts, and run/input provenance under Results/cancer/data.
#
# Run from the project root after matching inputs pass Gate G4:
#   Rscript Code/2_data_analysis/4.0.3_GLME_m1_df_aic_comparison.R
# Do not run this GLME comparison in Codex; Nemo runs the command above.
# Sourcing the file defines the runner but does not fit models.
# =============================================================================

.m1_df_aic_cohorts <- data.frame(
  analysis = c("4.1", "4.1", "4.2", "4.3", "4.5"),
  cohort = c(
    "Low/Moderate Burden Cohort",
    "High Burden Cohort",
    "Smoking-Related Cancer Cohort",
    "Obesity-Related Cancer Cohort",
    "All Cancer Cohort"
  ),
  input_file = c(
    "riskset_matched_analysis_long.rds",
    "riskset_matched_analysis_long.rds",
    "riskset_matched_smoking_long.rds",
    "riskset_matched_obesity_long.rds",
    "riskset_matched_overall_long.rds"
  ),
  stringsAsFactors = FALSE
)

.m1_df_aic_hash_keys <- function(keys) {
  path <- tempfile("m1_df_aic_row_keys_")
  on.exit(unlink(path), add = TRUE)
  saveRDS(as.character(keys), path, version = 3)
  unname(tools::md5sum(path))
}

.m1_df_aic_prepare_sample <- function(matched_long, cohort) {
  needed_cols <- c(
    "Cohort", "Group", "Age_Centered", "index_age_z", "base_race",
    "base_marital", "base_living", "id", "match_set", "role", "cycle",
    "fi_score_nocancer", "post_own_cancer", "age_at_cycle"
  )
  missing_cols <- setdiff(needed_cols, names(matched_long))
  if (length(missing_cols)) {
    stop("Matched dataset is missing required M1 fields: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!"matching_id" %in% names(matched_long)) {
    matched_long$matching_id <- paste(matched_long$Cohort,
                                     matched_long$match_set, sep = "__")
  }
  if (!"trajectory_id" %in% names(matched_long)) {
    matched_long$trajectory_id <- paste(matched_long$Cohort,
      matched_long$match_set, matched_long$id, matched_long$role, sep = "__")
  }
  required_identifiers <- c("id", "Group", "matching_id", "trajectory_id",
                            "cycle", "post_own_cancer")
  if (anyNA(matched_long[required_identifiers])) {
    stop("Missing identifiers, group, cycle, or own-cancer censoring flag.",
         call. = FALSE)
  }
  if (!all(matched_long$Group %in% c("Control", "Cancer Case"))) {
    stop("Unexpected Group values in matched data.", call. = FALSE)
  }
  observed_fi <- matched_long$fi_score_nocancer[!is.na(matched_long$fi_score_nocancer)]
  observed_time <- matched_long$Age_Centered[!is.na(matched_long$Age_Centered)]
  if (any(!is.finite(observed_fi) | observed_fi < 0 | observed_fi > 1) ||
      any(!is.finite(observed_time))) {
    stop("FI must be finite in [0,1] and observed relative time must be finite.",
         call. = FALSE)
  }

  d_base <- dplyr::mutate(
    matched_long,
    Cohort = droplevels(factor(Cohort)),
    Group = factor(Group, levels = c("Control", "Cancer Case")),
    id = factor(id),
    match_set = factor(match_set),
    matching_id = factor(matching_id),
    trajectory_id = factor(trajectory_id),
    base_race = factor(base_race),
    base_marital = factor(base_marital),
    base_living = factor(base_living)
  )
  d_cohort <- dplyr::filter(d_base, Cohort == cohort) |>
    droplevels() |>
    dplyr::filter(!post_own_cancer, !is.na(fi_score_nocancer),
                  !is.na(Age_Centered))
  if (!nrow(d_cohort)) return(NULL)

  d_model <- dplyr::filter(
    d_cohort,
    dplyr::if_all(dplyr::all_of(.primary_covars), ~ !is.na(.x))
  ) |>
    droplevels()
  if (!nrow(d_model) || length(unique(d_model$Group)) != 2L) {
    stop("The cohort has no two-arm M1 complete-case sample.", call. = FALSE)
  }
  if (any(!is.finite(d_model$age_at_cycle))) {
    stop("Finite age_at_cycle is required for the attained-age random slope.",
         call. = FALSE)
  }
  if (anyDuplicated(paste(d_model$trajectory_id, d_model$cycle, sep = "__"))) {
    stop("Duplicated trajectory_id x cycle rows remain after M1 filtering.",
         call. = FALSE)
  }
  d_model <- as.data.frame(d_model)
  row_keys <- paste(as.character(d_model$trajectory_id),
                    as.character(d_model$cycle), sep = "__")
  if (anyDuplicated(row_keys)) {
    stop("M1 row keys are not unique within the cohort.", call. = FALSE)
  }
  rownames(d_model) <- row_keys
  list(data = d_model, row_keys = row_keys,
       row_keys_md5 = .m1_df_aic_hash_keys(row_keys))
}

.m1_df_aic_prediction_window <- function(data, window_yrs = 20L,
                                         min_case_bin = 50L,
                                         min_ctrl_bin = 250L) {
  support <- add_relative_time_bin(data) %>%
    dplyr::group_by(rel_time_bin, .drop = FALSE) %>%
    dplyr::summarize(
      n_case = dplyr::n_distinct(id[Group == "Cancer Case"]),
      n_ctrl = dplyr::n_distinct(id[Group == "Control"]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(support_ok = n_case >= min_case_bin & n_ctrl >= min_ctrl_bin)
  window <- continuous_support_window(support, window_yrs)
  list(support = support, window = window)
}

.m1_df_aic_make_predictions <- function(fit, sample, spline, scaling,
                                        fixed_rhs, prediction_window,
                                        analysis, cohort, spline_df) {
  prediction_times <- seq(prediction_window[[1]], prediction_window[[2]], by = 0.25)
  reference_grid <- .make_ref_grid(
    sample$data,
    covars = .primary_covars,
    Age_Centered = prediction_times,
    Group = c("Control", "Cancer Case")
  ) %>%
    dplyr::mutate(Group = factor(Group, levels = c("Control", "Cancer Case")))
  prediction_basis <- stats::predict(spline$basis, reference_grid$Age_Centered)
  reference_grid <- dplyr::bind_cols(
    reference_grid,
    stats::setNames(as.data.frame(unclass(prediction_basis)), spline$terms)
  ) %>%
    .apply_model_scaling(scaling$parameters)

  beta <- lme4::fixef(fit)
  X <- stats::model.matrix(fixed_rhs, data = reference_grid)[,
    names(beta), drop = FALSE]
  V <- as.matrix(stats::vcov(fit))[names(beta), names(beta), drop = FALSE]
  predicted <- as.vector(X %*% beta)
  standard_error <- sqrt(pmax(0, rowSums((X %*% V) * X)))

  data.frame(
    analysis = analysis,
    Cohort = cohort,
    model_id = paste0("M1_df", spline_df),
    model_label = paste0("M1 natural spline, df = ", spline_df),
    spline_df = spline_df,
    Group = as.character(reference_grid$Group),
    Age_Centered = reference_grid$Age_Centered,
    pred = predicted,
    se = standard_error,
    lwr = predicted - 1.96 * standard_error,
    upr = predicted + 1.96 * standard_error,
    vcov_type = "Model-based covariance from ML fit (not CR2)",
    method = "ML",
    n_obs = nrow(sample$data),
    n_participants = dplyr::n_distinct(sample$data$id),
    n_assignments = dplyr::n_distinct(sample$data$trajectory_id),
    row_keys_md5 = sample$row_keys_md5,
    prediction_window_min = prediction_window[[1]],
    prediction_window_max = prediction_window[[2]],
    stringsAsFactors = FALSE
  )
}

.m1_df_aic_fit_one <- function(sample, cohort, analysis, spline_df) {
  row <- sample$data
  spline <- .make_spline_basis(row$Age_Centered, spline_df = spline_df,
                               prefix = "S")
  for (column in spline$terms) {
    row[[column]] <- spline$basis[, column]
  }
  scaling <- .scale_model_numeric_columns(row, spline$terms)
  row <- as.data.frame(scaling$data)
  rownames(row) <- sample$row_keys
  fixed_rhs <- .fixed_rhs(spline$terms, covars = .primary_covars)
  formula <- stats::as.formula(paste(
    "fi_score_nocancer", paste(deparse(fixed_rhs), collapse = ""),
    "+ (1 + .random_time | id)"
  ))
  model_id <- paste0("M1_df", spline_df)

  fit_info <- tryCatch(
    .fit_lmer_with_ladder(
      formula = formula,
      data = row,
      model_label = paste(analysis, cohort, model_id, sep = " / "),
      numerical_retries = TRUE,
      simplify = FALSE,
      random_clock = "attained_age",
      random_time_scale = 4,
      reml = FALSE
    ),
    error = function(e) e
  )
  attempts <- if (inherits(fit_info, "error")) {
    attr(fit_info, "attempts")
  } else {
    fit_info$attempts
  }
  if (!is.null(attempts) && nrow(attempts)) {
    attempts <- dplyr::mutate(
      attempts,
      analysis = analysis,
      Cohort = cohort,
      model_id = model_id,
      spline_df = spline_df,
      .before = 1
    )
  }

  message_text <- NA_character_
  status <- "fit"
  aic <- log_likelihood <- NA_real_
  nobs_fit <- NA_integer_
  fit_row_keys_match <- NA
  selected_rung <- selected_optimizer <- NA_character_
  ml_verified <- FALSE
  prediction_status <- "not_attempted"
  prediction_message <- NA_character_
  prediction_window <- c(NA_real_, NA_real_)
  predictions <- data.frame()
  if (inherits(fit_info, "error")) {
    status <- "failed_strict_convergence"
    message_text <- conditionMessage(fit_info)
  } else {
    fit <- fit_info$fit
    fit_keys <- rownames(stats::model.frame(fit))
    fit_row_keys_match <- identical(fit_keys, sample$row_keys)
    ml_verified <- !lme4::isREML(fit)
    nobs_fit <- as.integer(stats::nobs(fit))
    selected_attempt <- attempts[attempts$accepted %in% TRUE, , drop = FALSE]
    if (nrow(selected_attempt)) {
      selected_rung <- as.character(tail(selected_attempt$rung, 1L))
      selected_optimizer <- as.character(tail(selected_attempt$optimizer, 1L))
    }
    if (!fit_row_keys_match || nobs_fit != nrow(sample$data) || !ml_verified) {
      status <- "failed_postfit_validation"
      message_text <- paste(
        if (!fit_row_keys_match) "Fitted model row keys differ from the common M1 sample.",
        if (nobs_fit != nrow(sample$data)) "Fitted model observation count differs from the common M1 sample.",
        if (!ml_verified) "Fitted model log likelihood is not verified as ML.",
        collapse = " "
      )
    } else {
      aic <- as.numeric(stats::AIC(fit))
      log_likelihood <- as.numeric(stats::logLik(fit))
      if (!is.finite(aic) || !is.finite(log_likelihood)) {
        status <- "failed_postfit_validation"
        message_text <- "Fitted ML likelihood or AIC is not finite."
        aic <- log_likelihood <- NA_real_
      }
    }
    if (status == "fit") {
      support_result <- tryCatch(
        .m1_df_aic_prediction_window(sample$data),
        error = function(e) e
      )
      if (inherits(support_result, "error")) {
        prediction_status <- "failed_support_check"
        prediction_message <- conditionMessage(support_result)
      } else {
        prediction_window <- support_result$window
        if (any(!is.finite(prediction_window)) ||
            prediction_window[[1]] >= prediction_window[[2]]) {
          prediction_status <- "no_contiguous_supported_window"
          prediction_message <- "The reference bin or a contiguous two-arm-supported interval was unavailable."
        } else {
          prediction_result <- tryCatch(
            .m1_df_aic_make_predictions(
              fit, sample, spline, scaling, fixed_rhs, prediction_window,
              analysis, cohort, spline_df
            ),
            error = function(e) e
          )
          if (inherits(prediction_result, "error")) {
            prediction_status <- "failed_prediction"
            prediction_message <- conditionMessage(prediction_result)
          } else if (!nrow(prediction_result) ||
                     any(!is.finite(prediction_result$pred)) ||
                     any(!is.finite(prediction_result$se))) {
            prediction_status <- "failed_prediction_validation"
            prediction_message <- "Predictions or standard errors were empty or non-finite."
          } else {
            predictions <- prediction_result
            prediction_status <- "saved"
          }
        }
      }
    }
  }

  status_row <- data.frame(
    analysis = analysis,
    Cohort = cohort,
    model_id = model_id,
    spline_df = spline_df,
    status = status,
    method = "ML",
    reml = FALSE,
    ml_verified = ml_verified,
    ml_verification_basis = if (ml_verified) "lme4::isREML(fit) == FALSE" else "not verified",
    aic_source = "stats::AIC(fit)",
    n_obs = nrow(sample$data),
    n_obs_fit = nobs_fit,
    n_participants = dplyr::n_distinct(sample$data$id),
    n_assignments = dplyr::n_distinct(sample$data$trajectory_id),
    row_keys_md5 = sample$row_keys_md5,
    fit_row_keys_match = fit_row_keys_match,
    fixed_covariates = paste(.primary_covars, collapse = " + "),
    random_effects = "(1 + .random_time | id)",
    random_clock = "(age_at_cycle - 60) / 4",
    spline_knots = paste(signif(spline$knots, 10), collapse = ";"),
    spline_boundary_knots = paste(signif(spline$boundary_knots, 10), collapse = ";"),
    spline_scaling = paste(
      paste0(scaling$parameters$column, ":center=",
             signif(scaling$parameters$center, 8), ":scale=",
             signif(scaling$parameters$scale, 8)),
      collapse = " | "
    ),
    convergence_rung = selected_rung,
    optimizer = selected_optimizer,
    convergence_attempts = if (is.null(attempts)) 0L else nrow(attempts),
    prediction_status = prediction_status,
    prediction_window_min = prediction_window[[1]],
    prediction_window_max = prediction_window[[2]],
    n_predictions = nrow(predictions),
    prediction_message = prediction_message,
    AIC = aic,
    logLik = log_likelihood,
    formula = paste(deparse(formula), collapse = ""),
    message = message_text,
    stringsAsFactors = FALSE
  )
  list(status = status_row, attempts = attempts, predictions = predictions)
}

run_m1_df_aic_comparison <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  engine_path <- file.path(project_dir, "Code", "2_data_analysis",
                           "4.0_GLME_spline_functions.R")
  if (!file.exists(engine_path)) {
    stop("Run this script from the Frailty HPFS project root.", call. = FALSE)
  }
  source(engine_path, local = .GlobalEnv)

  data_dir <- file.path(project_dir, "Data")
  results_dir <- file.path(project_dir, "Results", "cancer", "data")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  input_provenance <- list()
  input_data <- list()
  for (input_file in unique(.m1_df_aic_cohorts$input_file)) {
    input_path <- file.path(data_dir, input_file)
    provenance <- validate_matching_provenance(input_path)
    input_data[[input_file]] <- readRDS(input_path)
    input_provenance[[input_file]] <- list(
      input_path = normalizePath(input_path, mustWork = TRUE),
      input_md5 = provenance$input_md5,
      gate_path = normalizePath(provenance$gate_path, mustWork = TRUE),
      gate_md5 = unname(tools::md5sum(provenance$gate_path)),
      matching_run_path = normalizePath(provenance$run_path, mustWork = TRUE),
      matching_run_md5 = unname(tools::md5sum(provenance$run_path)),
      matching_output_stem = provenance$output_stem
    )
  }

  fit_status <- list()
  convergence_attempts <- list()
  trajectory_predictions <- list()
  comparison_rows <- list()
  row_key_audit <- list()
  for (i in seq_len(nrow(.m1_df_aic_cohorts))) {
    spec <- .m1_df_aic_cohorts[i, , drop = FALSE]
    matched_long <- input_data[[spec$input_file]]
    sample_result <- tryCatch(
      .m1_df_aic_prepare_sample(matched_long, spec$cohort),
      error = function(e) e
    )
    if (inherits(sample_result, "error")) {
      for (spline_df in c(3L, 4L)) {
        fit_status[[length(fit_status) + 1L]] <- data.frame(
          analysis = spec$analysis, Cohort = spec$cohort,
          model_id = paste0("M1_df", spline_df), spline_df = spline_df,
          status = "failed_sample_validation", method = "ML", reml = FALSE,
          ml_verified = FALSE, ml_verification_basis = "not verified",
          aic_source = "not computed", n_obs = NA_integer_, n_obs_fit = NA_integer_,
          n_participants = NA_integer_, n_assignments = NA_integer_,
          row_keys_md5 = NA_character_, fit_row_keys_match = NA,
          fixed_covariates = paste(.primary_covars, collapse = " + "),
          random_effects = "(1 + .random_time | id)",
          random_clock = "(age_at_cycle - 60) / 4",
          spline_knots = NA_character_, spline_boundary_knots = NA_character_,
          spline_scaling = NA_character_, convergence_rung = NA_character_,
          optimizer = NA_character_, convergence_attempts = 0L,
          AIC = NA_real_, logLik = NA_real_, formula = NA_character_,
          message = conditionMessage(sample_result), stringsAsFactors = FALSE
        )
      }
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        analysis = spec$analysis, Cohort = spec$cohort,
        status = "failed_sample_validation", AIC_df3 = NA_real_,
        AIC_df4 = NA_real_, delta_AIC_df4_minus_df3 = NA_real_,
        n_obs = NA_integer_, n_participants = NA_integer_,
        n_assignments = NA_integer_, row_keys_md5 = NA_character_,
        row_keys_identical = NA, fixed_covariates = paste(.primary_covars,
          collapse = " + "), random_effects = "(1 + .random_time | id)",
        method = "ML", formula_comparison = "M1 fixed effects with natural spline df 3 versus 4",
        aic_source = "not computed",
        message = conditionMessage(sample_result), stringsAsFactors = FALSE
      )
      next
    }
    if (is.null(sample_result)) {
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        analysis = spec$analysis, Cohort = spec$cohort,
        status = "skipped_cohort_not_found", AIC_df3 = NA_real_,
        AIC_df4 = NA_real_, delta_AIC_df4_minus_df3 = NA_real_,
        n_obs = NA_integer_, n_participants = NA_integer_,
        n_assignments = NA_integer_, row_keys_md5 = NA_character_,
        row_keys_identical = NA, fixed_covariates = paste(.primary_covars,
          collapse = " + "), random_effects = "(1 + .random_time | id)",
        method = "ML", formula_comparison = "M1 fixed effects with natural spline df 3 versus 4",
        aic_source = "not computed",
        message = "Expected cohort label was not found in the matched input.",
        stringsAsFactors = FALSE
      )
      next
    }
    row_key_audit[[paste(spec$analysis, spec$cohort, sep = "__")]] <- list(
      analysis = spec$analysis, cohort = spec$cohort,
      n_row_keys = length(sample_result$row_keys),
      row_keys_md5 = sample_result$row_keys_md5,
      common_sample_used_for_both_df = TRUE,
      fitted_row_keys_verified = FALSE
    )

    pair_fits <- lapply(c(3L, 4L), function(spline_df) {
      .m1_df_aic_fit_one(sample_result, spec$cohort, spec$analysis,
                         spline_df)
    })
    for (fit_result in pair_fits) {
      fit_status[[length(fit_status) + 1L]] <- fit_result$status
      if (!is.null(fit_result$attempts) && nrow(fit_result$attempts)) {
        convergence_attempts[[length(convergence_attempts) + 1L]] <- fit_result$attempts
      }
    }
    row_key_audit[[paste(spec$analysis, spec$cohort, sep = "__")]]$fitted_row_keys_verified <-
      all(vapply(pair_fits, function(x)
        isTRUE(x$status$fit_row_keys_match[[1]]), logical(1)))
    status_pair <- vapply(pair_fits, function(x) x$status$status[[1]], character(1))
    aic_pair <- vapply(pair_fits, function(x) x$status$AIC[[1]], numeric(1))
    row_keys_match <- all(vapply(pair_fits, function(x)
      isTRUE(x$status$fit_row_keys_match[[1]]), logical(1)))
    pair_specification_validated <-
      identical(vapply(pair_fits, function(x) x$status$fixed_covariates[[1]],
                       character(1)), rep(paste(.primary_covars, collapse = " + "), 2L)) &&
      identical(vapply(pair_fits, function(x) x$status$random_effects[[1]],
                       character(1)), rep("(1 + .random_time | id)", 2L)) &&
      identical(vapply(pair_fits, function(x) x$status$method[[1]], character(1)),
                c("ML", "ML")) &&
      identical(vapply(pair_fits, function(x) x$status$reml[[1]], logical(1)),
                c(FALSE, FALSE)) &&
      identical(vapply(pair_fits, function(x) x$status$row_keys_md5[[1]],
                       character(1)), rep(sample_result$row_keys_md5, 2L)) &&
      all(vapply(pair_fits, function(x) x$status$n_obs[[1]], integer(1)) ==
            nrow(sample_result$data))
    pair_ok <- pair_specification_validated && row_keys_match &&
      all(status_pair == "fit") && all(is.finite(aic_pair))
    prediction_pair_ok <- pair_ok && all(vapply(
      pair_fits,
      function(x) identical(x$status$prediction_status[[1]], "saved"),
      logical(1)
    ))
    if (prediction_pair_ok) {
      trajectory_predictions[[length(trajectory_predictions) + 1L]] <-
        dplyr::bind_rows(lapply(pair_fits, `[[`, "predictions"))
    }
    pair_messages <- character()
    if (!pair_specification_validated) {
      pair_messages <- c(pair_messages,
        "The paired ML, row-count, covariate, or random-effects specification check failed.")
    }
    if (!row_keys_match) {
      pair_messages <- c(pair_messages,
        "At least one fitted model did not retain the common M1 row keys.")
    }
    failed_indices <- which(status_pair != "fit")
    if (length(failed_indices)) {
      pair_messages <- c(pair_messages, vapply(failed_indices, function(j) {
        msg <- pair_fits[[j]]$status$message[[1]]
        paste0("df ", c(3L, 4L)[[j]], " ", status_pair[[j]],
               if (!is.na(msg) && nzchar(msg)) paste0(": ", msg) else "")
      }, character(1)))
    }
    if (any(status_pair == "fit") && any(!is.finite(aic_pair))) {
      pair_messages <- c(pair_messages, "At least one fitted model has no finite AIC.")
    }
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      analysis = spec$analysis,
      Cohort = spec$cohort,
      status = if (pair_ok) "fit" else "failed_pair",
      AIC_df3 = aic_pair[[1]],
      AIC_df4 = aic_pair[[2]],
      delta_AIC_df4_minus_df3 = if (pair_ok) aic_pair[[2]] - aic_pair[[1]] else NA_real_,
      n_obs = nrow(sample_result$data),
      n_participants = dplyr::n_distinct(sample_result$data$id),
      n_assignments = dplyr::n_distinct(sample_result$data$trajectory_id),
      row_keys_md5 = sample_result$row_keys_md5,
      row_keys_identical = row_keys_match,
      pair_specification_validated = pair_specification_validated,
      fixed_covariates = paste(.primary_covars, collapse = " + "),
      random_effects = "(1 + .random_time | id)",
      method = "ML",
      formula_comparison = "M1 fixed effects with natural spline df 3 versus 4",
      aic_source = if (pair_ok) "stats::AIC(fit)" else "not computed",
      prediction_status = if (prediction_pair_ok) "saved" else "not_saved_for_complete_pair",
      message = if (pair_ok) NA_character_ else paste(pair_messages, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }

  comparison <- dplyr::bind_rows(comparison_rows)
  status_table <- dplyr::bind_rows(fit_status)
  attempts_table <- dplyr::bind_rows(convergence_attempts)
  prediction_table <- dplyr::bind_rows(trajectory_predictions)
  if (!nrow(prediction_table)) {
    prediction_table <- data.frame(
      analysis = character(), Cohort = character(), model_id = character(),
      model_label = character(), spline_df = integer(), Group = character(),
      Age_Centered = numeric(), pred = numeric(), se = numeric(), lwr = numeric(),
      upr = numeric(), vcov_type = character(), method = character(),
      n_obs = integer(), n_participants = integer(), n_assignments = integer(),
      row_keys_md5 = character(), prediction_window_min = numeric(),
      prediction_window_max = numeric(), stringsAsFactors = FALSE
    )
  }
  if (!nrow(attempts_table)) {
    attempts_table <- data.frame(
      note = "No optimizer attempts were recorded; inspect fit status for pre-fit failures.",
      stringsAsFactors = FALSE
    )
  }
  if (any(comparison$status == "fit" &
          abs(comparison$delta_AIC_df4_minus_df3 -
              (comparison$AIC_df4 - comparison$AIC_df3)) > 1e-10,
          na.rm = TRUE)) {
    stop("Internal check failed: delta AIC does not equal AIC(df 4) minus AIC(df 3).",
         call. = FALSE)
  }

  runner_path <- file.path(project_dir, "Code", "2_data_analysis",
                           "4.0.3_GLME_m1_df_aic_comparison.R")
  metadata <- list(
    schema_version = 2L,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    runner = normalizePath(runner_path, mustWork = TRUE),
    runner_md5 = unname(tools::md5sum(runner_path)),
    shared_engine = normalizePath(engine_path, mustWork = TRUE),
    shared_engine_md5 = unname(tools::md5sum(engine_path)),
    comparison = paste(
      "Within-cohort M1 natural spline df 3 versus df 4 by ML AIC",
      "with support-gated fitted trajectory predictions"
    ),
    estimator = "Maximum likelihood (REML = FALSE)",
    delta_definition = "AIC(df 4) - AIC(df 3); negative values favor df 4",
    random_effects = "(1 + (age_at_cycle - 60)/4 | id), correlated intercept and slope",
    fixed_covariates = .primary_covars,
    sample_rule = "Shared engine M1 rows: not post_own_cancer, observed fi_score_nocancer and Age_Centered, complete on all primary M1 covariates",
    spline_basis = "Shared engine natural-spline basis construction and per-fit numerical centering/scaling; only df differs",
    prediction_specification = paste(
      "Fixed-effect predictions at the shared M1 reference profile on a 0.25-year grid",
      "within the contiguous +/-20-year interval supported by >=50 distinct case",
      "and >=250 distinct control participants per 2-year bin; pointwise 95% intervals",
      "use conventional model-based covariance from each ML fit, not CR2"
    ),
    optimizer_policy = "bobyqa then nloptwrap retries; random-effects simplification disabled",
    cohorts = .m1_df_aic_cohorts,
    inputs = input_provenance,
    row_key_audit = row_key_audit,
    output_files = c(
      comparison = "4.0.3_m1_spline_df_aic_comparison.csv",
      status = "4.0.3_m1_spline_df_fit_status.csv",
      convergence_attempts = "4.0.3_m1_spline_df_convergence_attempts.csv",
      predictions = "4.0.3_m1_spline_df_predicted_trajectories.csv",
      metadata = "4.0.3_m1_spline_df_run_metadata.rds"
    )
  )
  write.csv(comparison,
            file.path(results_dir, "4.0.3_m1_spline_df_aic_comparison.csv"),
            row.names = FALSE, na = "")
  write.csv(status_table,
            file.path(results_dir, "4.0.3_m1_spline_df_fit_status.csv"),
            row.names = FALSE, na = "")
  write.csv(attempts_table,
            file.path(results_dir, "4.0.3_m1_spline_df_convergence_attempts.csv"),
            row.names = FALSE, na = "")
  write.csv(prediction_table,
            file.path(results_dir, "4.0.3_m1_spline_df_predicted_trajectories.csv"),
            row.names = FALSE, na = "")
  saveRDS(metadata,
          file.path(results_dir, "4.0.3_m1_spline_df_run_metadata.rds"))
  message("Saved paired M1 df AIC comparison and provenance under ", results_dir)
  invisible(list(comparison = comparison, fit_status = status_table,
                 convergence_attempts = attempts_table,
                 predictions = prediction_table, metadata = metadata))
}

if (sys.nframe() == 0L) {
  run_m1_df_aic_comparison()
}
