# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script: 4.0.1_GLME_sensitivity_analyses.R
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
# Purpose: Run explicitly selected overall M1 sensitivity analyses independently
# of the formal 4.0/4.5 model run. S3-S6 reuse a provenance-checked M1 cache;
# other selected analyses fit their own models. S7 is retired; no exact-cycle data are read.
# S11 retains its explicitly flagged weighted-inference limitation.
# S12 (two-visit support) is preserved as an optional existing analysis.
# No models run on source(), or on Rscript invocation without arguments.
# Nemo runs modeling commands from the project root, for example:
# Rscript Code/2_data_analysis/4.0.1_GLME_sensitivity_analyses.R --only=S1,S2
# Rscript Code/2_data_analysis/4.0.1_GLME_sensitivity_analyses.R --only=S3,S4,S5,S6
# --all explicitly selects S1-S6 and S8-S11; select --only=S12 for the additional analysis.
# Outputs: Results/cancer/data/4.5_sensitivities/<selection>/
# Each selection has separate summaries and provenance; other selections remain.
# =============================================================================
source("Code/2_data_analysis/4.0_GLME_spline_functions.R")

# Overall-cancer sensitivity grid ---------------------------------------------
# This deliberately fits the M1 primary-adjusted spline only. M0--M3 are
# already produced by run_spline_analysis(); each sensitivity
# changes one feature of M1 and writes harmonized, row-bindable outputs.
.filter_two_visit_sensitivity <- function(d) {
  support <- d %>%
    filter(!post_own_cancer, !is.na(fi_score_nocancer), !is.na(Age_Centered)) %>%
    distinct(trajectory_id, cycle) %>%
    count(trajectory_id, name = "n_distinct_fi_dates")
  keep_trajectory <- support$trajectory_id[support$n_distinct_fi_dates >= 2]
  d2 <- d %>% filter(trajectory_id %in% keep_trajectory)
  valid_sets <- d2 %>%
    distinct(matching_id, trajectory_id, Group) %>%
    group_by(matching_id) %>%
    summarize(n_case = sum(Group == "Cancer Case"),
              n_control = sum(Group == "Control"), .groups = "drop") %>%
    filter(n_case == 1, n_control >= 1)
  out <- d2 %>% filter(matching_id %in% valid_sets$matching_id)
  retained_integrity <- out %>%
    distinct(matching_id, trajectory_id, Group) %>%
    group_by(matching_id) %>%
    summarize(n_case = sum(Group == "Cancer Case"),
              n_control = sum(Group == "Control"), .groups = "drop")
  if (any(retained_integrity$n_case != 1 | retained_integrity$n_control < 1)) {
    stop("S12 matched-set integrity failed after the two-visit restriction.")
  }
  retained_visits <- out %>%
    distinct(trajectory_id, cycle) %>%
    count(trajectory_id, name = "n_distinct_fi_dates")
  if (any(retained_visits$n_distinct_fi_dates < 2)) {
    stop("S12 retained a trajectory with fewer than two distinct analytic FI dates.")
  }
  attr(out, "visit_support") <- support
  attr(out, "set_support") <- retained_integrity
  out
}

.fit_one_spline_sensitivity <- function(d, spec, cohort_label = "All Cancer Cohort", diagnostic_path = NULL) {
  required <- c("Group", "id", "matching_id", "trajectory_id", "cycle",
                "Age_Centered", "index_age_z", "base_race", "base_marital",
                "base_living", "post_own_cancer", spec$outcome)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols)) stop("missing columns: ", paste(missing_cols, collapse = ", "))

  d <- d %>%
    mutate(Group = factor(Group, levels = c("Control", "Cancer Case")),
           id = factor(id), matching_id = factor(matching_id),
           trajectory_id = factor(trajectory_id), cycle = factor(cycle),
           base_race = factor(base_race), base_marital = factor(base_marital),
           base_living = factor(base_living)) %>%
    filter(!post_own_cancer, !is.na(.data[[spec$outcome]]), !is.na(Age_Centered),
           if_all(all_of(.primary_covars), ~ !is.na(.x))) %>%
    droplevels()
  visit_support <- NULL
  set_support <- NULL
  if (identical(spec$id, "S12")) {
    d <- .filter_two_visit_sensitivity(d)
    visit_support <- attr(d, "visit_support")
    set_support <- attr(d, "set_support")
  }
  if (nrow(d) == 0 || n_distinct(d$Group) < 2) stop("insufficient two-arm support")

  min_case_bin <- if (is.null(spec$min_case_bin)) 50 else spec$min_case_bin
  min_ctrl_bin <- if (is.null(spec$min_ctrl_bin)) 250 else spec$min_ctrl_bin
  support <- add_relative_time_bin(d) %>%
    group_by(rel_time_bin, .drop = FALSE) %>%
    summarize(n_case = n_distinct(id[Group == "Cancer Case"]),
              n_ctrl = n_distinct(id[Group == "Control"]), .groups = "drop") %>%
    mutate(support_ok = n_case >= min_case_bin & n_ctrl >= min_ctrl_bin,
           spec = spec$id)
  supported_window <- continuous_support_window(support, spec$prediction_window)
  if (anyNA(supported_window) || supported_window[[1]] > -spec$prediction_window ||
      supported_window[[2]] < spec$prediction_window) {
    stop("requested +/-", spec$prediction_window,
         " prediction window lacks continuous two-arm support (requires >=",
         min_case_bin, " cases and >=", min_ctrl_bin, " controls per bin)")
  }

  y_name <- spec$outcome
  if (identical(spec$id, "S10")) {
    if (!("n_answered_nocancer" %in% names(d))) {
      stop("n_answered_nocancer is required for S10")
    }
    if (any(d$n_answered_nocancer <= 1, na.rm = TRUE)) {
      stop("S10 requires n_answered_nocancer > 1")
    }
    d <- d %>% mutate(
      .outcome_s10 = qlogis((.data[[y_name]] * (n_answered_nocancer - 1) + 0.5) /
                              n_answered_nocancer)
    )
    y_name <- ".outcome_s10"
  }

  B <- ns(d$Age_Centered, df = spec$spline_df)
  knots <- attr(B, "knots")
  boundary <- attr(B, "Boundary.knots")
  sterms <- paste0("S", seq_len(ncol(B)))
  d <- bind_cols(d, setNames(as.data.frame(B), sterms))
  covars <- .primary_covars
  if (identical(spec$id, "S8")) covars <- c(covars, "cycle")
  rhs <- .fixed_rhs(sterms, covars)
  form <- as.formula(paste(y_name, paste(deparse(rhs), collapse = ""),
                           "+ (1 + .random_time | id)"))
  fit_weights <- NULL
  if (identical(spec$id, "S11")) {
    assignment_weights <- d %>%
      distinct(matching_id, trajectory_id, Group) %>%
      group_by(matching_id) %>%
      mutate(n_controls = sum(Group == "Control"),
             match_weight = if_else(Group == "Cancer Case", 1, 1 / n_controls)) %>%
      ungroup() %>% select(matching_id, trajectory_id, match_weight)
    d <- d %>% left_join(assignment_weights, by = c("matching_id", "trajectory_id"))
    fit_weights <- d$match_weight
  }
  fit_info <- .fit_lmer_with_ladder(form, d, paste("sensitivity", spec$id),
                                    weights = fit_weights)
  fit <- fit_info$fit
  beta <- fixef(fit)

  reference_profile <- .make_ref_grid(d, covars = covars, .reference = 1)
  make_grid <- function(t, group) {
    g <- .make_ref_grid(reference_profile, covars = covars, Age_Centered = t, Group = group) %>%
      mutate(Group = factor(Group, levels = c("Control", "Cancer Case")))
    bg <- ns(g$Age_Centered, knots = knots, Boundary.knots = boundary)
    bind_cols(g, setNames(as.data.frame(bg), sterms))
  }
  design <- function(t, group) {
    model.matrix(rhs, make_grid(t, group))[, names(beta), drop = FALSE]
  }
  diff_design <- function(t) design(t, "Cancer Case") - design(t, "Control")
  deriv_design <- function(t, h = 1e-4) {
    (diff_design(t + h) - diff_design(t - h)) / (2 * h)
  }
  constraints <- .inference_constraints(names(beta), sterms, diff_design, spec$theta_window)
  Vobj <- get_primary_vcov(fit, model.frame(fit)$id, constraints = constraints,
                            diagnostic_path = diagnostic_path)
  times <- seq(-spec$prediction_window, spec$prediction_window, by = .25)
  Xc <- design(times, "Cancer Case")
  Xr <- design(times, "Control")
  Xd <- Xc - Xr
  mu_c <- as.vector(Xc %*% beta)
  mu_r <- as.vector(Xr %*% beta)
  est <- as.vector(Xd %*% beta)
  se <- sqrt(rowSums((Xd %*% Vobj$V) * Xd))
  response_diff <- rep(NA_real_, length(times))
  response_se <- rep(NA_real_, length(times))
  response_lwr <- rep(NA_real_, length(times))
  response_upr <- rep(NA_real_, length(times))
  if (identical(spec$id, "S10")) {
    response_diff <- plogis(mu_c) - plogis(mu_r)
    response_gradient <- Xc * (plogis(mu_c) * (1 - plogis(mu_c))) -
      Xr * (plogis(mu_r) * (1 - plogis(mu_r)))
    response_se <- sqrt(rowSums((response_gradient %*% Vobj$V) * response_gradient))
    response_lwr <- response_diff - 1.96 * response_se
    response_upr <- response_diff + 1.96 * response_se
  }
  curve <- data.frame(
    spec = spec$id, Age_Centered = times, estimate = est, se = se,
    lwr = est - 1.96 * se, upr = est + 1.96 * se,
    case_prediction = mu_c, control_prediction = mu_r,
    back_transformed_difference = response_diff,
    back_transformed_se = response_se,
    back_transformed_lwr = response_lwr,
    back_transformed_upr = response_upr,
    vcov_type = Vobj$type
  )

  w <- spec$theta_window
  exact <- .exact_slope_contrasts(diff_design, w)
  c_post <- exact$post; c_pre <- exact$pre; c_theta <- exact$theta
  theta <- as.numeric(c_theta %*% beta)
  theta_closed <- (as.numeric(diff_design(w) %*% beta) -
                     as.numeric(diff_design(0) %*% beta)) / w -
    (as.numeric(diff_design(0) %*% beta) -
       as.numeric(diff_design(-w) %*% beta)) / w
  theta_difference <- theta - theta_closed
  if (!is.finite(theta_difference) || abs(theta_difference) > 1e-6) {
    stop(spec$id, " exact contrast and endpoint theta differ by more than 1e-6")
  }
  theta_se <- sqrt(as.numeric(t(c_theta) %*% Vobj$V %*% c_theta))
  theta_wald <- wald_with_vcov(fit, Vobj, matrix(c_theta, nrow = 1), "theta")
  theta_df <- theta_wald$df_denom[[1]]
  crit <- if (is.na(theta_df)) NA_real_ else qt(.975, theta_df)
  theta_out <- data.frame(
    spec = spec$id, theta = theta, se = theta_se,
    lwr = theta - crit * theta_se, upr = theta + crit * theta_se,
    p_value = theta_wald$p_value[[1]], df = theta_df,
    pre_window = paste0("[-", w, ",0)"), post_window = paste0("(0,+", w, "]"),
    theta_closed_form = theta_closed,
    theta_closed_form_difference = theta_difference,
    infer_method = theta_wald$infer_method[[1]],
    inference_status = Vobj$inference_status,
    transformed_scale = identical(spec$id, "S10"), vcov_type = Vobj$type,
    rung = fit_info$rung
  )
  status <- data.frame(
    spec = spec$id, status = "fit", n_obs = nrow(d), n_id = n_distinct(d$id),
    n_trajectory_id = n_distinct(d$trajectory_id),
    n_matching_id = n_distinct(d$matching_id),
    fitting_min_time = min(d$Age_Centered), fitting_max_time = max(d$Age_Centered),
    prediction_window = spec$prediction_window, spline_df = spec$spline_df,
    outcome = spec$outcome, rung = fit_info$rung, singular = fit_info$singular,
    convergence = fit_info$converged, optimizer_code = fit_info$optimizer_code,
    scaled_gradient = fit_info$scaled_gradient, vcov_type = Vobj$type,
    inference_status = Vobj$inference_status, message = NA_character_
  )
  list(curve = curve, theta = theta_out, status = status,
       convergence = fit_info$attempts %>% mutate(spec = spec$id, .before = 1),
       variance = .variance_components(fit, cohort_label, spec$id),
       support = support, visit_support = visit_support, set_support = set_support)
}

# S3--S6 change only the display or contrast window: use the exact M1 fit.
.summarize_cached_sensitivity <- function(cache, spec) {
  required_window <- max(spec$prediction_window, spec$theta_window)
  if (cache$window[1] > -required_window || cache$window[2] < required_window)
    stop("Requested sensitivity window lacks continuous M1 support")
  beta <- fixef(cache$fit); Vobj <- cache$Vobj
  times <- seq(-spec$prediction_window, spec$prediction_window, by = .25)
  Xc <- cache$grids$design(times, "Cancer Case")
  Xr <- cache$grids$design(times, "Control")
  Xd <- Xc - Xr
  est <- as.vector(Xd %*% beta)
  se <- sqrt(rowSums((Xd %*% Vobj$V) * Xd))
  curve <- data.frame(spec = spec$id, Age_Centered = times, estimate = est, se = se,
    lwr = est - 1.96 * se, upr = est + 1.96 * se,
    case_prediction = as.vector(Xc %*% beta), control_prediction = as.vector(Xr %*% beta),
    back_transformed_difference = NA_real_, back_transformed_se = NA_real_,
    back_transformed_lwr = NA_real_, back_transformed_upr = NA_real_, vcov_type = Vobj$type)
  C <- .exact_slope_contrasts(cache$grids$difference, spec$theta_window)$theta
  theta <- as.numeric(C %*% beta)
  theta_se <- sqrt(as.numeric(C %*% Vobj$V %*% C))
  wt <- wald_with_vcov(cache$fit, Vobj, matrix(C, nrow = 1), paste0(spec$id, " theta"))
  crit <- if (is.na(wt$df_denom)) NA_real_ else qt(.975, wt$df_denom)
  theta_out <- data.frame(spec = spec$id, theta = theta, se = theta_se,
    lwr = theta - crit * theta_se, upr = theta + crit * theta_se,
    p_value = wt$p_value, df = wt$df_denom,
    pre_window = paste0("[-", spec$theta_window, ",0)"),
    post_window = paste0("(0,+", spec$theta_window, "]"),
    theta_closed_form = theta, theta_closed_form_difference = 0,
    transformed_scale = FALSE, vcov_type = Vobj$type, rung = cache$rung)
  status <- data.frame(spec = spec$id, status = "fit", reused_primary_fit = TRUE,
    n_obs = nobs(cache$fit), n_id = length(unique(model.frame(cache$fit)$id)),
    prediction_window = spec$prediction_window, spline_df = 3,
    outcome = spec$outcome, rung = cache$rung, convergence = TRUE,
    singular = isSingular(cache$fit), vcov_type = Vobj$type,
    inference_status = Vobj$inference_status, message = "Reused M1; no model refit")
  list(curve = curve, theta = theta_out, status = status,
    convergence = NULL, variance = NULL,
    support = mutate(cache$support, spec = spec$id), visit_support = NULL, set_support = NULL)
}

run_overall_glme_sensitivities <- function(matched_path,
                                           sensitivity_dir, primary_result = NULL, selected) {
  if (!missing(selected) && "S7" %in% selected)
    stop("S7 was retired: exact-questionnaire-cycle eligibility is no longer an active sensitivity.")
  if (missing(selected) || !length(selected) || any(!selected %in% paste0("S", c(1:6, 8:12))))
    stop("Select sensitivity IDs explicitly (S1-S6 or S8-S12)")
  selected <- unique(selected)
  if (any(selected %in% c("S3", "S4", "S5", "S6")) && is.null(primary_result))
    stop("S3-S6 require a verified saved M1 cache; run 4.5 first")
  if (!dir.exists(sensitivity_dir)) dir.create(sensitivity_dir, recursive = TRUE)
  base_provenance <- validate_matching_provenance(matched_path)
  base <- readRDS(matched_path)
  if (!("matching_id" %in% names(base))) {
    base <- base %>% mutate(matching_id = paste(Cohort, match_set, sep = "__"))
  }
  if (!("trajectory_id" %in% names(base))) {
    base <- base %>% mutate(trajectory_id = paste(Cohort, match_set, id, role, sep = "__"))
  }
  specs <- list(
    list(id = "S1", spline_df = 2, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S2", spline_df = 4, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S3", spline_df = 3, prediction_window = 8, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S4", spline_df = 3, prediction_window = 12, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S5", spline_df = 3, prediction_window = 20, theta_window = 4,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S6", spline_df = 3, prediction_window = 20, theta_window = 12,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S8", spline_df = 3, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S9", spline_df = 3, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer_nocarry", data = base),
    list(id = "S10", spline_df = 3, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S11", spline_df = 3, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base),
    list(id = "S12", spline_df = 3, prediction_window = 20, theta_window = 8,
         outcome = "fi_score_nocancer", data = base)
  )
  specs <- specs[vapply(specs, function(x) x$id %in% selected, logical(1))]

  ans <- lapply(specs, function(spec) {
    tryCatch({
      if (spec$id %in% c("S3", "S4", "S5", "S6") && !is.null(primary_result)) {
        if (!identical(primary_result$matching_provenance$input_md5, base_provenance$input_md5) ||
            !identical(primary_result$matching_provenance$assignment_md5, base_provenance$assignment_md5))
          stop("Cached primary fit uses a different matched dataset")
        cache <- primary_result$sp_contexts[["All Cancer Cohort"]]
        if (is.null(cache)) stop("Primary M1 fit unavailable for reuse")
        .summarize_cached_sensitivity(cache, spec)
      } else {
        .fit_one_spline_sensitivity(spec$data, spec,
          diagnostic_path = file.path(sensitivity_dir, paste0("4.5_", spec$id, "_covariance.csv")))
      }
    }, error = function(e) {
      list(curve = NULL, theta = NULL,
           status = data.frame(spec = spec$id, status = "skipped_or_failed",
                               n_obs = NA_integer_, n_id = NA_integer_, n_trajectory_id = NA_integer_,
                               n_matching_id = NA_integer_, fitting_min_time = NA_real_,
                               fitting_max_time = NA_real_, prediction_window = spec$prediction_window,
                               spline_df = spec$spline_df, outcome = spec$outcome,
                               rung = NA_character_, singular = NA, convergence = FALSE,
                               optimizer_code = NA_integer_, scaled_gradient = NA_real_,
                               vcov_type = NA_character_, message = conditionMessage(e)),
           convergence = attr(e, "attempts"), variance = NULL,
           support = NULL, visit_support = NULL, set_support = NULL)
    })
  })
  names(ans) <- vapply(specs, `[[`, character(1), "id")
  write.csv(bind_rows(lapply(ans, `[[`, "curve")),
            file.path(sensitivity_dir, "4.5_sensitivity_difference_curves.csv"), row.names = FALSE)
  write.csv(bind_rows(lapply(ans, `[[`, "theta")),
            file.path(sensitivity_dir, "4.5_sensitivity_theta.csv"), row.names = FALSE)
  write.csv(bind_rows(lapply(ans, `[[`, "status")),
            file.path(sensitivity_dir, "4.5_sensitivity_status.csv"), row.names = FALSE)
  write.csv(bind_rows(lapply(ans, `[[`, "convergence")),
            file.path(sensitivity_dir, "4.5_sensitivity_convergence_attempts.csv"), row.names = FALSE)
  write.csv(bind_rows(lapply(ans, `[[`, "variance")),
            file.path(sensitivity_dir, "4.5_sensitivity_variance_components.csv"), row.names = FALSE)
  write.csv(bind_rows(lapply(ans, `[[`, "support")),
            file.path(sensitivity_dir, "4.5_sensitivity_support_by_time_bin.csv"), row.names = FALSE)
  s12 <- ans[["S12"]]
  if (!is.null(s12$visit_support)) {
    write.csv(s12$visit_support, file.path(sensitivity_dir, "4.5_S12_visit_support.csv"), row.names = FALSE)
    write.csv(s12$set_support, file.path(sensitivity_dir, "4.5_S12_matched_set_support.csv"), row.names = FALSE)
  }
  invisible(ans)
}


run_selected_glme_sensitivities <- function(selected, project_dir = getwd()) {
  if ("S7" %in% selected)
    stop("S7 was retired: exact-questionnaire-cycle eligibility is no longer an active sensitivity.")
  engine <- file.path(project_dir, "Code/2_data_analysis/4.0_GLME_spline_functions.R")
  matched <- file.path(project_dir, "Data/riskset_matched_overall_long.rds")
  results <- file.path(project_dir, "Results/cancer/data")
  if (!length(selected) || any(!selected %in% paste0("S", c(1:6, 8:12)))) stop("Invalid sensitivity selection")
  selected <- paste0("S", sort(unique(as.integer(sub("S", "", selected)))))
  primary <- NULL
  if (any(selected %in% c("S3", "S4", "S5", "S6"))) {
    cache_path <- file.path(results, "4.5_primary_cache.rds")
    if (!file.exists(cache_path)) stop("Run the updated 4.5 formal analysis first to save its M1 cache")
    primary <- readRDS(cache_path)
    if (!identical(primary$engine_md5, unname(tools::md5sum(engine))) ||
        !identical(primary$matching_provenance$input_md5, unname(tools::md5sum(matched))))
      stop("Saved M1 cache is stale: engine or matched-input hash differs; rerun 4.5")
  }
  output <- file.path(results, "4.5_sensitivities", paste(selected, collapse = "_"))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(selected = selected, engine_md5 = unname(tools::md5sum(engine)),
    sensitivity_script_md5 = unname(tools::md5sum(file.path(project_dir,
      "Code/2_data_analysis/4.0.1_GLME_sensitivity_analyses.R"))),
    input = validate_matching_provenance(matched), generated_at = Sys.time()),
    file.path(output, "run_configuration.rds"))
  ans <- run_overall_glme_sensitivities(matched,
    output, primary_result = primary, selected = selected)
  status <- bind_rows(lapply(ans, `[[`, "status"))
  if (any(status$status != "fit" | !status$convergence))
    stop("Selected sensitivities have failed or unavailable results; inspect ", output)
  message("Selected sensitivity results saved under: ", output)
  invisible(ans)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    cat("No analyses run. Use --only=S1,S2 (or another selection); --all selects S1-S6 and S8-S11.\n")
  } else if (identical(args, "--all")) {
    run_selected_glme_sensitivities(paste0("S", c(1:6, 8:11)))
  } else if (length(args) == 1L && grepl("^--only=", args)) {
    run_selected_glme_sensitivities(strsplit(sub("^--only=", "", args), ",", fixed = TRUE)[[1]])
  } else stop("Use --only=S1,S2 or --all; no analyses run")
}
