# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.0_GLME_spline_functions.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
#
# Purpose: Shared natural-spline Gaussian LME engine for wrappers 4.1--4.5.
# Fits M0--M3 on all eligible matched rows; support bins govern predictions only.
# Event-study fitting was removed. Original FI construction and matching remain
# upstream. The primary estimand is the post-minus-pre average differential slope.
# Numerical scaling, exact contrasts, covariance diagnostics and failure logs are
# shared with the separately invoked 4.0.1 sensitivity module. Fixed splines use relative time;
# participant slopes use (attained age - 60)/4. Numerical failure alone never
# permits structural simplification; valid CR2 survives small-sample-test failure.
# No analysis runs when sourced.
# Dataset/plot-ready outputs: Results/cancer/data. Render reports with knitr under
# Results/cancer/visuals. Non-results diagnostics belong under Codex.
# The shared fitter defaults to REML; runner 4.0.3 opts into ML for paired AIC
# and its separate, support-gated M1 trajectory display.
# See Documents/Methods/GLME_Natural_Spline_Trajectory_Analysis.md.
# =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(splines)

# A GLME run is valid only when its matched input has passed the matching
# engine's Gate G4 and the provenance record still hashes the exact RDS being
# fitted.  This check is shared by the subgroup wrappers (4.1-4.4) as well as
# the overall wrapper; it prevents silently fitting stale pre-repair cohorts.
source("/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS/Code/2_data_analysis/2.0_matching_provenance.R")

# Prediction-support bin breaks -------------------------------------------
.support_inner_breaks <- seq(-20, 20, by = 2)
.support_rel_breaks   <- c(-Inf, .support_inner_breaks, Inf)
.support_rel_labels   <- c(
  "<= -20 years",
  paste0(
    .support_inner_breaks[-length(.support_inner_breaks)],
    " to ",
    if_else(.support_inner_breaks[-1] > 0,
            paste0("+", .support_inner_breaks[-1]),
            as.character(.support_inner_breaks[-1]))
  ),
  "> +20"
)
.support_ref_label <- "-2 to 0"

add_relative_time_bin <- function(data, rel_time_col = "Age_Centered") {
  data %>%
    mutate(
      rel_time_bin = cut(
        .data[[rel_time_col]],
        breaks = .support_rel_breaks,
        labels = .support_rel_labels,
        right = TRUE,
        include.lowest = TRUE
      ),
      rel_time_bin = factor(rel_time_bin, levels = .support_rel_labels)
    )
}

continuous_support_window <- function(support, window_yrs) {
  if (anyDuplicated(as.character(support$rel_time_bin))) stop("Duplicate support bins")
  support <- data.frame(rel_time_bin = .support_rel_labels) %>%
    left_join(mutate(support, rel_time_bin = as.character(rel_time_bin)), by = "rel_time_bin") %>%
    mutate(support_ok = !is.na(support_ok) & support_ok)
  ref_i <- match(.support_ref_label, as.character(support$rel_time_bin))
  if (is.na(ref_i) || !support$support_ok[[ref_i]]) return(c(NA_real_, NA_real_))
  left_i <- ref_i
  right_i <- ref_i
  while (left_i > 1 && support$support_ok[[left_i - 1]]) left_i <- left_i - 1
  while (right_i < nrow(support) && support$support_ok[[right_i + 1]]) right_i <- right_i + 1
  left_break_i <- match(as.character(support$rel_time_bin[[left_i]]), .support_rel_labels)
  right_break_i <- match(as.character(support$rel_time_bin[[right_i]]), .support_rel_labels)
  lo <- max(-window_yrs, .support_rel_breaks[[left_break_i]])
  hi <- min(window_yrs, .support_rel_breaks[[right_break_i + 1]])
  c(lo, hi)
}

# Small pure helpers keep prediction and validation contracts independently testable.
.bind_model_rows <- function(x) {
  out <- bind_rows(x)
  if (!"Cohort" %in% names(out)) out$Cohort <- character(nrow(out))
  if (!"model_id" %in% names(out)) out$model_id <- character(nrow(out))
  out
}

.exact_slope_contrasts <- function(difference_design, window) {
  stopifnot(length(window) == 1L, is.finite(window), window > 0)
  x <- difference_design(c(-window, 0, window))
  pre <- (x[2, ] - x[1, ]) / window
  post <- (x[3, ] - x[2, ]) / window
  list(pre = pre, post = post, theta = post - pre)
}

.spline_grid_factory <- function(reference, spec, bases, scaling, beta_names) {
  rhs <- .fixed_rhs(.model_time_terms(spec, bases), spec$covars)
  grid <- function(time, group) {
    g <- .make_ref_grid(reference, covars = spec$covars,
                        Age_Centered = time, Group = group) %>%
      mutate(Group = factor(Group, levels = c("Control", "Cancer Case")))
    .apply_model_scaling(.add_model_time_terms(g, spec, bases), scaling)
  }
  design <- function(time, group) model.matrix(rhs, grid(time, group))[, beta_names, drop = FALSE]
  difference <- function(time) design(time, "Cancer Case") - design(time, "Control")
  list(grid = grid, design = design, difference = difference)
}

.inference_constraints <- function(beta_names, time_terms, difference, window = 8,
                                   theta_supported = TRUE) {
  terms <- .spline_interaction_terms(beta_names, time_terms)
  out <- list(omnibus = diag(length(beta_names))[match(terms, beta_names), , drop = FALSE])
  if (theta_supported) {
    exact <- .exact_slope_contrasts(difference, window)
    out$theta <- matrix(exact$theta, nrow = 1)
    out$pre_slope <- matrix(exact$pre, nrow = 1)
  }
  out[vapply(out, nrow, integer(1)) > 0L]
}

.validate_covariance <- function(V, beta_names) {
  if (!identical(dim(V), rep(length(beta_names), 2L)) || !all(is.finite(V)))
    stop("Covariance dimensions or finite values are invalid")
  if (!identical(rownames(V), beta_names) || !identical(colnames(V), beta_names))
    stop("Covariance coefficient names/order differ from fixef")
  if (!isTRUE(all.equal(V, t(V), tolerance = 1e-8))) stop("Covariance is not symmetric")
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8 * max(abs(ev), .Machine$double.eps)) stop("Covariance is not positive semidefinite")
  if (any(diag(V) <= 0)) stop("Nonpositive coefficient variance")
  invisible(TRUE)
}

.validate_wald <- function(x) {
  if (!all(c("Fstat", "df_num", "df_denom", "p_val") %in% names(x)) ||
      !nrow(x) || any(!is.finite(x$Fstat)) || any(x$Fstat < 0) ||
      any(!is.finite(x$df_num) | x$df_num <= 0) ||
      any(!is.finite(x$df_denom) | x$df_denom <= 0) ||
      any(!is.finite(x$p_val) | x$p_val < 0 | x$p_val > 1)) stop("Invalid Wald inference values")
  invisible(TRUE)
}

.nested_covariance_engine <- function(fit, cluster) {
  stopifnot(as.character(packageVersion('clubSandwich')) == '0.6.2',
            is.null(fit@call$weights))
  cluster <- droplevels(factor(cluster))
  groups <- getME(fit, 'flist')
  stopifnot(all(vapply(groups, function(g) identical(as.character(g), as.character(cluster)), logical(1))))
  rows <- split(seq_along(cluster), cluster)
  frame <- model.frame(fit)
  vc <- VarCorr(fit, sigma = 1)
  designs <- lapply(vc, function(G) {
    terms <- colnames(G)
    Z <- matrix(1, nrow(frame), length(terms), dimnames = list(NULL, terms))
    for (v in setdiff(terms, '(Intercept)')) Z[, v] <- frame[[v]]
    Z
  })
  target <- lapply(rows, function(i) {
    T <- diag(length(i))
    for (j in seq_along(vc)) {
      Z <- designs[[j]][i, , drop = FALSE]
      T <- T + Z %*% as.matrix(vc[[j]]) %*% t(Z)
    }
    T
  })
  sparse <- Matrix::tcrossprod(getME(fit, 'Z'), getME(fit, 'Lambdat'))
  probe <- unique(c(1L, length(rows), which.max(lengths(rows)),
                    round(seq(1, length(rows), length.out = 10))))
  for (j in probe) {
    Z <- sparse[rows[[j]], , drop = FALSE]
    original <- as.matrix(Matrix::tcrossprod(Z)) + diag(length(rows[[j]]))
    stopifnot(isTRUE(all.equal(unname(target[[j]]), unname(original), tolerance = 1e-10)))
  }
  rm(sparse)
  weights <- lapply(target, function(T) chol2inv(chol(T)))
  core <- getFromNamespace('vcov_CR', 'clubSandwich')
  scope <- new.env(parent = environment(core))
  scope$weightMatrix <- function(obj, cluster) weights
  scope$targetVariance <- function(obj, cluster) target
  # Reuse row indices instead of scanning all observations once per cluster.
  scope$matrix_list <- function(x, fac, dim) {
    stopifnot(identical(fac, cluster), dim == 'row', is.matrix(x))
    lapply(rows, function(i) x[i, , drop = FALSE])
  }
  environment(core) <- scope
  function(type, form = 'sandwich') core(fit, cluster, type = type,
                                         inverse_var = TRUE, form = form)
}

# Use the validated optimization only for unweighted participant-nested fits.
# Crossed/weighted models keep the package's explicit support checks.
.participant_vcov <- function(model, cluster, type, form = "sandwich") {
  nested <- inherits(model, "lmerMod") && is.null(model@call$weights) &&
    as.character(utils::packageVersion("clubSandwich")) == "0.6.2" &&
    all(vapply(lme4::getME(model, "flist"), function(g)
      identical(as.character(g), as.character(cluster)), logical(1)))
  if (nested) return(.nested_covariance_engine(model, cluster)(type, form))
  clubSandwich::vcovCR(model, cluster = cluster, type = type, form = form)
}

# Covariance selection is independent of small-sample inference feasibility.
# The bound deliberately guards the J-by-J, coefficient-sized intermediate
# arrays in small-sample tests; it is conservative, not a memory guarantee.
.small_sample_feasible <- function(model, cluster) {
  bytes <- 8 * length(unique(cluster))^2 * length(lme4::fixef(model))^2
  is.finite(bytes) && bytes <= getOption("hpfs.small_sample_max_bytes", 256 * 1024^2)
}

get_primary_vcov <- function(model, cluster, cluster_name = "id", constraints = list(),
                             diagnostic_path = NULL) {
  fe <- names(lme4::fixef(model))
  records <- list()
  add <- function(type, stage, status, message = NA_character_) {
    records[[length(records) + 1L]] <<- data.frame(
      type = type, stage = stage, status = status, message = message,
      n_rows = nobs(model), n_clusters = length(unique(cluster)),
      clubSandwich_version = if (requireNamespace("clubSandwich", quietly = TRUE))
        as.character(utils::packageVersion("clubSandwich")) else NA_character_,
      generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    if (!is.null(diagnostic_path)) {
      dir.create(dirname(diagnostic_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(bind_rows(records), diagnostic_path, row.names = FALSE)
    }
  }
  attempt <- function(type, stage, expr) {
    tryCatch(withCallingHandlers({
      value <- force(expr); add(type, stage, "success"); list(ok = TRUE, value = value)
    }, warning = function(w) { add(type, stage, "warning", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { add(type, stage, "error", conditionMessage(e)); list(ok = FALSE) })
  }
  frame <- model.frame(model)
  if (inherits(model, "lmerMod")) {
    if (!is.null(model@call$weights)) add(NA_character_, "model_support", "weighted_requires_review",
      "Weighted lmer covariance is a separate limitation (S11); inspect package errors")
    if (length(lme4::getME(model, "flist")) > 1L) add(NA_character_, "model_support", "crossed_requires_review",
      "Multiple grouping factors: participant-clustered covariance may be unsupported (M3)")
  }
  aligned <- length(cluster) == nrow(frame) && !anyNA(cluster) &&
    cluster_name %in% names(frame) &&
    identical(as.character(cluster), as.character(frame[[cluster_name]]))
  add(NA_character_, "cluster_alignment", if (aligned) "success" else "error",
      if (aligned) NA_character_ else "Cluster vector does not match fitted model rows")
  if (!aligned) stop("Cluster vector does not match fitted model rows; see covariance diagnostics")
  if (aligned && requireNamespace("clubSandwich", quietly = TRUE)) {
    for (type in c("CR2", "CR0")) {
      cv <- attempt(type, "vcovCR", .participant_vcov(model, cluster, type))
      if (!cv$ok) next
      V <- as.matrix(cv$value)
      valid <- attempt(type, "covariance_validation", .validate_covariance(V, fe))
      if (!valid$ok) next
      for (name in names(constraints)) attempt(type, paste0("constraint_validation_", name), {
        C <- constraints[[name]]
        VV <- C %*% V %*% t(C)
        if (any(!is.finite(VV)) || qr(VV)$rank < nrow(C) || any(diag(VV) <= 0))
          stop("Constraint covariance is singular or invalid; this test requires review")
        TRUE
      })
      small_sample <- .small_sample_feasible(model, cluster)
      ct <- list(ok = FALSE, value = NULL)
      if (small_sample) ct <- attempt(type, "coef_test", {
        tab <- as.data.frame(clubSandwich::coef_test(model, vcov = cv$value, test = "Satterthwaite"))
        if (nrow(tab) != length(fe) || !all(c("SE", "df_Satt", "p_Satt") %in% names(tab)) ||
            any(!is.finite(tab$SE) | tab$SE <= 0) ||
            any(!is.finite(tab$df_Satt) | tab$df_Satt <= 0) ||
            any(!is.finite(tab$p_Satt) | tab$p_Satt < 0 | tab$p_Satt > 1)) stop("Invalid coefficient inference values")
        tab
      })
      if (!small_sample) add(type, "small_sample", "skipped_memory_guard",
        "Conservative allocation bound exceeds hpfs.small_sample_max_bytes; retaining robust covariance")
      if (!ct$ok) add(type, "inference", "asymptotic",
        "Using normal/chi-square inference with the selected robust covariance")
      add(type, "selection", "selected")
      return(list(V = V, V_cs = cv$value, robust = TRUE, coefficient_test = ct$value,
                  small_sample = small_sample && ct$ok,
                  type = paste0(type, " (cluster-robust on ", cluster_name, ")"),
                  inference_status = if (ct$ok) "robust_small_sample" else "robust_asymptotic", diagnostics = bind_rows(records), record = add))
    }
  } else if (!requireNamespace("clubSandwich", quietly = TRUE)) add(NA_character_, "package", "error", "clubSandwich not installed")
  V <- as.matrix(vcov(model))
  valid <- attempt("model-based", "covariance_validation", .validate_covariance(V, fe))
  if (!valid$ok) stop("Model-based covariance also invalid; see covariance diagnostics")
  add("model-based", "selection", "fallback_requires_review")
  warning("Robust inference unavailable; selected model-based covariance. Inspect covariance diagnostics.", call. = FALSE)
  list(V = V, V_cs = NULL, robust = FALSE, type = "model-based",
       inference_status = "model_based_fallback_requires_review", diagnostics = bind_rows(records), record = add)
}

coef_table_with_vcov <- function(model, Vobj) {
  if (!is.null(Vobj$coefficient_test)) {
    ct <- Vobj$coefficient_test
  } else {
    beta <- lme4::fixef(model)
    se <- sqrt(diag(Vobj$V))
    ct <- data.frame(Coef. = names(beta), Estimate = unname(beta), SE = se,
                     `d.f. (Satt)` = Inf, `p-val (Satt)` = 2 * pnorm(-abs(beta / se)),
                     check.names = FALSE)
  }
  names(ct)[names(ct) %in% c("Coef.", "Coef")] <- "Term"
  names(ct)[names(ct) %in% c("Estimate", "beta")] <- "Estimate"
  names(ct)[names(ct) %in% c("d.f. (Satt)", "df_Satt")] <- "df"
  names(ct)[names(ct) %in% c("p-val (Satt)", "p_Satt")] <- "p_value"
  ct %>% transmute(
    Term, Estimate, SE, df = as.numeric(df), p_value = as.numeric(p_value),
    CI_low = Estimate - qt(0.975, df = df) * SE,
    CI_high = Estimate + qt(0.975, df = df) * SE,
    infer_method = if (!is.null(Vobj$coefficient_test)) paste(Vobj$type, "Satterthwaite") else paste(Vobj$type, "normal")
  )
}

wald_with_vcov <- function(model, Vobj, constraints, label) {
  if (!is.null(Vobj$V_cs) && isTRUE(Vobj$small_sample)) {
    ans <- tryCatch({
      x <- as.data.frame(clubSandwich::Wald_test(model, constraints = constraints,
                                               vcov = Vobj$V_cs, test = "HTZ"))
      .validate_wald(x); x
    }, error = function(e) e)
    if (inherits(ans, "error") && is.function(Vobj$record))
      Vobj$record(Vobj$type, paste0("downstream_Wald_", label), "error", conditionMessage(ans))
    if (!inherits(ans, "error")) return(ans %>% transmute(
      Test = label, Fstat, df_num, df_denom, p_value = p_val,
      infer_method = paste(Vobj$type, "HTZ")
    ))
    if (is.function(Vobj$record)) Vobj$record(Vobj$type,
      paste0("downstream_Wald_", label), "asymptotic", "Retained covariance after optional HTZ failure")
    # A failed optional test never changes the selected covariance.
  }
  beta <- lme4::fixef(model)
  est <- as.vector(constraints %*% beta)
  VV <- constraints %*% Vobj$V %*% t(constraints)
  q <- qr(VV)$rank
  if (q < nrow(constraints) || any(diag(VV) <= 0)) return(data.frame(
    Test = label, Fstat = NA_real_, df_num = q, df_denom = NA_real_, p_value = NA_real_,
    infer_method = paste(Vobj$type, "constraint covariance invalid")))
  chisq <- as.numeric(t(est) %*% MASS::ginv(VV) %*% est)
  data.frame(Test = label, Fstat = chisq / q, df_num = q,
             df_denom = Inf, p_value = pchisq(chisq, q, lower.tail = FALSE),
             infer_method = paste(Vobj$type, "asymptotic Wald chi-square"))
}

# Joint Wald test that a named subset of fixed-effect coefficients are all 0. ---
joint_wald <- function(model, Vobj, terms, label) {
  beta <- lme4::fixef(model)
  terms <- intersect(terms, names(beta))
  if (length(terms) == 0) {
    return(data.frame(Test = label, Fstat = NA_real_, df_num = 0,
                      df_denom = NA_real_, p_value = NA_real_,
                      infer_method = Vobj$type))
  }
  constraints <- diag(length(beta))[match(terms, names(beta)), , drop = FALSE]
  wald_with_vcov(model, Vobj, constraints, label)
}

# Build a per-cohort prediction grid with covariates held at reference levels. --
.primary_covars <- c("index_age_z", "base_race", "base_marital", "base_living")
.full_extra_covars <- c("base_pckgr", "base_calor", "base_sat",
                        "base_diet_chol", "base_alco")
.full_covars <- c(.primary_covars, .full_extra_covars)

.trajectory_model_specs <- list(
  M0_raw = list(
    label = "M0 raw",
    role = "Minimally adjusted matched contrast",
    covars = "index_age_z",
    time_structure = "natural spline",
    spline_df = 3L,
    basis_key = "m0_spline",
    random = "(1 + .random_time | id)",
    matching_set_random = FALSE
  ),
  M1_primary_spline = list(
    label = "M1 primary spline (df = 3)",
    role = "Primary-adjusted natural-spline model",
    covars = .primary_covars,
    time_structure = "natural spline",
    spline_df = 3L,
    basis_key = "adjusted_spline",
    random = "(1 + .random_time | id)",
    matching_set_random = FALSE
  ),
  M2_full_spline = list(
    label = "M2 full spline (df = 3)",
    role = "Expanded baseline-covariate natural-spline sensitivity",
    covars = .full_covars,
    time_structure = "natural spline",
    spline_df = 3L,
    basis_key = "adjusted_spline",
    random = "(1 + .random_time | id)",
    matching_set_random = FALSE
  ),
  M3_primary_matching_spline = list(
    label = "M3 primary spline + matching_id (df = 3)",
    role = "Primary-adjusted natural spline plus matched-set random intercept",
    covars = .primary_covars,
    time_structure = "natural spline",
    spline_df = 3L,
    basis_key = "adjusted_spline",
    random = "(1 + .random_time | id) + (1 | matching_id)",
    matching_set_random = TRUE
  )
)

.most_common_participant_level <- function(data, column, id_col = "id") {
  factor_levels <- levels(data[[column]])
  if (is.null(factor_levels) || !length(factor_levels)) {
    stop("Factor prediction variable has no defined levels: ", column,
         call. = FALSE)
  }

  values <- data.frame(
    id = if (id_col %in% names(data)) as.character(data[[id_col]]) else seq_len(nrow(data)),
    value = as.character(data[[column]]),
    stringsAsFactors = FALSE
  )
  values <- unique(values[!is.na(values$value), , drop = FALSE])
  counts <- tabulate(match(values$value, factor_levels), nbins = length(factor_levels))
  if (!any(counts > 0L)) return(factor_levels[[1]])

  # which.max uses the factor-level order to break ties deterministically.
  factor_levels[[which.max(counts)]]
}

.make_ref_grid <- function(d, covars = .primary_covars, ...) {
  base_cols <- list()
  for (cv in covars) {
    if (cv == "index_age_z") {
      base_cols[[cv]] <- 0
    } else if (is.factor(d[[cv]])) {
      modal_level <- .most_common_participant_level(d, cv)
      base_cols[[cv]] <- factor(modal_level, levels = levels(d[[cv]]))
    } else if (is.numeric(d[[cv]])) {
      base_cols[[cv]] <- mean(d[[cv]], na.rm = TRUE)
    } else {
      base_cols[[cv]] <- d[[cv]][which(!is.na(d[[cv]]))[1]]
    }
  }
  expand.grid(..., stringsAsFactors = FALSE) %>%
    mutate(!!!base_cols)
}

.fixed_rhs <- function(spline_terms, covars = .primary_covars) {
  if (!length(spline_terms)) stop("Spline terms are required")
  group_time <- paste0("Group * (", paste(spline_terms, collapse = " + "), ")")
  as.formula(paste("~", paste(c(group_time, covars), collapse = " + ")))
}

# Center and scale numeric model columns without changing the fitted model
# space.  Scaling is estimated on the exact model-specific complete-case
# fitting data and must subsequently be applied to its prediction grids.
.scale_model_numeric_columns <- function(data, columns) {
  columns <- intersect(columns, names(data))
  if (!length(columns)) {
    return(list(data = data,
                parameters = data.frame(column = character(), center = numeric(),
                                        scale = numeric(), stringsAsFactors = FALSE)))
  }
  if (!all(vapply(data[columns], is.numeric, logical(1)))) {
    stop("Only numeric columns may be standardized for model fitting.", call. = FALSE)
  }
  centers <- vapply(data[columns], mean, numeric(1), na.rm = TRUE)
  scales <- vapply(data[columns], stats::sd, numeric(1), na.rm = TRUE)
  if (any(!is.finite(centers)) || any(!is.finite(scales)) || any(scales <= 0)) {
    bad <- columns[!is.finite(centers) | !is.finite(scales) | scales <= 0]
    stop("Cannot standardize non-varying or non-finite model columns: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  data[columns] <- Map(function(x, center, scale) (x - center) / scale,
                       data[columns], centers, scales)
  list(data = data,
       parameters = data.frame(column = columns, center = unname(centers),
                               scale = unname(scales), stringsAsFactors = FALSE))
}

.apply_model_scaling <- function(data, parameters) {
  if (!nrow(parameters)) return(data)
  missing_cols <- setdiff(parameters$column, names(data))
  if (length(missing_cols)) {
    stop("Prediction grid is missing scaled model columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  for (i in seq_len(nrow(parameters))) {
    column <- parameters$column[[i]]
    data[[column]] <- (data[[column]] - parameters$center[[i]]) /
      parameters$scale[[i]]
  }
  data
}

.make_spline_basis <- function(time, spline_df, prefix) {
  boundary_knots <- range(time, na.rm = TRUE)
  interior_knots <- as.numeric(stats::quantile(
    time,
    probs = seq(0, 1, length.out = spline_df + 1)[-c(1, spline_df + 1)],
    na.rm = TRUE
  ))
  basis <- splines::ns(time, knots = interior_knots,
                       Boundary.knots = boundary_knots)
  terms <- paste0(prefix, seq_len(ncol(basis)))
  colnames(basis) <- terms
  list(
    basis = basis, terms = terms, spline_df = spline_df,
    knots = attr(basis, "knots"),
    boundary_knots = attr(basis, "Boundary.knots")
  )
}

.model_time_terms <- function(spec, spline_bases) {
  spline_bases[[spec$basis_key]]$terms
}

.add_model_time_terms <- function(data, spec, spline_bases) {
  basis_spec <- spline_bases[[spec$basis_key]]
  values <- predict(basis_spec$basis, data$Age_Centered)
  bind_cols(data, setNames(as.data.frame(values), basis_spec$terms))
}

.assess_lmer_fit <- function(fit, allow_singular = FALSE, gradient_tol = 0.002) {
  opt_code <- fit@optinfo$conv$opt
  opt_code <- if (length(opt_code) == 0 || is.null(opt_code)) 0L else as.integer(opt_code[[1]])
  conv_messages <- fit@optinfo$conv$lme4$messages
  substantive_messages <- conv_messages[!grepl("very different scales|consider rescaling", conv_messages, ignore.case = TRUE)]
  if (allow_singular && length(substantive_messages)) {
    substantive_messages <- substantive_messages[
      !grepl("^boundary \\(singular\\) fit", substantive_messages)
    ]
  }
  conv_messages <- if (length(conv_messages)) paste(conv_messages, collapse = " | ") else NA_character_
  derivs <- fit@optinfo$derivs
  scaled_gradient <- NA_real_
  if (!is.null(derivs$gradient) && !is.null(derivs$Hessian)) {
    scaled_gradient <- tryCatch({
      chol_h <- chol(derivs$Hessian)
      max(abs(backsolve(chol_h, derivs$gradient)))
    }, error = function(e) Inf)
  }
  singular <- lme4::isSingular(fit)
  finite <- all(is.finite(lme4::fixef(fit))) &&
    all(is.finite(unlist(lme4::VarCorr(fit)))) &&
    is.finite(stats::sigma(fit))
  converged <- identical(opt_code, 0L) && length(substantive_messages) == 0 &&
    is.finite(scaled_gradient) && scaled_gradient <= gradient_tol && finite &&
    (allow_singular || !singular)
  list(converged = converged, optimizer_code = opt_code,
       convergence_messages = conv_messages, scaled_gradient = scaled_gradient,
       singular = singular, finite = finite)
}

.participant_boundary <- function(fit, tolerance = 1e-4) {
  blocks <- lme4::VarCorr(fit, sigma = 1)
  blocks <- blocks[grepl("^id($|\\.)", names(blocks))]
  slope <- unlist(lapply(blocks, function(G) diag(G)[names(diag(G)) == ".random_time"]))
  singular <- any(vapply(blocks, function(G)
    min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) <= tolerance^2, logical(1)))
  list(slope_boundary = length(slope) > 0 && all(slope <= tolerance^2),
       correlation_boundary = singular && any(lengths(blocks) > 1L),
       participant_supported = length(blocks) > 0 && !singular)
}

.matching_zero_variance <- function(fit, tolerance = 1e-7) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  matching_var <- vc$vcov[vc$grp == "matching_id" &
                            vc$var1 == "(Intercept)" & is.na(vc$var2)]
  id_var <- vc$vcov[vc$grp == "id" &
                      vc$var1 == "(Intercept)" & is.na(vc$var2)]
  list(
    permitted = length(matching_var) == 1L && length(id_var) == 1L &&
      is.finite(matching_var) && is.finite(id_var) &&
      matching_var <= tolerance && id_var > tolerance &&
      .participant_boundary(fit)$participant_supported,
    matching_id_variance = if (length(matching_var) == 1L) matching_var else NA_real_
  )
}

.lmer_attempts <- function(formula, numerical_retries = FALSE, simplify = TRUE) {
  formula_text <- paste(deparse(formula), collapse = "")
  slope_uncorrelated <- as.formula(gsub("\\(1 \\+ .random_time \\| id\\)",
                                        "(1 + .random_time || id)", formula_text))
  intercept_only <- as.formula(gsub("\\(1 \\+ .random_time \\|\\| id\\)|\\(1 \\+ .random_time \\| id\\)",
                                    "(1 | id)", formula_text))
  attempts <- list(
    list(
      rung = "correlated random slope",
      formula = formula,
      control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE,
                            optCtrl = list(maxfun = 2e5))
    ),
    list(
      rung = "uncorrelated random slope",
      formula = slope_uncorrelated,
      control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE,
                            optCtrl = list(maxfun = 2e5))
    ),
    list(
      rung = "random intercept only",
      formula = intercept_only,
      control = lmerControl(optimizer = "bobyqa", calc.derivs = TRUE,
                            optCtrl = list(maxfun = 2e5))
    )
  )

  if (!simplify) attempts <- attempts[1]
  formula_keys <- vapply(attempts, function(x) paste(deparse(x$formula), collapse = ""), character(1))
  attempts <- attempts[!duplicated(formula_keys)]
  for (i in seq_along(attempts)) {
    txt <- paste(deparse(attempts[[i]]$formula), collapse = "")
    attempts[[i]]$rung <- if (!grepl(".random_time", txt, fixed = TRUE)) "random intercept only" else
      if (grepl("||", txt, fixed = TRUE)) "uncorrelated random slope" else "correlated random slope"
    attempts[[i]]$optimizer <- "bobyqa"
  }
  if (numerical_retries) attempts <- unlist(lapply(attempts, function(at) {
    retry <- at
    retry$optimizer <- "nloptwrap"
    retry$control <- lmerControl(optimizer = "nloptwrap", calc.derivs = TRUE,
      optCtrl = list(maxeval = 2e5, ftol_abs = 1e-8, xtol_abs = 1e-8))
    list(at, retry)
  }), recursive = FALSE)
  attempts
}

.set_random_clock <- function(data, random_clock = "attained_age", random_time_scale = 4) {
  if (random_clock == "attained_age") {
    if (!"age_at_cycle" %in% names(data) || any(!is.finite(data$age_at_cycle)))
      stop("Finite age_at_cycle is required for attained-age random slopes")
    if (all(c("id", "cycle") %in% names(data))) {
      copies <- split(data$age_at_cycle, interaction(data$id, data$cycle, drop = TRUE))
      if (any(vapply(copies, function(x) diff(range(x)) > 1e-8, logical(1))))
        stop("Attained age differs across copies of the same participant-cycle")
    }
    data$.random_time <- (data$age_at_cycle - 60) / random_time_scale
  } else if (random_clock == "relative_time") {
    data$.random_time <- data$Age_Centered / random_time_scale
  } else stop("Unknown random clock")
  data
}

# Isolated call boundary permits orchestration tests without running a model.
.fit_lmer_candidate <- function(formula, data, control, weights = NULL,
                                reml = TRUE) {
  if (is.null(weights)) lme4::lmer(formula, data = data, REML = reml,
                                  control = control, na.action = na.fail)
  else lme4::lmer(formula, data = data, REML = reml, control = control,
                  weights = weights, na.action = na.fail)
}

.fit_lmer_with_ladder <- function(formula, data, model_label, weights = NULL,
                                  allow_matching_boundary = FALSE,
                                  numerical_retries = TRUE, simplify = TRUE,
                                  random_clock = "attained_age", random_time_scale = 4,
                                  diagnostic_path = NULL, reml = TRUE) {
  stopifnot(is.finite(random_time_scale), random_time_scale > 0)
  data <- .set_random_clock(data, random_clock, random_time_scale)
  formula_text <- paste(deparse(formula), collapse = "")
  formula_text <- gsub("Age_Centered |", ".random_time |", formula_text, fixed = TRUE)
  formula <- as.formula(formula_text, env = environment(formula))
  attempts <- .lmer_attempts(formula, numerical_retries, simplify)
  last_error <- NULL
  attempt_log <- list()
  save_attempts <- function() {
    if (!is.null(diagnostic_path)) {
      dir.create(dirname(diagnostic_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(bind_rows(attempt_log), diagnostic_path, row.names = FALSE)
    }
  }
  on.exit(save_attempts(), add = TRUE)
  supported_fallbacks <- character()
  initial_rung <- attempts[[1]]$rung
  for (at in attempts) {
    if (at$rung != initial_rung && !at$rung %in% supported_fallbacks) next
    warnings_seen <- character()
    fit_call <- function() .fit_lmer_candidate(at$formula, data, at$control,
                                                weights, reml = reml)
    fit <- tryCatch(withCallingHandlers(
      fit_call(),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(fit)) {
      attempt_log[[length(attempt_log) + 1L]] <- data.frame(
        model_label = model_label, rung = at$rung, optimizer = at$optimizer,
        random_clock = random_clock, random_time_center = if (random_clock == "attained_age") 60 else 0,
        random_time_scale = random_time_scale, accepted = FALSE,
        optimizer_code = NA_integer_, convergence_messages = NA_character_,
        scaled_gradient = NA_real_, singular = NA, finite = NA,
        warnings = paste(warnings_seen, collapse = " | "), error = last_error
      )
      next
    }
    boundary_matching <- if (isTRUE(allow_matching_boundary)) {
      .matching_zero_variance(fit)
    } else {
      list(permitted = FALSE, matching_id_variance = NA_real_)
    }
    check <- .assess_lmer_fit(fit, allow_singular = boundary_matching$permitted)
    # Scaling warnings alone do not justify changing the random structure.
    structural <- .participant_boundary(fit)
    numerical <- .assess_lmer_fit(fit, allow_singular = TRUE)
    if (numerical$converged) {
      if (structural$correlation_boundary || structural$slope_boundary)
        supported_fallbacks <- union(supported_fallbacks, "uncorrelated random slope")
      if (at$rung == "uncorrelated random slope" && structural$slope_boundary)
        supported_fallbacks <- union(supported_fallbacks, "random intercept only")
    }
    substantive_warnings <- warnings_seen[!grepl("very different scales|consider rescaling", warnings_seen, ignore.case = TRUE)]
    if (boundary_matching$permitted && length(substantive_warnings)) {
      substantive_warnings <- substantive_warnings[
        !grepl("^boundary \\(singular\\) fit", substantive_warnings)
      ]
    }
    check$converged <- check$converged && length(substantive_warnings) == 0
    attempt_log[[length(attempt_log) + 1L]] <- data.frame(
      model_label = model_label, rung = at$rung, optimizer = at$optimizer,
        random_clock = random_clock, random_time_center = if (random_clock == "attained_age") 60 else 0,
        random_time_scale = random_time_scale, accepted = check$converged,
      optimizer_code = check$optimizer_code,
      convergence_messages = check$convergence_messages,
      scaled_gradient = check$scaled_gradient, singular = check$singular,
      finite = check$finite, warnings = paste(warnings_seen, collapse = " | "),
      fallback_evidence = paste(supported_fallbacks, collapse = " | "),
      error = NA_character_, logLik = as.numeric(logLik(fit)),
      variance_components = paste(capture.output(print(VarCorr(fit))), collapse = " | "),
      fixed_effects = paste(names(fixef(fit)), signif(fixef(fit), 10), collapse = " | ")
    )
    save_attempts()
    if (check$converged) {
      return(c(list(fit = fit, rung = at$rung, error = NA_character_,
                    attempts = bind_rows(attempt_log),
                    boundary_matching_variance_zero = boundary_matching$permitted,
                    matching_id_variance = boundary_matching$matching_id_variance), check))
    }
    last_error <- paste0(at$rung, " failed strict acceptance checks")
  }

  err <- simpleError(paste0("Could not fit ", model_label, ": ", last_error))
  attr(err, "attempts") <- bind_rows(attempt_log)
  stop(err)
}

.variance_components <- function(model, cohort, model_id) {
  vc <- as.data.frame(lme4::VarCorr(model))
  vc %>% transmute(
    Cohort = cohort, model_id = model_id, grouping_factor = grp,
    term_1 = var1, term_2 = var2, variance = vcov, sd_or_correlation = sdcor
  )
}

.model_diagnostics <- function(model, data, cohort, model_id) {
  rr <- residuals(model)
  ff <- fitted(model)
  residual_frame <- data.frame(Group = as.character(data$Group), residual = rr,
                               fitted = ff, Age_Centered = data$Age_Centered)
  summarize_residuals <- function(x, group_label) {
    qq <- quantile(x$residual, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1),
                   na.rm = TRUE)
    data.frame(
      Cohort = cohort, model_id = model_id, Group = group_label,
      n = nrow(x), residual_mean = mean(x$residual), residual_sd = sd(x$residual),
      residual_fitted_correlation = cor(x$residual, x$fitted),
      residual_time_correlation = cor(x$residual, x$Age_Centered),
      q0 = qq[[1]], q01 = qq[[2]], q05 = qq[[3]], q25 = qq[[4]],
      q50 = qq[[5]], q75 = qq[[6]], q95 = qq[[7]], q99 = qq[[8]], q100 = qq[[9]]
    )
  }
  residual_summary <- bind_rows(
    summarize_residuals(residual_frame, "All"),
    lapply(split(residual_frame, residual_frame$Group), function(x) {
      summarize_residuals(x, unique(x$Group))
    })
  )
  random_effects <- ranef(model, condVar = FALSE)
  re_summary <- bind_rows(lapply(names(random_effects), function(grp) {
    re <- random_effects[[grp]]
    bind_rows(lapply(names(re), function(term) {
      z <- re[[term]]
      data.frame(Cohort = cohort, model_id = model_id,
                 grouping_factor = grp, term = term, n = length(z),
                 mean = mean(z), sd = sd(z), q01 = quantile(z, .01),
                 q50 = median(z), q99 = quantile(z, .99))
    }))
  }))
  list(residual = residual_summary, random_effect = re_summary)
}

.spline_interaction_terms <- function(beta_names, sterms) {
  beta_names[grepl(":", beta_names) & grepl("GroupCancer Case", beta_names) &
               Reduce(`|`, lapply(sterms, grepl, x = beta_names))]
}

run_spline_analysis <- function(matched_path,
                                           results_dir,
                                           visuals_dir,
                                           out_prefix,
                                           analysis_title,
                                           builder_script,
                                           window_yrs   = 20,
                                           spline_df    = 3,
                                           min_case_bin = 50,
                                           min_ctrl_bin = 250,
                                           run_influence = FALSE) {

  needed_cols <- c("Cohort", "Group", "Age_Centered", "index_age_z",
                   "base_race", "base_marital", "base_living",
                   "id", "match_set", "role", "cycle", "fi_score_nocancer")

  if (!file.exists(matched_path)) {
    stop("Matched dataset not found at ", matched_path, ". Run ", builder_script, " first.")
  }

  matching_provenance <- validate_matching_provenance(matched_path)
  matched_long <- readRDS(matched_path)
  assignment_ledger <- readRDS(matching_provenance$assignment_path)
  missing_cols <- setdiff(needed_cols, names(matched_long))
  if (length(missing_cols) > 0) {
    stop("Matched dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!("matching_id" %in% names(matched_long))) {
    matched_long <- matched_long %>%
      mutate(matching_id = paste(Cohort, match_set, sep = "__"))
  }
  if (!("trajectory_id" %in% names(matched_long))) {
    matched_long <- matched_long %>%
      mutate(trajectory_id = paste(Cohort, match_set, id, role, sep = "__"))
  }
  if (!("post_own_cancer" %in% names(matched_long))) {
    stop("Matched dataset is missing post_own_cancer. Rebuild it with 2.0_riskset_matching_functions.R before fitting GLME models.")
  }

  if ("base_pckgr" %in% names(matched_long)) matched_long$base_pckgr <- factor(matched_long$base_pckgr)
  if (anyNA(matched_long[c("id", "Group", "matching_id", "trajectory_id", "cycle", "post_own_cancer")])) {
    stop("Missing identifiers, group, cycle or own-cancer censoring flag")
  }
  if (!all(matched_long$Group %in% c("Control", "Cancer Case"))) stop("Unexpected Group values")
  if (any(!is.na(matched_long$fi_score_nocancer) &
          (!is.finite(matched_long$fi_score_nocancer) | matched_long$fi_score_nocancer < 0 |
             matched_long$fi_score_nocancer > 1)) ||
      any(!is.na(matched_long$Age_Centered) & !is.finite(matched_long$Age_Centered)))
    stop("FI must be finite in [0,1] and observed relative time must be finite")
  n_post_own_cancer_rows <- sum(matched_long$post_own_cancer %in% TRUE, na.rm = TRUE)
  censoring_log <- matched_long %>%
    group_by(Cohort, Group, role) %>%
    summarize(
      n_rows_before_own_cancer_censoring = n(),
      n_rows_removed_post_own_cancer = sum(post_own_cancer),
      n_rows_after_own_cancer_censoring = sum(!post_own_cancer),
      n_assignments_with_rows_removed = n_distinct(trajectory_id[post_own_cancer]),
      .groups = "drop"
    )

  matched_long <- matched_long %>%
    mutate(
      Cohort    = droplevels(factor(Cohort)),
      Group     = factor(Group, levels = c("Control", "Cancer Case")),
      id        = factor(id),
      match_set = factor(match_set),
      matching_id = factor(matching_id),
      trajectory_id = factor(trajectory_id),
      base_race = factor(base_race),
      base_marital = factor(base_marital),
      base_living = factor(base_living)
    ) %>%
    # Future cases remain eligible controls at an earlier index, but primary
    # follow-up ends at that control's own subsequent cancer diagnosis.
    filter(!post_own_cancer, !is.na(fi_score_nocancer), !is.na(Age_Centered)) %>%
    add_relative_time_bin()

  if (anyDuplicated(paste(matched_long$trajectory_id, matched_long$cycle, sep = "__"))) {
    stop("Duplicated trajectory_id x cycle rows remain after primary filtering.")
  }

  message("Primary own-cancer censoring removed ", n_post_own_cancer_rows,
          " control rows at/after a later own cancer diagnosis.")

  cohorts <- levels(matched_long$Cohort)
  group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")
  model_specification <- bind_rows(lapply(names(.trajectory_model_specs), function(model_id) {
    spec <- .trajectory_model_specs[[model_id]]
    data.frame(
      model_id = model_id, model_label = spec$label, model_role = spec$role,
      time_structure = spec$time_structure, spline_df = spec$spline_df,
      random_clock = "attained_age", random_time_center = 60, random_time_scale = 4,
      covariates = paste(spec$covars, collapse = " + "),
      random_effects = spec$random,
      stringsAsFactors = FALSE
    )
  }))

  support_all   <- list()
  filter_all    <- list()
  metadata_all  <- list()
  sp_pred_all   <- list()
  sp_diff_all   <- list()
  sp_deriv_all  <- list()
  sp_coef_all   <- list()
  sp_theta_all  <- list()
  sp_status_all <- list()
  convergence_all <- list()
  variance_all <- list()
  residual_diag_all <- list()
  random_effect_diag_all <- list()
  bounded_diag_all <- list()
  omnibus_all <- list()
  influence_all <- list()
  scaling_all <- list()
  sp_models     <- list()
  sp_contexts <- list()

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  assignment_flow <- bind_rows(lapply(names(.trajectory_model_specs), function(model_id) {
    spec <- .trajectory_model_specs[[model_id]]
    missing_covars <- setdiff(spec$covars, names(matched_long))
    model_rows <- if (length(missing_covars)) matched_long[FALSE, ] else
      matched_long %>% filter(if_all(all_of(spec$covars), ~ !is.na(.x)))
    assignment_ledger %>% group_by(Cohort, Group) %>% summarize(
      n_matched_assignments = n(),
      n_fi_contributing_assignments = sum(trajectory_id %in% matched_long$trajectory_id),
      n_model_complete_assignments = if (length(missing_covars)) NA_integer_ else
        sum(trajectory_id %in% model_rows$trajectory_id),
      n_without_usable_fi = n_matched_assignments - n_fi_contributing_assignments,
      n_removed_model_completeness = n_fi_contributing_assignments - n_model_complete_assignments,
      model_id = model_id,
      status = if (length(missing_covars)) "missing_covariate_columns" else "assessed",
      .groups = "drop")
  }))
  write.csv(assignment_flow, file.path(results_dir, paste0(out_prefix, "_assignment_flow.csv")), row.names = FALSE)
  checkpoint <- function() {
    write.csv(bind_rows(unlist(sp_status_all, recursive = FALSE)),
              file.path(results_dir, paste0(out_prefix, "_spline_model_status.csv")), row.names = FALSE)
    write.csv(bind_rows(convergence_all),
              file.path(results_dir, paste0(out_prefix, "_convergence_attempts.csv")), row.names = FALSE)
    write.csv(bind_rows(filter_all),
              file.path(results_dir, paste0(out_prefix, "_model_filter_log.csv")), row.names = FALSE)
    write.csv(bind_rows(scaling_all),
              file.path(results_dir, paste0(out_prefix, "_spline_scaling.csv")), row.names = FALSE)
    saveRDS(metadata_all, file.path(results_dir, paste0(out_prefix, "_spline_metadata.rds")))
  }
  on.exit(checkpoint(), add = TRUE)

  for (ch in cohorts) {

    d_base <- matched_long %>%
      filter(Cohort == ch) %>%
      droplevels()

    filter_all[[ch]] <- bind_rows(lapply(names(.trajectory_model_specs), function(model_id) {
      spec <- .trajectory_model_specs[[model_id]]
      if (!all(spec$covars %in% names(d_base))) return(NULL)
      d_model_check <- d_base %>% filter(if_all(all_of(spec$covars), ~ !is.na(.x)))
      before <- d_base %>% distinct(Group, trajectory_id) %>%
        count(Group, name = "n_before_complete_case")
      after <- d_model_check %>% distinct(Group, trajectory_id) %>%
        count(Group, name = "n_after_complete_case")
      before %>%
        left_join(after, by = "Group") %>%
        mutate(
          Cohort = ch, model_id = model_id,
          n_after_complete_case = if_else(is.na(n_after_complete_case), 0L,
                                          n_after_complete_case),
          n_dropped_complete_case = n_before_complete_case - n_after_complete_case,
          pct_dropped_complete_case = 100 * n_dropped_complete_case / n_before_complete_case
        )
    }))
    if (any(filter_all[[ch]]$pct_dropped_complete_case > 5, na.rm = TRUE)) {
      warning("Complete-case filtering removes more than 5% of assignments in at least one arm/model for ",
              ch, ". See the model filter log.", call. = FALSE)
    }

    # ---------------------------- 2) trajectory GLME/LME model set -------------
    # M0-M3 share a df-3 spline basis estimated on the primary complete-case
    # sample.
    d_spline_basis <- d_base %>%
      filter(if_all(all_of(.primary_covars), ~ !is.na(.x)))
    if (!nrow(d_spline_basis)) stop("No primary complete-case rows for cohort: ", ch)
    adjusted_spline_df <- 3L
    spline_bases <- list(
      m0_spline = .make_spline_basis(d_spline_basis$Age_Centered,
                                      spline_df = spline_df, prefix = "M0S"),
      adjusted_spline = .make_spline_basis(d_spline_basis$Age_Centered,
                                            spline_df = adjusted_spline_df, prefix = "S")
    )
    metadata_all[[ch]] <- list(
      schema_version = 3L,
      input_md5 = matching_provenance$input_md5,
      assignment_md5 = matching_provenance$assignment_md5,
      inference_policy = "CR2 then CR0 then labeled model-based fallback; policy review pending",
      event_study = FALSE,
      influence_refits_requested = run_influence,
      random_clock = "attained_age", random_time_center = 60, random_time_scale = 4,
      generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      cohort = ch,
      outcome = "fi_score_nocancer",
      primary_model_id = "M1_primary_spline",
      primary_time_structure = "natural spline",
      spline_reference_model_id = "M1_primary_spline",
      model_specification = lapply(.trajectory_model_specs, function(x) {
        x[c("label", "role", "covars", "time_structure", "spline_df", "basis_key",
            "random", "matching_set_random")]
      }),
      spline_df = adjusted_spline_df,
      knots = spline_bases$adjusted_spline$knots,
      boundary_knots = spline_bases$adjusted_spline$boundary_knots,
      spline_bases = lapply(spline_bases, function(x) x[c("terms", "spline_df", "knots", "boundary_knots")]),
      fitting_time_range = range(d_spline_basis$Age_Centered, na.rm = TRUE),
      n_obs = nrow(d_spline_basis),
      n_id = n_distinct(d_spline_basis$id),
      n_trajectory_id = n_distinct(d_spline_basis$trajectory_id),
      n_matching_id = n_distinct(d_spline_basis$matching_id),
      prediction_window_requested = c(-window_yrs, window_yrs),
      theta_windows = list(pre = c(-8, 0), post = c(0, 8)),
      primary_covariates = .primary_covars,
      own_cancer_censored_rows = n_post_own_cancer_rows
    )
    d_sp_all <- d_base %>%
      .add_model_time_terms(.trajectory_model_specs$M0_raw, spline_bases) %>%
      # M1, M2, and M3 all use the same adjusted_spline basis (S1-S3).
      # Add it once; appending it again for M2 creates duplicate names and
      # causes tibble name repair to hide the columns expected by the models.
      .add_model_time_terms(.trajectory_model_specs$M1_primary_spline, spline_bases)

    sp_models[[ch]] <- list()
    sp_pred_all[[ch]] <- list()
    sp_diff_all[[ch]] <- list()
    sp_deriv_all[[ch]] <- list()
    sp_coef_all[[ch]] <- list()
    sp_theta_all[[ch]] <- list()
    sp_status_all[[ch]] <- list()
    metadata_all[[ch]]$model_scaling <- list()

    for (model_id in names(.trajectory_model_specs)) {
      spec <- .trajectory_model_specs[[model_id]]
      time_terms <- .model_time_terms(spec, spline_bases)
      missing_model_cols <- setdiff(c(spec$covars, time_terms), names(d_sp_all))
      if (length(missing_model_cols) > 0) {
        msg <- paste("Skipped", spec$label, "for", ch,
                     "- missing columns:",
                     paste(missing_model_cols, collapse = ", "),
                     ". Rerun upstream cleaning/matching if baseline nutrition variables are expected.")
        message(msg)
        sp_status_all[[ch]][[model_id]] <- data.frame(
          Cohort = ch, model_id = model_id, model_label = spec$label,
          model_role = spec$role, status = "skipped_missing_columns",
          n_obs = NA_integer_, n_id = NA_integer_, n_matching_id = NA_integer_,
          fitting_min_time = NA_real_, fitting_max_time = NA_real_,
          rung = NA_character_, singular = NA, vcov_type = NA_character_,
          message = msg
        )
        next
      }

      d_model <- d_sp_all %>%
        filter(if_all(all_of(spec$covars), ~ !is.na(.x))) %>%
        droplevels()

      if (length(unique(d_model$Group)) < 2 || nrow(d_model) == 0) {
        msg <- paste("Skipped", spec$label, "for", ch,
                     "- insufficient two-arm support after model-specific complete-case filtering.")
        message(msg)
        sp_status_all[[ch]][[model_id]] <- data.frame(
          Cohort = ch, model_id = model_id, model_label = spec$label,
          model_role = spec$role, status = "skipped_insufficient_support",
          n_obs = nrow(d_model), n_id = n_distinct(d_model$id),
          n_matching_id = n_distinct(d_model$matching_id),
          fitting_min_time = suppressWarnings(min(d_model$Age_Centered, na.rm = TRUE)),
          fitting_max_time = suppressWarnings(max(d_model$Age_Centered, na.rm = TRUE)),
          rung = NA_character_, singular = NA, vcov_type = NA_character_,
          message = msg
        )
        next
      }

      reference_profile <- .make_ref_grid(d_model, covars = spec$covars, .reference = 1)
      # This is a purely numerical reparameterization: centering/scaling the
      # spline basis retains the same fixed-effect column space, while M2's
      # continuous dietary covariates need commensurate scales with the spline
      # and factor columns.  Store every parameter for exact grid reconstruction.
      scaling_columns <- time_terms
      if (identical(model_id, "M2_full_spline")) {
        scaling_columns <- c(scaling_columns,
                             c("base_calor", "base_sat", "base_diet_chol", "base_alco"))
      }
      scaling <- .scale_model_numeric_columns(d_model, scaling_columns)
      d_model <- scaling$data
      metadata_all[[ch]]$model_scaling[[model_id]] <- scaling$parameters
      scaling_all[[paste(ch, model_id, sep = "__")]] <-
        scaling$parameters %>%
        mutate(Cohort = ch, model_id = model_id, model_label = spec$label,
               .before = 1)

      model_support <- d_model %>%
        mutate(rel_time_bin = factor(rel_time_bin, levels = .support_rel_labels)) %>%
        group_by(rel_time_bin, .drop = FALSE) %>%
        summarize(
          n_case = n_distinct(id[Group == "Cancer Case"]),
          n_ctrl = n_distinct(id[Group == "Control"]),
          .groups = "drop"
        ) %>%
        mutate(
          support_ok = n_case >= min_case_bin & n_ctrl >= min_ctrl_bin,
          Cohort = ch, analysis = paste0("spline_", model_id)
        )
      support_all[[paste(ch, model_id, sep = "__")]] <- model_support
      prediction_window <- continuous_support_window(model_support, window_yrs)
      if (anyNA(prediction_window) || prediction_window[[1]] >= prediction_window[[2]]) {
        msg <- paste("Skipped", spec$label, "for", ch,
                     "- no continuous two-arm support meeting the configured thresholds.")
        message(msg)
        sp_status_all[[ch]][[model_id]] <- data.frame(
          Cohort = ch, model_id = model_id, model_label = spec$label,
          model_role = spec$role, status = "skipped_insufficient_supported_window",
          n_obs = nrow(d_model), n_id = n_distinct(d_model$id),
          n_matching_id = n_distinct(d_model$matching_id),
          fitting_min_time = min(d_model$Age_Centered, na.rm = TRUE),
          fitting_max_time = max(d_model$Age_Centered, na.rm = TRUE),
          prediction_min_time = NA_real_, prediction_max_time = NA_real_,
          rung = NA_character_, singular = NA, vcov_type = NA_character_,
          message = msg
        )
        next
      }

      sp_fixed_rhs <- .fixed_rhs(spline_terms = time_terms, covars = spec$covars)
      random_rhs <- spec$random
      model_formula <- as.formula(paste(
        "fi_score_nocancer",
        paste(deparse(sp_fixed_rhs), collapse = ""),
        "+",
        random_rhs
      ))

      fit_info <- tryCatch(
        .fit_lmer_with_ladder(
          model_formula, d_model, spec$label,
          allow_matching_boundary = isTRUE(spec$matching_set_random)
        ),
        error = function(e) e
      )
      if (inherits(fit_info, "error")) {
        convergence_all[[paste(ch, model_id, sep = "__")]] <-
          attr(fit_info, "attempts") %>%
          mutate(Cohort = ch, model_id = model_id, .before = 1)
        msg <- conditionMessage(fit_info)
        sp_status_all[[ch]][[model_id]] <- data.frame(
          Cohort = ch, model_id = model_id, model_label = spec$label,
          model_role = spec$role, status = "failed_strict_convergence",
          n_obs = nrow(d_model), n_id = n_distinct(d_model$id),
          n_matching_id = n_distinct(d_model$matching_id),
          fitting_min_time = min(d_model$Age_Centered, na.rm = TRUE),
          fitting_max_time = max(d_model$Age_Centered, na.rm = TRUE),
          prediction_min_time = prediction_window[[1]],
          prediction_max_time = prediction_window[[2]],
          rung = NA_character_, singular = NA, convergence = FALSE,
          optimizer_code = NA_integer_, scaled_gradient = NA_real_,
          convergence_messages = NA_character_, vcov_type = NA_character_,
          message = msg
        )
        next
      }
      sp_status_all[[ch]][[model_id]] <- data.frame(
        Cohort = ch, model_id = model_id, model_label = spec$label,
        status = "fit_pending_inference", convergence = fit_info$converged,
        rung = fit_info$rung, message = NA_character_)
      m_sp <- fit_info$fit
      convergence_all[[paste(ch, model_id, sep = "__")]] <-
        fit_info$attempts %>% mutate(Cohort = ch, model_id = model_id, .before = 1)
      sp_models[[ch]][[model_id]] <- m_sp

      beta_sp <- fixef(m_sp)
      grids <- .spline_grid_factory(reference_profile, spec, spline_bases,
                                     scaling$parameters, names(beta_sp))
      constraints <- .inference_constraints(names(beta_sp), time_terms, grids$difference,
        theta_supported = prediction_window[[1]] <= -8 && prediction_window[[2]] >= 8)
      Vobj <- get_primary_vcov(m_sp, model.frame(m_sp)$id, constraints = constraints,
        diagnostic_path = file.path(results_dir, paste0(out_prefix, "_covariance_",
          gsub("[^A-Za-z0-9]+", "_", ch), "_", model_id, ".csv")))
      if (identical(model_id, "M1_primary_spline")) sp_contexts[[ch]] <- list(
        fit = m_sp, grids = grids, Vobj = Vobj, window = prediction_window,
        support = model_support, rung = fit_info$rung)
      Vsp <- Vobj$V
      se_sp <- sqrt(diag(Vsp))

      variance_all[[paste(ch, model_id, sep = "__")]] <-
        .variance_components(m_sp, ch, model_id)
      dg <- .model_diagnostics(m_sp, d_model, ch, model_id)
      residual_diag_all[[paste(ch, model_id, sep = "__")]] <- dg$residual
      random_effect_diag_all[[paste(ch, model_id, sep = "__")]] <- dg$random_effect

      sp_coef_all[[ch]][[model_id]] <- coef_table_with_vcov(m_sp, Vobj) %>%
        mutate(Cohort = ch, model_id = model_id, model_label = spec$label,
               model_role = spec$role, vcov_type = Vobj$type, .before = 1)

      make_sp_grid <- grids$grid

      pred_age <- seq(prediction_window[[1]], prediction_window[[2]], by = 0.25)
      grid_sp <- .make_ref_grid(
        reference_profile,
        covars = spec$covars,
        Age_Centered = pred_age,
        Group = c("Control", "Cancer Case")
      ) %>%
        mutate(Group = factor(Group, levels = c("Control", "Cancer Case")))
      grid_sp <- grid_sp %>%
        .add_model_time_terms(spec, spline_bases) %>%
        .apply_model_scaling(scaling$parameters)

      Xsp <- model.matrix(sp_fixed_rhs, data = grid_sp)[, names(beta_sp), drop = FALSE]
      ps  <- as.vector(Xsp %*% beta_sp)
      sds <- sqrt(rowSums((Xsp %*% Vsp) * Xsp))
      grid_sp$pred <- ps
      grid_sp$se <- sds
      grid_sp$lwr <- ps - 1.96 * sds
      grid_sp$upr <- ps + 1.96 * sds
      grid_sp$Cohort <- ch
      grid_sp$model_id <- model_id
      grid_sp$model_label <- spec$label
      grid_sp$model_role <- spec$role
      grid_sp$vcov_type <- Vobj$type
      bounded_diag_all[[paste(ch, model_id, sep = "__")]] <- data.frame(
        Cohort = ch, model_id = model_id, n_predictions = nrow(grid_sp),
        n_point_outside_0_1 = sum(grid_sp$pred < 0 | grid_sp$pred > 1),
        n_ci_outside_0_1 = sum(grid_sp$lwr < 0 | grid_sp$upr > 1),
        min_prediction = min(grid_sp$pred), max_prediction = max(grid_sp$pred),
        min_ci = min(grid_sp$lwr), max_ci = max(grid_sp$upr)
      )
      sp_pred_all[[ch]][[model_id]] <- grid_sp[, c(
        "Cohort", "model_id", "model_label", "model_role", "Group",
        "Age_Centered", "pred", "se", "lwr", "upr", "vcov_type"
      )]

      grid_case <- make_sp_grid(pred_age, "Cancer Case")
      grid_ctrl <- make_sp_grid(pred_age, "Control")
      X_case <- model.matrix(sp_fixed_rhs, data = grid_case)[, names(beta_sp), drop = FALSE]
      X_ctrl <- model.matrix(sp_fixed_rhs, data = grid_ctrl)[, names(beta_sp), drop = FALSE]
      X_diff <- X_case - X_ctrl
      diff_est <- as.vector(X_diff %*% beta_sp)
      diff_se <- sqrt(rowSums((X_diff %*% Vsp) * X_diff))
      sp_diff_all[[ch]][[model_id]] <- data.frame(
        Cohort = ch, model_id = model_id, model_label = spec$label,
        model_role = spec$role, Age_Centered = pred_age,
        diff = diff_est, se = diff_se,
        lwr = diff_est - 1.96 * diff_se,
        upr = diff_est + 1.96 * diff_se,
        vcov_type = Vobj$type
      )

      deriv_row <- function(t, h = 1e-4) {
        c_plus <- make_sp_grid(t + h, "Cancer Case")
        c_minus <- make_sp_grid(t - h, "Cancer Case")
        r_plus <- make_sp_grid(t + h, "Control")
        r_minus <- make_sp_grid(t - h, "Control")
        (model.matrix(sp_fixed_rhs, data = c_plus)[, names(beta_sp), drop = FALSE] -
           model.matrix(sp_fixed_rhs, data = c_minus)[, names(beta_sp), drop = FALSE] -
           model.matrix(sp_fixed_rhs, data = r_plus)[, names(beta_sp), drop = FALSE] +
           model.matrix(sp_fixed_rhs, data = r_minus)[, names(beta_sp), drop = FALSE]) /
          (2 * h)
      }
      deriv_X <- deriv_row(pred_age)
      deriv_est <- as.vector(deriv_X %*% beta_sp)
      deriv_se <- sqrt(rowSums((deriv_X %*% Vsp) * deriv_X))
      sp_deriv_all[[ch]][[model_id]] <- data.frame(
        Cohort = ch, model_id = model_id, model_label = spec$label,
        model_role = spec$role, Age_Centered = pred_age,
        derivative = deriv_est, se = deriv_se,
        lwr = deriv_est - 1.96 * deriv_se,
        upr = deriv_est + 1.96 * deriv_se,
        vcov_type = Vobj$type,
        infer_method = paste(Vobj$type, "normal approximation")
      )
      theta_supported <- prediction_window[[1]] <= -8 && prediction_window[[2]] >= 8
      difference_design <- function(t) {
        model.matrix(sp_fixed_rhs, make_sp_grid(t, "Cancer Case"))[, names(beta_sp), drop = FALSE] -
          model.matrix(sp_fixed_rhs, make_sp_grid(t, "Control"))[, names(beta_sp), drop = FALSE]
      }
      exact <- .exact_slope_contrasts(difference_design, 8)
      c_post <- exact$post; c_pre <- exact$pre; c_theta <- exact$theta
      theta <- if (theta_supported) as.numeric(c_theta %*% beta_sp) else NA_real_
      theta_se <- if (theta_supported) sqrt(as.numeric(t(c_theta) %*% Vsp %*% c_theta)) else NA_real_
      theta_test <- if (theta_supported) {
        wald_with_vcov(m_sp, Vobj, matrix(c_theta, nrow = 1), "theta")
      } else {
        data.frame(Test = "theta", Fstat = NA_real_, df_num = 1,
                   df_denom = NA_real_, p_value = NA_real_, infer_method = Vobj$type)
      }
      theta_df <- theta_test$df_denom[[1]]
      theta_crit <- if (is.na(theta_df)) NA_real_ else qt(0.975, df = theta_df)
      pre_theta <- if (theta_supported) as.numeric(c_pre %*% beta_sp) else NA_real_
      pre_theta_se <- if (theta_supported) sqrt(as.numeric(t(c_pre) %*% Vsp %*% c_pre)) else NA_real_
      endpoint_diff <- function(t) {
        xx <- model.matrix(sp_fixed_rhs, data = make_sp_grid(t, "Cancer Case"))[, names(beta_sp), drop = FALSE] -
          model.matrix(sp_fixed_rhs, data = make_sp_grid(t, "Control"))[, names(beta_sp), drop = FALSE]
        as.numeric(xx %*% beta_sp)
      }
      theta_closed_form <- if (theta_supported) {
        (endpoint_diff(8) - endpoint_diff(0)) / 8 -
          (endpoint_diff(0) - endpoint_diff(-8)) / 8
      } else NA_real_
      theta_difference <- theta - theta_closed_form
      if (theta_supported && (!is.finite(theta_difference) || abs(theta_difference) > 1e-6)) {
        stop("Exact contrast and endpoint theta differ by more than 1e-6 for ",
             ch, " / ", model_id, ": ", format(theta_difference, digits = 12), call. = FALSE)
      }
      pre_constraint <- matrix(c_pre, nrow = 1)
      pre_test <- if (theta_supported) {
        wald_with_vcov(m_sp, Vobj, pre_constraint, "pre-index average differential slope")
      } else data.frame(p_value = NA_real_, df_denom = NA_real_)
      pre_df <- pre_test$df_denom[[1]]
      pre_crit <- if (is.na(pre_df)) NA_real_ else qt(0.975, df = pre_df)
      interaction_terms <- .spline_interaction_terms(names(beta_sp), time_terms)
      omnibus_all[[paste(ch, model_id, sep = "__")]] <-
        joint_wald(m_sp, Vobj, interaction_terms,
                   "trajectory-shape difference (secondary)") %>%
        mutate(Cohort = ch, model_id = model_id, .before = 1)
      sp_theta_all[[ch]][[model_id]] <- data.frame(
        Cohort = ch, model_id = model_id, model_label = spec$label,
        model_role = spec$role, theta = theta, se = theta_se,
        df = theta_df, p_value = theta_test$p_value[[1]],
        lwr = theta - theta_crit * theta_se,
        upr = theta + theta_crit * theta_se,
        infer_method = theta_test$infer_method[[1]],
        inference_status = Vobj$inference_status,
        theta_closed_form = theta_closed_form,
        theta_closed_form_difference = theta_difference,
        pre_index_slope = pre_theta, pre_index_slope_se = pre_theta_se,
        pre_index_slope_lwr = pre_theta - pre_crit * pre_theta_se,
        pre_index_slope_upr = pre_theta + pre_crit * pre_theta_se,
        pre_index_slope_df = pre_df,
        pre_index_slope_p = pre_test$p_value[[1]],
        pre_index_slope_infer_method = if ("infer_method" %in% names(pre_test)) pre_test$infer_method[[1]] else Vobj$type,
        pre_window = "[-8, 0)", post_window = "(0, +8]",
        theta_supported = theta_supported,
        prediction_min_time = prediction_window[[1]],
        prediction_max_time = prediction_window[[2]],
        vcov_type = Vobj$type
      )

      # Deterministic grouped jackknife: sort matched sets, assign ten nearly
      # equal groups, and refit the primary-adjusted spline (M1) after dropping
      # each group in turn.
      if (run_influence && identical(model_id, "M1_primary_spline") && theta_supported) {
        set_ids <- sort(unique(as.character(d_model$matching_id)))
        decile_map <- setNames(pmin(10L, ceiling(seq_along(set_ids) * 10 / length(set_ids))),
                               set_ids)
        influence_all[[ch]] <- bind_rows(lapply(seq_len(10), function(decile) {
          d_jk <- d_model[decile_map[as.character(d_model$matching_id)] != decile, , drop = FALSE]
          jk <- tryCatch(.fit_lmer_with_ladder(model_formula, d_jk,
                                               paste0(spec$label, " influence decile ", decile)),
                         error = function(e) e)
          if (inherits(jk, "error")) {
            return(data.frame(Cohort = ch, model_id = model_id,
                              dropped_matching_id_decile = decile,
                              n_obs = nrow(d_jk), theta = NA_real_,
                              theta_change = NA_real_, pct_theta_change = NA_real_,
                              rung = NA_character_, status = conditionMessage(jk)))
          }
          theta_jk <- as.numeric(c_theta %*% fixef(jk$fit))
          data.frame(Cohort = ch, model_id = model_id,
                     dropped_matching_id_decile = decile,
                     n_obs = nrow(d_jk), theta = theta_jk,
                     theta_change = theta_jk - theta,
                     pct_theta_change = if (theta == 0) NA_real_ else 100 * (theta_jk - theta) / abs(theta),
                     rung = jk$rung, status = "fit")
        }))
      }

      sp_status_all[[ch]][[model_id]] <- data.frame(
        Cohort = ch, model_id = model_id, model_label = spec$label,
        model_role = spec$role, time_structure = spec$time_structure,
        spline_df = spec$spline_df,
        status = if (isTRUE(fit_info$boundary_matching_variance_zero)) {
          "fit_boundary_matching_variance_zero"
        } else "fit",
        n_obs = nrow(d_model), n_id = n_distinct(d_model$id),
        n_matching_id = n_distinct(d_model$matching_id),
          fitting_min_time = min(d_model$Age_Centered, na.rm = TRUE),
          fitting_max_time = max(d_model$Age_Centered, na.rm = TRUE),
          prediction_min_time = prediction_window[[1]],
          prediction_max_time = prediction_window[[2]],
        rung = fit_info$rung, singular = fit_info$singular,
        convergence = fit_info$converged,
        optimizer_code = fit_info$optimizer_code,
        scaled_gradient = fit_info$scaled_gradient,
        convergence_messages = fit_info$convergence_messages,
        vcov_type = Vobj$type,
        inference_status = Vobj$inference_status,
        matching_id_variance = fit_info$matching_id_variance,
        message = if (isTRUE(fit_info$boundary_matching_variance_zero)) {
          "matching_id intercept variance was on the permitted zero boundary"
        } else NA_character_
      )
      checkpoint()
    }
  }

  support   <- bind_rows(support_all)
  filter_log <- bind_rows(filter_all)
  sp_pred <- .bind_model_rows(unlist(sp_pred_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  sp_diff <- .bind_model_rows(unlist(sp_diff_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  sp_deriv <- .bind_model_rows(unlist(sp_deriv_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  sp_coef <- .bind_model_rows(unlist(sp_coef_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  sp_theta <- .bind_model_rows(unlist(sp_theta_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  sp_status <- .bind_model_rows(unlist(sp_status_all, recursive = FALSE)) %>%
    mutate(Cohort = factor(Cohort, levels = cohorts))
  convergence_log <- bind_rows(convergence_all)
  variance_components <- bind_rows(variance_all)
  residual_diagnostics <- bind_rows(residual_diag_all)
  random_effect_diagnostics <- bind_rows(random_effect_diag_all)
  bounded_diagnostics <- bind_rows(bounded_diag_all)
  spline_omnibus <- bind_rows(omnibus_all)
  influence_diagnostics <- if (length(influence_all)) bind_rows(influence_all) else
    data.frame(Cohort = cohorts, model_id = "M1_primary_spline",
               status = "not_run; use the separate diagnostic runner")
  spline_scaling <- bind_rows(scaling_all)
  sp_pred_primary <- sp_pred %>% filter(model_id == "M1_primary_spline")

  # Persist the actionable diagnostics when the primary M1 spline has no
  # predictions rather than
  # allowing ggplot's generic empty-facet error to hide the failed strict fit.
  if (!nrow(sp_pred_primary)) {
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
    write.csv(sp_status, file.path(results_dir, paste0(out_prefix, "_spline_model_status.csv")),
              row.names = FALSE)
    write.csv(convergence_log, file.path(results_dir, paste0(out_prefix, "_convergence_attempts.csv")),
              row.names = FALSE)
    write.csv(spline_scaling, file.path(results_dir, paste0(out_prefix, "_spline_scaling.csv")),
              row.names = FALSE)
    write.csv(model_specification,
              file.path(results_dir, paste0(out_prefix, "_model_specification.csv")),
              row.names = FALSE)
    write.csv(filter_log, file.path(results_dir, paste0(out_prefix, "_model_filter_log.csv")),
              row.names = FALSE)
    saveRDS(metadata_all, file.path(results_dir, paste0(out_prefix, "_spline_metadata.rds")))
    stop(
      "No M1 primary spline prediction rows were produced. Inspect its saved ",
      "model status and convergence attempts: ",
      file.path(results_dir, paste0(out_prefix, "_spline_model_status.csv")),
      " and ",
      file.path(results_dir, paste0(out_prefix, "_convergence_attempts.csv")),
      call. = FALSE
    )
  }

  p_spline <- ggplot(sp_pred_primary, aes(x = Age_Centered, y = pred, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    facet_wrap(~ Cohort) +
    scale_color_manual(values = group_cols) +
    scale_fill_manual(values = group_cols) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    labs(title = analysis_title,
         subtitle = "M1 primary natural-spline Gaussian LME; centered on own attained age at index",
         caption = paste("Pointwise normal-approximation intervals; covariance:",
                         paste(unique(sp_pred_primary$vcov_type), collapse = "; ")),
         x = "Years relative to index",
         y = "Predicted frailty index (no-cancer FI)", color = NULL, fill = NULL) +
    scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 4))

  # --------------------------------- outputs ----------------------------------
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  if (!dir.exists(visuals_dir)) dir.create(visuals_dir, recursive = TRUE)

  write.csv(support, file.path(results_dir, paste0(out_prefix, "_support_by_time_bin.csv")),
            row.names = FALSE)
  write.csv(filter_log, file.path(results_dir, paste0(out_prefix, "_model_filter_log.csv")),
            row.names = FALSE)
  write.csv(censoring_log, file.path(results_dir, paste0(out_prefix, "_own_cancer_censoring_log.csv")),
            row.names = FALSE)
  write.csv(sp_coef, file.path(results_dir, paste0(out_prefix, "_spline_fixed_effects_CI.csv")),
            row.names = FALSE)
  write.csv(sp_pred, file.path(results_dir, paste0(out_prefix, "_spline_predicted_trajectories.csv")),
            row.names = FALSE)
  write.csv(sp_diff, file.path(results_dir, paste0(out_prefix, "_spline_group_difference.csv")),
            row.names = FALSE)
  write.csv(sp_deriv, file.path(results_dir, paste0(out_prefix, "_spline_derivative.csv")),
            row.names = FALSE)
  write.csv(sp_theta, file.path(results_dir, paste0(out_prefix, "_spline_theta.csv")),
            row.names = FALSE)
  write.csv(sp_status, file.path(results_dir, paste0(out_prefix, "_spline_model_status.csv")),
            row.names = FALSE)
  write.csv(convergence_log, file.path(results_dir, paste0(out_prefix, "_convergence_attempts.csv")),
            row.names = FALSE)
  write.csv(variance_components, file.path(results_dir, paste0(out_prefix, "_spline_variance_components.csv")),
            row.names = FALSE)
  write.csv(residual_diagnostics, file.path(results_dir, paste0(out_prefix, "_residual_diagnostics.csv")),
            row.names = FALSE)
  write.csv(random_effect_diagnostics,
            file.path(results_dir, paste0(out_prefix, "_random_effect_diagnostics.csv")), row.names = FALSE)
  write.csv(bounded_diagnostics, file.path(results_dir, paste0(out_prefix, "_spline_bounded_prediction_checks.csv")),
            row.names = FALSE)
  write.csv(spline_omnibus, file.path(results_dir, paste0(out_prefix, "_spline_omnibus_tests.csv")),
            row.names = FALSE)
  write.csv(influence_diagnostics, file.path(results_dir, paste0(out_prefix, "_spline_influence_refits.csv")),
            row.names = FALSE)
  write.csv(spline_scaling, file.path(results_dir, paste0(out_prefix, "_spline_scaling.csv")),
            row.names = FALSE)
  write.csv(model_specification,
            file.path(results_dir, paste0(out_prefix, "_model_specification.csv")),
            row.names = FALSE)
  saveRDS(metadata_all, file.path(results_dir, paste0(out_prefix, "_spline_metadata.rds")))
  if (identical(out_prefix, "4.5")) {
    saveRDS(metadata_all[["All Cancer Cohort"]],
            file.path(results_dir, "glme_spline_metadata.rds"))
  }

  cat("\nSaved spline outputs (prefix '", out_prefix, "') to: ",
      results_dir, "\n", sep = "")


  invisible(list(
    sp_contexts = sp_contexts, sp_models = sp_models, sp_coef = sp_coef, sp_pred = sp_pred,
    sp_diff = sp_diff, sp_deriv = sp_deriv, sp_theta = sp_theta, sp_status = sp_status,
    convergence_log = convergence_log, variance_components = variance_components,
    residual_diagnostics = residual_diagnostics,
    random_effect_diagnostics = random_effect_diagnostics,
    bounded_diagnostics = bounded_diagnostics, spline_omnibus = spline_omnibus,
    influence_diagnostics = influence_diagnostics, spline_scaling = spline_scaling,
    model_specification = model_specification,
    matching_provenance = matching_provenance,
    figures = list(spline = p_spline)
  ))
}

render_overall_glme_html_report <- function(results_dir, sensitivity_dir, visuals_dir,
                                            report_rmd = NULL) {
  if (is.null(report_rmd) || !file.exists(report_rmd)) {
    warning("Overall GLME HTML report skipped because the active R Markdown file is unavailable: ",
            report_rmd, call. = FALSE)
    return(invisible(FALSE))
  }
  if (!dir.exists(visuals_dir)) dir.create(visuals_dir, recursive = TRUE)
  output_html <- file.path(visuals_dir, "4.5_GLME_overall.html")
  ok <- tryCatch({
    if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
      rmarkdown::render(report_rmd, output_file = basename(output_html),
                        output_dir = visuals_dir, quiet = TRUE,
                        envir = new.env(parent = globalenv()))
    } else {
      if (!requireNamespace("knitr", quietly = TRUE) ||
          !requireNamespace("markdown", quietly = TRUE)) {
        stop("Rendering requires either rmarkdown with Pandoc or both knitr and markdown.")
      }
      knitted_md <- tempfile(fileext = ".md")
      knitr::knit(report_rmd, output = knitted_md, quiet = TRUE,
                  envir = new.env(parent = globalenv()))
      markdown::markdownToHTML(
        knitted_md,
        output = output_html,
        options = c("+embed_resources", "+toc", "+table", "+auto_identifiers"),
        title = "Overall Incident Cancer Risk-Set Matching and GLME"
      )
    }
    if (!file.exists(output_html) || file.info(output_html)$size <= 0) {
      stop("The report renderer did not create a nonempty HTML file.")
    }
    TRUE
  }, error = function(e) {
    warning("Overall GLME HTML report render failed: ", conditionMessage(e), call. = FALSE)
    FALSE
  })
  invisible(ok)
}
