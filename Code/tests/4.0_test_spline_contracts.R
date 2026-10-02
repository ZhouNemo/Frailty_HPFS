# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
# Purpose: Fast regression checks for scaling, support gaps, exact contrasts,
# covariance/test separation, attained-age clock consistency and failure-safe orchestration. Never fits any GLME:
# the fitting function is replaced by an error stub for orchestration tests.
# Run from project root: Rscript Code/tests/4.0_test_spline_contracts.R
# Temporary fixture data/outputs are removed on exit; no production data writes.
# =============================================================================
source("Code/2_data_analysis/4.0_GLME_spline_functions.R")
expect_error <- function(expr, pattern = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(e, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(e), fixed = TRUE))
}
# An entirely absent bin must break the continuous interval.
support <- data.frame(rel_time_bin = c("-2 to 0", "2 to +4"), support_ok = TRUE)
stopifnot(identical(continuous_support_window(support, 20), c(-2, 0)))
support <- rbind(support, data.frame(rel_time_bin = "0 to +2", support_ok = TRUE))
stopifnot(identical(continuous_support_window(support, 20), c(-2, 4)))
stopifnot(all(is.na(continuous_support_window(support[-1, ], 20))))
# Retain df-3 basis construction and apply raw dietary means exactly once.
x <- data.frame(id = 1:30, Age_Centered = seq(-20, 20, length.out = 30),
  base_calor = seq(1800, 2200, length.out = 30))
spec <- list(covars = "base_calor", basis_key = "adjusted_spline")
bases <- list(adjusted_spline = .make_spline_basis(x$Age_Centered, 3, "S"))
stopifnot(ncol(bases$adjusted_spline$basis) == 3)
ref <- .make_ref_grid(x, covars = spec$covars, .reference = 1)
scaled <- .scale_model_numeric_columns(.add_model_time_terms(x, spec, bases), c("S1", "S2", "S3", "base_calor"))
g <- .apply_model_scaling(.add_model_time_terms(transform(ref, Age_Centered = 0), spec, bases), scaled$parameters)
stopifnot(abs(g$base_calor) < 1e-12)
# Exact average-slope contrasts for a quadratic difference curve.
D <- function(t) cbind(1, t, t^2)
C <- .exact_slope_contrasts(D, 8)
stopifnot(isTRUE(all.equal(unname(C$theta), c(0, 0, 16))))
# Covariance guards distinguish valid, misaligned and indefinite matrices.
V <- diag(2); dimnames(V) <- list(c("a", "b"), c("a", "b"))
stopifnot(.validate_covariance(V, c("a", "b")))
expect_error(.validate_covariance(V, c("b", "a")), "names/order")
V[1, 2] <- V[2, 1] <- 2
expect_error(.validate_covariance(V, c("a", "b")), "semidefinite")
expect_error(.validate_wald(data.frame(Fstat = 1, df_num = 1, df_denom = NA_real_, p_val = .5)))
stopifnot(identical(names(.bind_model_rows(list())), c("Cohort", "model_id")))
# Candidate construction is testable without calling lmer.
candidates <- .lmer_attempts(y ~ S1 + (1 + .random_time | id))
stopifnot(length(candidates) == 3L,
  grepl(".random_time || id", paste(deparse(candidates[[2]]$formula), collapse = ""), fixed = TRUE),
  candidates[[3]]$rung == "random intercept only")
stopifnot(length(.lmer_attempts(y ~ S1 + (1 | id))) == 1L,
  length(.lmer_attempts(y ~ S1 + (1 + .random_time | id), TRUE, FALSE)) == 2L)
# Exercise logged covariance failure with an S3 fixture, never a fitted model.
fixef.review_fixture <- function(object, ...) c(a = 1, b = 2)
model.frame.review_fixture <- function(formula, ...) data.frame(id = factor(1:4))
nobs.review_fixture <- function(object, ...) 4L
vcov.review_fixture <- function(object, ...) {
  m <- diag(2); dimnames(m) <- list(c("a", "b"), c("a", "b")); m
}
fixture <- structure(list(), class = "review_fixture")
log_path <- tempfile(fileext = ".csv")
v <- suppressWarnings(get_primary_vcov(fixture, factor(1:4), diagnostic_path = log_path))
log <- read.csv(log_path)
stopifnot(!v$robust, v$inference_status == "model_based_fallback_requires_review",
  any(log$stage == "vcovCR" & log$status == "error"),
  any(log$status == "fallback_requires_review"))
unlink(log_path)
# Successful covariance must survive optional coefficient/HTZ failure.
original_covariance <- .participant_vcov
selected_types <- character()
.participant_vcov <- function(model, cluster, type, ...) {
  selected_types <<- c(selected_types, type)
  vcov(model)
}
v <- get_primary_vcov(fixture, factor(1:4))
stopifnot(v$robust, identical(selected_types, "CR2"),
          v$inference_status == "robust_asymptotic",
          all(is.finite(coef_table_with_vcov(fixture, v)$CI_low)))
v$small_sample <- TRUE # force a downstream optional-test error with the fixture
w <- wald_with_vcov(fixture, v, matrix(c(1, 0), nrow = 1), "theta")
stopifnot(is.finite(w$p_value), is.infinite(w$df_denom), grepl("CR2", w$infer_method))
.participant_vcov <- original_covariance
expect_error(get_primary_vcov(fixture, factor(4:1)), "Cluster vector")
# Clock consistency includes every copy of a participant-cycle and preserves
# fixed diagnosis-relative time. Missing/inconsistent attained age fails closed.
x <- data.frame(id = c(1, 1, 1), cycle = c(2000, 2000, 2004),
                age_at_cycle = c(60, 60, 64), Age_Centered = c(-4, 0, 4))
y <- .set_random_clock(x)
stopifnot(identical(y$Age_Centered, x$Age_Centered), identical(y$.random_time, c(0, 0, 1)))
x$age_at_cycle[2] <- 61
expect_error(.set_random_clock(x), "differs across copies")
expect_error(.set_random_clock(x[, names(x) != "age_at_cycle"]), "required")
old_budget <- getOption("hpfs.small_sample_max_bytes")
options(hpfs.small_sample_max_bytes = 0)
stopifnot(!.small_sample_feasible(fixture, factor(1:4)))
options(hpfs.small_sample_max_bytes = old_budget)
# Missing M2 variables must skip M2; all other fit failures must still save logs.
run_fixture <- function(success = FALSE) {
  root <- tempfile("glme_contract_")
  dir.create(file.path(root, "Data"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  diagnostic_dir <- file.path(root, "Results/cancer/data/matching_diagnostics")
  dir.create(diagnostic_dir, recursive = TRUE)
  d <- expand.grid(id = 1:12, cycle = as.character(seq(1986, 2018, by = 4)))
  d$Age_Centered <- rep(seq(-8, 8, by = 2), each = 12)
  d$Cohort <- "fixture"; d$Group <- ifelse(d$id <= 6, "Cancer Case", "Control")
  d$role <- ifelse(d$id <= 6, "Case", "Control")
  d$match_set <- (d$id - 1) %% 6 + 1
  d$index_age_z <- (d$id - 6) / 3
  d$base_race <- d$base_marital <- d$base_living <- as.character(d$id %% 2)
  d$fi_score_nocancer <- .1 + d$id / 100
  d$post_own_cancer <- FALSE
  if (success) {
    d$base_pckgr <- as.character(d$id %% 3)
    d$base_calor <- 2000 + 20 * d$id
    d$base_sat <- 15 + d$id
    d$base_diet_chol <- 100 + d$id
    d$base_alco <- d$id / 2
  }
  path <- file.path(root, "Data", "fixture_long.rds"); saveRDS(d, path)
  write.csv(data.frame(gate_pass = TRUE), file.path(diagnostic_dir, "fixture_long_gate_g4.csv"), row.names = FALSE)
  ledger <- unique(d[c("id","Cohort","Group","role","match_set")])
  extra <- ledger[c(1,7),]; extra$id <- c(13,14); extra$match_set <- 7
  ledger <- rbind(ledger,extra)
  ledger$trajectory_id <- with(ledger,paste(Cohort,match_set,id,role,sep="__"))
  ap <- file.path(root,"Data","fixture_assignments.rds"); saveRDS(ledger,ap)
  saveRDS(list(eligibility_version = "no_fi_requirement_v1", fi_requirement = "none",
    cohort_entry_rule = "first_participated_analytic_cycle_return", control_entry_on_or_before_index = TRUE,
    assignment_md5=unname(tools::md5sum(ap)), output_md5 = unname(tools::md5sum(path))), file.path(diagnostic_dir, "fixture_long_run_metadata.rds"))
  original <- .fit_lmer_with_ladder
  calls <- 0L
  original_diagnostics <- .model_diagnostics
  original_variance <- .variance_components
  if (success) {
    assign(".model_diagnostics", function(...) list(residual = data.frame(check = "stub"),
      random_effect = data.frame(check = "stub")), envir = .GlobalEnv)
    assign(".variance_components", function(...) data.frame(check = "stub"), envir = .GlobalEnv)
    on.exit(assign(".model_diagnostics", original_diagnostics, envir = .GlobalEnv), add = TRUE)
    on.exit(assign(".variance_components", original_variance, envir = .GlobalEnv), add = TRUE)
  }
  assign(".fit_lmer_with_ladder", function(formula, data, ...) {
    calls <<- calls + 1L
    if (success) {
      X <- model.matrix(lme4::nobars(formula), data)
      beta <- setNames(rep(.001, ncol(X)), colnames(X)); beta[1] <- .2
      fit <- structure(list(beta = beta, data = data), class = "engine_fixture")
      return(list(fit = fit, converged = TRUE, rung = "stub", singular = FALSE,
        optimizer_code = 0L, scaled_gradient = 0, convergence_messages = NA_character_,
        matching_id_variance = NA_real_, boundary_matching_variance_zero = FALSE,
        attempts = data.frame(rung = "stub", accepted = TRUE)))
    }
    e <- simpleError("Intentional fit stub: no model executed")
    attr(e, "attempts") <- data.frame(rung = "stub", accepted = FALSE)
    stop(e)
  }, envir = .GlobalEnv)
  on.exit(assign(".fit_lmer_with_ladder", original, envir = .GlobalEnv), add = TRUE)
  results <- file.path(root, "Results/cancer/data")
  run <- function() suppressWarnings(run_spline_analysis(path, results, file.path(root, "visuals"),
    "test", "Fixture", "none", min_case_bin = 1, min_ctrl_bin = 1))
  if (success) {
    ans <- run()
    stopifnot(calls == 4L, nrow(ans$sp_pred) > 0, all(is.finite(ans$sp_theta$theta)),
      all(ans$sp_status$inference_status == "model_based_fallback_requires_review"))
    stopifnot(all(ans$sp_theta$theta_closed_form_difference < 1e-10))
  } else expect_error(run(), "No M1 primary")
  flow <- read.csv(file.path(results,"test_assignment_flow.csv"))
  stopifnot(all(flow$n_matched_assignments==7), all(flow$n_without_usable_fi==1),
    all(flow$n_fi_contributing_assignments==6))
  status <- read.csv(file.path(results, "test_spline_model_status.csv"))
  stopifnot(nrow(status) == 4)
  if (!success) stopifnot(calls == 3L,
    status$status[status$model_id == "M2_full_spline"] == "skipped_missing_columns")
  stopifnot(file.exists(file.path(results, "test_convergence_attempts.csv")))
  stopifnot(!any(grepl("eventstudy", list.files(results))))
}
fixef.engine_fixture <- function(object, ...) object$beta
model.frame.engine_fixture <- function(formula, ...) formula$data
nobs.engine_fixture <- function(object, ...) nrow(object$data)
vcov.engine_fixture <- function(object, ...) {
  m <- diag(1e-6, length(object$beta))
  dimnames(m) <- list(names(object$beta), names(object$beta)); m
}
run_fixture()
run_fixture(success = TRUE)
# Optional saved-model fixtures exercise structural gating without any refit.
# The fit-call boundary is replaced throughout; only cached objects are returned.
check_saved_ladder <- function() {
  cache <- "Codex/glme_diagnostics/riskset_matched_overall_long/364e02ea_0cf6c5ee_f1fa5b8b/All_Cancer_Cohort"
  path <- file.path(cache, "M1_attained_age_fit.rds")
  if (!file.exists(path)) return(invisible(NULL))
  good <- readRDS(path)$fit
  old <- .fit_lmer_candidate
  on.exit(assign(".fit_lmer_candidate", old, .GlobalEnv), add = TRUE)
  d <- data.frame(id = 1:3, cycle = 2000, age_at_cycle = 60:62, Age_Centered = -1:1)
  calls <- 0L
  assign(".fit_lmer_candidate", function(...) {
    calls <<- calls + 1L
    warning("Some predictor variables are on very different scales: consider rescaling")
    good
  }, .GlobalEnv)
  ans <- .fit_lmer_with_ladder(y ~ (1 + .random_time | id), d, "scaling fixture")
  stopifnot(calls == 1L, ans$rung == "correlated random slope")
  bad <- good
  bad@optinfo$conv$opt <- 1L
  assign(".fit_lmer_candidate", function(...) { calls <<- calls + 1L; bad }, .GlobalEnv)
  calls <- 0L
  expect_error(.fit_lmer_with_ladder(y ~ (1 + .random_time | id), d, "numerical fixture"), "Could not fit")
  stopifnot(calls == 2L) # Numerical failure alone cannot drop the slope.
  boundary <- good
  boundary@theta[] <- c(1, 0, 0)
  boundary@optinfo$derivs <- list(gradient = rep(0, 3), Hessian = diag(3))
  uncorrelated <- readRDS(file.path(cache, "M1_uncorrelated_scaled_fit.rds"))$fit
  calls <- 0L
  assign(".fit_lmer_candidate", function(...) {
    calls <<- calls + 1L
    if (calls <= 2L) boundary else uncorrelated
  }, .GlobalEnv)
  ans <- .fit_lmer_with_ladder(y ~ (1 + .random_time | id), d, "boundary fixture")
  stopifnot(calls == 3L, ans$rung == "uncorrelated random slope")
}
check_saved_ladder()
# Sensitivity selection is isolated from the formal engine; S7 is retired.
stopifnot(!exists("run_overall_glme_sensitivities", inherits = FALSE))
source("Code/2_data_analysis/4.0.1_GLME_sensitivity_analyses.R")
check_sensitivity_selection <- function() {
  root <- tempfile("glme_selection_"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "base.rds")
  saveRDS(data.frame(Cohort = "All Cancer Cohort", match_set = 1, id = 1,
                     role = "Case", Group = "Cancer Case"), path)
  old_provenance <- validate_matching_provenance
  old_fit <- .fit_one_spline_sensitivity
  old_summary <- .summarize_cached_sensitivity
  on.exit({
    assign("validate_matching_provenance", old_provenance, .GlobalEnv)
    assign(".fit_one_spline_sensitivity", old_fit, .GlobalEnv)
    assign(".summarize_cached_sensitivity", old_summary, .GlobalEnv)
  }, add = TRUE)
  provenance <- list(input_md5 = unname(tools::md5sum(path)))
  assign("validate_matching_provenance", function(...) provenance, .GlobalEnv)
  calls <- character()
  stub <- function(spec) list(status = data.frame(spec = spec$id, status = "fit", convergence = TRUE))
  assign(".fit_one_spline_sensitivity", function(d, spec, ...) {
    calls <<- c(calls, spec$id); stub(spec)
  }, .GlobalEnv)
  ans <- run_overall_glme_sensitivities(path, root, selected = "S1")
  stopifnot(identical(names(ans), "S1"), identical(calls, "S1"))
  expect_error(run_overall_glme_sensitivities(path, root, selected = "S3"), "require")
  assign(".summarize_cached_sensitivity", function(cache, spec) stub(spec), .GlobalEnv)
  ans <- run_overall_glme_sensitivities(path, root,
    primary_result = list(matching_provenance = provenance,
      sp_contexts = list("All Cancer Cohort" = list())), selected = c("S3", "S5"))
  stopifnot(identical(names(ans), c("S3", "S5")), identical(calls, "S1"))
  expect_error(run_overall_glme_sensitivities(path, root, selected = "bad"), "Select")
  expect_error(run_overall_glme_sensitivities(path, root, selected = "S7"), "retired")
}
check_sensitivity_selection()
cat("PASS: spline contracts; no matching or GLME fitting executed.\n")
