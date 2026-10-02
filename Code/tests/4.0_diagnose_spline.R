# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-20
# Purpose: Diagnose participant reuse, random-slope numerics, matching-set
# boundaries and covariance failures without changing production results.
# Default --prepare mode reads data and writes descriptive diagnostics only.
# Nemo must run --fit to perform GLME refits. --influence additionally runs
# participant-deletion diagnostics for the five most reused participants.
# All diagnostic files and resumable fits are saved under Codex/glme_diagnostics.
# Usage from project root:
#   Rscript Code/tests/4.0_diagnose_spline.R --prepare
#   Rscript Code/tests/4.0_diagnose_spline.R --fit
#   Rscript Code/tests/4.0_diagnose_spline.R --fit --influence
# Final-plan checks: --fit --plan --influence (Nemo runs; includes refits).
# Optional: --dataset=riskset_matched_analysis_long.rds (runs each cohort).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(getwd())
engine <- file.path(project_dir, "Code/2_data_analysis/4.0_GLME_spline_functions.R")
if (!file.exists(engine)) stop("Run this script from the Frailty HPFS project root")
source(engine)
dataset_arg <- grep("^--dataset=", args, value = TRUE)
dataset <- if (length(dataset_arg)) sub("^--dataset=", "", dataset_arg[[1]]) else "riskset_matched_overall_long.rds"
if (basename(dataset) != dataset) stop("--dataset must be a filename under Data")
matched_path <- file.path(project_dir, "Data", dataset)
provenance <- validate_matching_provenance(matched_path)
engine_md5 <- unname(tools::md5sum(engine))
runner_md5 <- unname(tools::md5sum(file.path(project_dir, "Code/tests/4.0_diagnose_spline.R")))
root <- file.path(project_dir, "Codex", "glme_diagnostics", tools::file_path_sans_ext(dataset),
                  paste0(substr(provenance$input_md5, 1, 8), "_", substr(engine_md5, 1, 8), "_", substr(runner_md5, 1, 8)))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
if (requireNamespace("clubSandwich", quietly = TRUE)) {
  method <- getS3method("vcovCR", "lmerMod", envir = asNamespace("clubSandwich"))
  writeLines(capture.output(print(method)), file.path(root, "installed_vcovCR_lmerMod.txt"))
}
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
saveRDS(list(input = provenance, engine_md5 = engine_md5, runner_md5 = runner_md5, arguments = args), file.path(root, "provenance.rds"))
d <- readRDS(matched_path)
if (!"matching_id" %in% names(d)) d$matching_id <- paste(d$Cohort, d$match_set, sep = "__")
if (!"trajectory_id" %in% names(d)) d$trajectory_id <- paste(d$Cohort, d$match_set, d$id, d$role, sep = "__")
required <- c(.primary_covars, "fi_score_nocancer", "Age_Centered", "id", "cycle", "Group", "post_own_cancer")
if (!all(required %in% names(d))) stop("Missing diagnostic input columns")
d <- d %>% filter(!post_own_cancer, !is.na(fi_score_nocancer), !is.na(Age_Centered),
                  if_all(all_of(.primary_covars), ~ !is.na(.x))) %>%
  mutate(Group = factor(Group, levels = c("Control", "Cancer Case")),
         across(all_of(c("id", "matching_id", "base_race", "base_marital", "base_living")), factor))
for (cohort in unique(d$Cohort)) {
  out <- file.path(root, gsub("[^A-Za-z0-9]+", "_", cohort))
  dir.create(out, showWarnings = FALSE)
  dat <- droplevels(d[d$Cohort == cohort, ])
  if ("--plan" %in% args) invisible(.set_random_clock(dat))
  reuse <- dat %>% group_by(id) %>% summarize(
    n_rows = n(), n_assignments = n_distinct(trajectory_id), n_cycles = n_distinct(cycle),
    .groups = "drop") %>% arrange(desc(n_assignments), id)
  repeated <- dat %>% group_by(id, cycle) %>% summarize(
    n_copies = n(), n_relative_times = n_distinct(Age_Centered),
    relative_time_range = max(Age_Centered) - min(Age_Centered), .groups = "drop")
  write.csv(reuse, file.path(out, "participant_reuse.csv"), row.names = FALSE)
  write.csv(data.frame(cohort = cohort, n_rows = nrow(dat), n_participants = nrow(reuse),
    n_reused_participants = sum(reuse$n_assignments > 1),
    max_assignments_per_participant = max(reuse$n_assignments),
    n_unique_observations = nrow(repeated),
    n_observations_with_multiple_relative_times = sum(repeated$n_relative_times > 1),
    maximum_index_shift_years = max(repeated$relative_time_range)),
    file.path(out, "preparation_summary.csv"), row.names = FALSE)
  support <- add_relative_time_bin(dat) %>% group_by(rel_time_bin, .drop = FALSE) %>%
    summarize(n_case = n_distinct(id[Group == "Cancer Case"]),
              n_ctrl = n_distinct(id[Group == "Control"]), .groups = "drop") %>%
    mutate(support_ok = n_case >= 50 & n_ctrl >= 250)
  write.csv(support, file.path(out, "support.csv"), row.names = FALSE)
  if (!"--fit" %in% args) next

  spec <- .trajectory_model_specs$M1_primary_spline
  bases <- list(adjusted_spline = .make_spline_basis(dat$Age_Centered, 3, "S"))
  reference <- .make_ref_grid(dat, covars = spec$covars, .reference = 1)
  scaled <- .scale_model_numeric_columns(.add_model_time_terms(dat, spec, bases), bases$adjusted_spline$terms)
  rhs <- .fixed_rhs(bases$adjusted_spline$terms, spec$covars)
  configurations <- list(
    M1_intercept = list(random = "(1 | id)", scale = 4),
    M1_correlated_unscaled = list(random = "(1 + Age_Centered | id)", scale = 1),
    M1_correlated_scaled = list(random = "(1 + Age_Centered | id)", scale = 4),
    M1_uncorrelated_scaled = list(random = "(1 + Age_Centered || id)", scale = 4),
    M3_correlated_scaled = list(random = "(1 + Age_Centered | id) + (1 | matching_id)", scale = 4),
    M3_intercept = list(random = "(1 | id) + (1 | matching_id)", scale = 4))
  if ("age_at_cycle" %in% names(dat) && all(is.finite(dat$age_at_cycle)))
    configurations$M1_attained_age <- list(random = "(1 + Age_Centered | id)", scale = 4, attained = TRUE)
  if ("--plan" %in% args) configurations <- list(
    M1_attained_age = list(random = "(1 + .random_time | id)", scale = 4, attained = TRUE),
    M1_attained_age_uncorrelated = list(random = "(1 + .random_time || id)", scale = 4, attained = TRUE),
    M1_attained_age_intercept = list(random = "(1 | id)", scale = 4, attained = TRUE),
    M3_attained_age = list(random = "(1 + .random_time | id) + (1 | matching_id)", scale = 4, attained = TRUE))
  if ("--influence" %in% args && !"--plan" %in% args) for (i in seq_len(min(5L, nrow(reuse)))) {
    configurations[[paste0("M1_drop_reused_", i)]] <- list(
      random = "(1 + Age_Centered | id)", scale = 4, drop_id = as.character(reuse$id[[i]]))
  }
  summary <- list()
  config_index <- 0L
  while (config_index < length(configurations)) {
    config_index <- config_index + 1L
    name <- names(configurations)[[config_index]]
    cfg <- configurations[[name]]
    model_data <- scaled$data
    if (!is.null(cfg$drop_id)) model_data <- model_data[!as.character(model_data$id) %in% cfg$drop_id, ]
    # Fixed splines and original relative time remain unchanged for every clock.
    form <- as.formula(paste("fi_score_nocancer", paste(deparse(rhs), collapse = ""), "+", cfg$random))
    cache_path <- file.path(out, paste0(name, "_fit.rds"))
    info <- if (file.exists(cache_path)) readRDS(cache_path) else tryCatch(
      .fit_lmer_with_ladder(form, model_data, name, simplify = FALSE, numerical_retries = TRUE,
        random_time_scale = cfg$scale, random_clock = if (isTRUE(cfg$attained)) "attained_age" else "relative_time",
        allow_matching_boundary = startsWith(name, "M3_"),
        diagnostic_path = file.path(out, paste0(name, "_attempts.csv"))), error = function(e) e)
    if (inherits(info, "error")) {
      summary[[name]] <- data.frame(configuration = name, status = "failed", message = conditionMessage(info))
    } else {
      saveRDS(info, cache_path)
      fit <- info$fit
      grids <- .spline_grid_factory(reference, spec, bases, scaled$parameters, names(fixef(fit)))
      window <- continuous_support_window(support, 20)
      theta_supported <- all(is.finite(window)) && window[1] <= -8 && window[2] >= 8
      C <- .inference_constraints(names(fixef(fit)), bases$adjusted_spline$terms,
                                  grids$difference, theta_supported = theta_supported)
      cov <- get_primary_vcov(fit, model.frame(fit)$id, constraints = C,
                              diagnostic_path = file.path(out, paste0(name, "_covariance.csv")))
      exact <- .exact_slope_contrasts(grids$difference, 8)$theta
      wt <- if (theta_supported) wald_with_vcov(fit, cov, matrix(exact, nrow = 1), "theta") else NULL
      summary[[name]] <- data.frame(configuration = name, status = "fit", message = NA_character_,
        n_rows = nobs(fit), logLik = as.numeric(logLik(fit)), AIC = AIC(fit),
        matching_id_variance = info$matching_id_variance,
        matching_boundary = info$boundary_matching_variance_zero,
        theta = if (theta_supported) as.numeric(exact %*% fixef(fit)) else NA_real_,
        theta_se = if (theta_supported) sqrt(as.numeric(exact %*% cov$V %*% exact)) else NA_real_,
        theta_p = if (theta_supported) wt$p_value else NA_real_,
        theta_lwr = if (theta_supported) as.numeric(exact %*% fixef(fit)) - qt(.975, wt$df_denom) * sqrt(as.numeric(exact %*% cov$V %*% exact)) else NA_real_,
        theta_upr = if (theta_supported) as.numeric(exact %*% fixef(fit)) + qt(.975, wt$df_denom) * sqrt(as.numeric(exact %*% cov$V %*% exact)) else NA_real_,
        infer_method = if (theta_supported) wt$infer_method else NA_character_,
        covariance = cov$type, inference_status = cov$inference_status,
        dropped_id = if (is.null(cfg$drop_id)) NA_character_ else paste(cfg$drop_id, collapse = ";"))
      write.csv(coef_table_with_vcov(fit, cov), file.path(out, paste0(name, "_coefficients.csv")), row.names = FALSE)
      if (name == "M1_attained_age" && all(c("--plan", "--influence") %in% args)) {
        # Rank participants using this fit, then remove all copies of their rows.
        # Keep the original fixed-effect basis/scaling and random structure.
        E <- .participant_vcov(fit, model.frame(fit)$id, "CR2", "estfun")
        stopifnot(cov$robust, startsWith(cov$type, "CR2"),
          isTRUE(all.equal(unname(tcrossprod(E)), unname(cov$V), tolerance = 1e-8)))
        score <- as.vector(exact %*% E)
        ranking <- data.frame(id = levels(droplevels(model.frame(fit)$id)), contribution = score^2)
        ranking <- ranking[order(ranking$contribution, decreasing = TRUE), ]
        ranking$variance_share <- ranking$contribution / sum(ranking$contribution)
        write.csv(ranking, file.path(out, "M1_influence_ranking.csv"), row.names = FALSE)
        for (j in seq_len(min(5L, nrow(ranking)))) configurations[[paste0("M1_drop_influential_", j)]] <-
          list(random = "(1 + .random_time | id)", scale = 4, attained = TRUE, drop_id = ranking$id[j])
        configurations$M1_drop_top_1pct <- list(random = "(1 + .random_time | id)", scale = 4,
          attained = TRUE, drop_id = head(ranking$id, ceiling(.01 * nrow(ranking))))
      }
      write.csv(.variance_components(fit, cohort, name), file.path(out, paste0(name, "_variance.csv")), row.names = FALSE)
    }
    write.csv(bind_rows(summary), file.path(out, "fit_comparison.csv"), row.names = FALSE)
  }
}
cat("Diagnostics saved under:", root, "\n")
if (!"--fit" %in% args) cat("Preparation only: no GLME models fitted. Nemo can run this script with --fit.\n")
