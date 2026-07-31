# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.2.2_GLME_M2_time_varying_other_cohorts.R
# Author:  Nemo Zhou
# Date started:      2026-07-28
# Date last updated: 2026-07-28
#
# Purpose:
#   Fit the same restricted M2 natural-spline GLME as
#   4.2.2_GLME_M2_time_varying.R for the other active cancer cohorts:
#   Low/Moderate Burden, High Burden, Obesity-Related Cancer, and All Cancer.
#   Each cohort is fit separately; the burden dataset supplies two cohort
#   levels and the other datasets supply one level each.
#
#   Model:
#     fi_score_nocancer ~ Group * natural spline(Age_Centered, df = 4)
#       + index_age_z + cycle-level race, marital status, living arrangement,
#         pack-year group, smoking status, pack-years, BMI, and physical activity
#       + (1 + Age_Centered | id)
#
#   Cycle-level covariates are joined from the canonical longitudinal panel by
#   participant id and questionnaire cycle. Current matched analytic rows have
#   no nonmissing cycle-level calories, saturated fat, dietary cholesterol, or
#   alcohol values, so those nutrition fields are omitted rather than imputed
#   or carried forward. `index_age_z` is retained as the index-age adjustment.
#
#   This script intentionally fits no event study, M0, M1, M3, sensitivity,
#   derivative, theta, or CR covariance analysis. Standard errors and 95% CIs
#   use the model-based covariance matrix from `vcov(fit)` only.
#
# Inputs:
#   Data/riskset_matched_analysis_long.rds
#   Data/riskset_matched_obesity_long.rds
#   Data/riskset_matched_overall_long.rds
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Results/cancer/data/matching_diagnostics/
#
# Outputs:
#   Results/cancer/data/4.2.2_other_m2_time_varying_fixed_effects_CI.csv
#   Results/cancer/data/4.2.2_other_m2_time_varying_predicted_trajectories.csv
#   Results/cancer/data/4.2.2_other_m2_time_varying_models.rds
# =============================================================================

library(dplyr)
library(lme4)
library(splines)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

panel_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
out_prefix <- "4.2.2_other_m2_time_varying"
spline_df <- 4L
prediction_window <- c(-20, 20)

cohort_specs <- list(
  burden = file.path(data_dir, "riskset_matched_analysis_long.rds"),
  obesity = file.path(data_dir, "riskset_matched_obesity_long.rds"),
  overall = file.path(data_dir, "riskset_matched_overall_long.rds")
)

if (!file.exists(panel_path)) {
  stop("Canonical longitudinal panel not found: ", panel_path, call. = FALSE)
}
missing_inputs <- names(cohort_specs)[!file.exists(unlist(cohort_specs))]
if (length(missing_inputs)) {
  stop("Required matched cohort files are missing: ",
       paste(missing_inputs, collapse = ", "), call. = FALSE)
}
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

validate_matching_provenance <- function(matched_path) {
  matching_stem <- tools::file_path_sans_ext(basename(matched_path))
  diagnostics_dir <- file.path(results_dir, "matching_diagnostics")
  gate_path <- file.path(diagnostics_dir, paste0(matching_stem, "_gate_g4.csv"))
  run_metadata_path <- file.path(
    diagnostics_dir, paste0(matching_stem, "_run_metadata.rds")
  )
  if (!file.exists(gate_path) || !file.exists(run_metadata_path)) {
    stop("Gate G4 or matching provenance is missing for ", matching_stem,
         ". Rerun its corresponding 2.x matching builder.", call. = FALSE)
  }
  gate <- read.csv(gate_path, stringsAsFactors = FALSE)
  if (!"gate_pass" %in% names(gate) ||
      !all(as.logical(gate$gate_pass) %in% TRUE)) {
    stop("Gate G4 did not pass for ", matching_stem, call. = FALSE)
  }
  run_metadata <- readRDS(run_metadata_path)
  expected_md5 <- unname(as.character(run_metadata$output_md5))
  observed_md5 <- unname(tools::md5sum(matched_path))
  if (length(expected_md5) != 1L || !nzchar(expected_md5) ||
      !identical(expected_md5, observed_md5)) {
    stop("Matched RDS MD5 does not match Gate G4 provenance for ",
         matching_stem, ". Rerun its matching builder.", call. = FALSE)
  }
  list(
    matching_stem = matching_stem,
    gate_path = gate_path,
    observed_md5 = observed_md5
  )
}

as_missing_factor <- function(x) {
  value <- trimws(as.character(x))
  ifelse(is.na(value) | value == "", "Missing", value)
}

panel <- readRDS(panel_path)
panel_required <- c(
  "id", "cycle", "race", "marital_status", "living_arr", "pckgr",
  "smoke", "pckyr", "bmi", "act"
)
missing_panel <- setdiff(panel_required, names(panel))
if (length(missing_panel)) {
  stop("Canonical panel is missing required time-varying columns: ",
       paste(missing_panel, collapse = ", "), call. = FALSE)
}

panel_key_check <- panel %>%
  transmute(id = as.character(id), cycle = as.character(cycle)) %>%
  count(id, cycle, name = "n_panel_rows") %>%
  filter(n_panel_rows > 1L)
if (nrow(panel_key_check)) {
  stop("Canonical panel has duplicate participant-cycle rows; refusing a many-to-many covariate join.",
       call. = FALSE)
}

tv_panel <- panel %>%
  transmute(
    id = as.character(id),
    cycle = as.character(cycle),
    panel_row_found = TRUE,
    tv_race = as_missing_factor(race),
    tv_marital = as_missing_factor(marital_status),
    tv_living = as_missing_factor(living_arr),
    tv_pckgr = as_missing_factor(pckgr),
    tv_smoke = as_missing_factor(smoke),
    tv_pckyr = as.numeric(as.character(pckyr)),
    tv_bmi = as.numeric(as.character(bmi)),
    tv_act = as.numeric(as.character(act))
  )

tv_factor_covars <- c(
  "tv_race", "tv_marital", "tv_living", "tv_pckgr", "tv_smoke"
)
tv_numeric_covars <- c("tv_pckyr", "tv_bmi", "tv_act")
tv_covars <- c(tv_factor_covars, tv_numeric_covars)

modal_level <- function(x, factor_levels) {
  counts <- tabulate(match(as.character(x), factor_levels),
                      nbins = length(factor_levels))
  if (!any(counts > 0L)) return(factor_levels[[1]])
  factor_levels[[which.max(counts)]]
}

fit_one_cohort <- function(matched_path, cohort_name, matched_provenance) {
  matched_long <- readRDS(matched_path)
  matched_required <- c(
    "Cohort", "Group", "Age_Centered", "index_age_z", "id", "cycle",
    "fi_score_nocancer", "post_own_cancer"
  )
  missing_matched <- setdiff(matched_required, names(matched_long))
  if (length(missing_matched)) {
    stop("Matched dataset is missing required columns: ",
         paste(missing_matched, collapse = ", "), call. = FALSE)
  }
  if (!cohort_name %in% unique(as.character(matched_long$Cohort))) {
    stop("Cohort '", cohort_name, "' is not present in ", basename(matched_path),
         call. = FALSE)
  }

  matched_long <- matched_long %>%
    mutate(id_join = as.character(id), cycle_join = as.character(cycle)) %>%
    left_join(tv_panel, by = c("id_join" = "id", "cycle_join" = "cycle"))
  if (any(!matched_long$panel_row_found %in% TRUE)) {
    stop("Some matched participant-cycle rows could not be linked to the canonical panel in ",
         basename(matched_path), call. = FALSE)
  }

  analysis_data <- matched_long %>%
    filter(as.character(Cohort) == cohort_name) %>%
    mutate(
      Cohort = droplevels(factor(Cohort)),
      Group = factor(Group, levels = c("Control", "Cancer Case")),
      id = factor(id),
      tv_race = factor(tv_race),
      tv_marital = factor(tv_marital),
      tv_living = factor(tv_living),
      tv_pckgr = factor(tv_pckgr),
      tv_smoke = factor(tv_smoke),
      index_age_z = as.numeric(index_age_z),
      Age_Centered = as.numeric(Age_Centered),
      fi_score_nocancer = as.numeric(fi_score_nocancer)
    ) %>%
    filter(
      !post_own_cancer,
      !is.na(fi_score_nocancer),
      !is.na(Age_Centered),
      !is.na(index_age_z),
      if_all(all_of(tv_numeric_covars), ~ !is.na(.x))
    ) %>%
    droplevels()

  if (nrow(analysis_data) == 0L || length(unique(analysis_data$Group)) < 2L) {
    stop("The time-varying M2 model requires nonempty rows in both arms for ",
         cohort_name, call. = FALSE)
  }

  prediction_min <- max(prediction_window[[1]], min(analysis_data$Age_Centered))
  prediction_max <- min(prediction_window[[2]], max(analysis_data$Age_Centered))
  if (!is.finite(prediction_min) || !is.finite(prediction_max) ||
      prediction_min >= prediction_max) {
    stop("The requested prediction window has no observed support for ",
         cohort_name, call. = FALSE)
  }

  spline_basis <- splines::ns(analysis_data$Age_Centered, df = spline_df)
  spline_terms <- paste0("S", seq_len(ncol(spline_basis)))
  analysis_data_sp <- bind_cols(
    analysis_data,
    setNames(as.data.frame(unclass(spline_basis)), spline_terms)
  )
  fixed_rhs <- as.formula(paste0(
    "~ Group * (", paste(spline_terms, collapse = " + "), ") + ",
    "index_age_z + ", paste(tv_covars, collapse = " + ")
  ))
  model_formula <- update(
    fixed_rhs,
    fi_score_nocancer ~ . + (1 + Age_Centered | id)
  )

  # This is the only model fit for this cohort. No CR covariance is attempted.
  fit <- lme4::lmer(
    model_formula,
    data = analysis_data_sp,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      calc.derivs = TRUE,
      optCtrl = list(maxfun = 2e5)
    )
  )

  beta <- lme4::fixef(fit)
  V <- as.matrix(stats::vcov(fit))
  dimnames(V) <- list(names(beta), names(beta))
  se <- sqrt(diag(V))
  coefficient_output <- data.frame(
    Cohort = cohort_name,
    model_id = "M2_time_varying",
    model_label = "M2 time-varying covariates (natural spline df = 4)",
    spline_df = spline_df,
    Term = names(beta),
    Estimate = as.numeric(beta),
    SE = as.numeric(se),
    CI_low = as.numeric(beta) - 1.96 * as.numeric(se),
    CI_high = as.numeric(beta) + 1.96 * as.numeric(se),
    vcov_type = "model-based vcov(fit)",
    row.names = NULL,
    check.names = FALSE
  )

  make_prediction_grid <- function(age_values, group_values) {
    grid <- expand.grid(
      Age_Centered = age_values,
      Group = factor(group_values, levels = c("Control", "Cancer Case")),
      stringsAsFactors = FALSE
    ) %>% mutate(index_age_z = 0)
    for (covar in tv_factor_covars) {
      grid[[covar]] <- factor(
        modal_level(analysis_data_sp[[covar]], levels(analysis_data_sp[[covar]])),
        levels = levels(analysis_data_sp[[covar]])
      )
    }
    for (covar in tv_numeric_covars) {
      grid[[covar]] <- mean(analysis_data_sp[[covar]], na.rm = TRUE)
    }
    grid
  }

  prediction_grid <- make_prediction_grid(
    seq(prediction_min, prediction_max, by = 0.1),
    c("Control", "Cancer Case")
  )
  prediction_basis <- predict(spline_basis, prediction_grid$Age_Centered)
  prediction_grid <- bind_cols(
    prediction_grid,
    setNames(as.data.frame(unclass(prediction_basis)), spline_terms)
  )
  X <- model.matrix(fixed_rhs, data = prediction_grid)[, names(beta), drop = FALSE]
  prediction_grid$pred <- as.vector(X %*% beta)
  prediction_grid$se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  prediction_grid$lwr <- prediction_grid$pred - 1.96 * prediction_grid$se
  prediction_grid$upr <- prediction_grid$pred + 1.96 * prediction_grid$se
  prediction_output <- prediction_grid %>%
    transmute(
      Cohort = cohort_name,
      model_id = "M2_time_varying",
      model_label = "M2 time-varying covariates (natural spline df = 4)",
      Group,
      Age_Centered,
      pred,
      se,
      lwr,
      upr,
      spline_df = spline_df,
      vcov_type = "model-based vcov(fit)"
    )

  list(
    model = fit,
    coefficients = coefficient_output,
    predictions = prediction_output,
    run_metadata = list(
      cohort = cohort_name,
      matched_path = normalizePath(matched_path, mustWork = TRUE),
      matched_md5 = matched_provenance$observed_md5,
      matching_gate_path = normalizePath(matched_provenance$gate_path, mustWork = TRUE),
      spline_df = spline_df,
      fitting_rows = nrow(analysis_data_sp),
      prediction_window = c(prediction_min, prediction_max),
      fixed_covariate = "index_age_z",
      time_varying_covariates = tv_covars,
      omitted_available_panel_fields = c("calor", "sat", "diet_chol", "alco"),
      random_effects = "(1 + Age_Centered | id)",
      vcov_type = "model-based vcov(fit)",
      analysis_scope = "Only M2 GLME; no event study, other models, sensitivity, CR covariance, or secondary contrasts"
    )
  )
}

model_results <- list()
for (spec_name in names(cohort_specs)) {
  matched_path <- cohort_specs[[spec_name]]
  provenance <- validate_matching_provenance(matched_path)
  matched_cohorts <- unique(as.character(readRDS(matched_path)$Cohort))
  for (cohort_name in matched_cohorts) {
    key <- paste(spec_name, cohort_name, sep = "__")
    model_results[[key]] <- fit_one_cohort(
      matched_path = matched_path,
      cohort_name = cohort_name,
      matched_provenance = provenance
    )
  }
}

coefficient_output <- bind_rows(lapply(model_results, `[[`, "coefficients"))
prediction_output <- bind_rows(lapply(model_results, `[[`, "predictions"))
run_metadata_output <- list(
  schema_version = 1L,
  analysis = "4.2.2 other active cancer cohorts M2 GLME with time-varying covariates",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  cohort_files = unname(cohort_specs),
  cohorts = lapply(model_results, `[[`, "run_metadata"),
  spline_df = spline_df,
  time_varying_covariates = tv_covars,
  omitted_available_panel_fields = c("calor", "sat", "diet_chol", "alco"),
  vcov_type = "model-based vcov(fit)"
)

write.csv(
  coefficient_output,
  file.path(results_dir, paste0(out_prefix, "_fixed_effects_CI.csv")),
  row.names = FALSE
)
write.csv(
  prediction_output,
  file.path(results_dir, paste0(out_prefix, "_predicted_trajectories.csv")),
  row.names = FALSE
)
saveRDS(
  list(
    models = lapply(model_results, `[[`, "model"),
    coefficients = coefficient_output,
    predictions = prediction_output,
    run_metadata = run_metadata_output
  ),
  file.path(results_dir, paste0(out_prefix, "_models.rds"))
)

message("Saved only the other-cohort M2 time-varying df-4 GLME outputs with prefix '",
        out_prefix, "' to ", results_dir)
