# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.4_ipcw_glme_spline_subtype_cohorts.R
# Author:  Nemo Zhou
# Date started:      2026-06-30
# Date last updated: 2026-07-17 (subtype cohorts now related cancer vs controls)
#
# Purpose:
#   Fits IPCW-weighted natural-spline GLME trajectory models for the active
#   cancer subtype cohort families, excluding prostate:
#     * high-burden and low/moderate-burden cancer cohorts,
#     * smoking-related cancer versus cancer-free controls, and
#     * obesity-related cancer versus cancer-free controls.
#
#   Each cohort is fit separately as cancer cases versus its risk-set matched
#   controls, using the IPCW-augmented matched datasets from 7.3. Survival-
#   stratified cohorts are not included because they classify cases by post-index
#   survival/prognosis rather than cancer subtype.
#
#   Model within each cohort:
#     fi_score_nocancer ~ Group * ns(Age_Centered, df = 4)
#                       + index_age_z + base_race + base_marital + base_pckgr
#                       + (1 + Age_Centered | id),
#     fit with weights = sw_ipcw.
#
#   The spline basis is materialized before model fitting so the same basis can
#   be used for prediction and cancer-minus-control contrasts. Confidence
#   intervals use model-based SEs because `clubSandwich` does not currently
#   support CR2/CR0 covariance estimation for prior-weighted `lmer` models.
#
# Inputs:
#   Data/riskset_matched_analysis_long_ipcw.rds
#   Data/riskset_matched_smoking_long_ipcw.rds
#   Data/riskset_matched_obesity_long_ipcw.rds
# Outputs:
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_fixed_effects_CI.csv
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_predicted_trajectories.csv
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_group_difference.csv
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_models.rds
# =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(splines)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

out_prefix <- "7.4_ipcw_glme_spline_subtypes"
window_yrs <- 16
spline_df  <- 4

needed_cols <- c(
  "Cohort", "Group", "Age_Centered", "index_age_z",
  "base_race", "base_marital", "base_pckgr", "id",
  "fi_score_nocancer", "sw_ipcw"
)

subtype_specs <- list(
  list(
    analysis = "Cancer Burden",
    path = file.path(data_dir, "riskset_matched_analysis_long_ipcw.rds"),
    builder_script = "Code/2_data_analysis/7.3_ipcw_subtype_cohorts.R"
  ),
  list(
    analysis = "Smoking-Related Cancer",
    path = file.path(data_dir, "riskset_matched_smoking_long_ipcw.rds"),
    builder_script = "Code/2_data_analysis/7.3_ipcw_subtype_cohorts.R"
  ),
  list(
    analysis = "Obesity-Related Cancer",
    path = file.path(data_dir, "riskset_matched_obesity_long_ipcw.rds"),
    builder_script = "Code/2_data_analysis/7.3_ipcw_subtype_cohorts.R"
  )
)

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

get_vcov <- function(model) {
  fe <- names(lme4::fixef(model))
  message("  Prior weights detected; using model-based covariance because ",
          "clubSandwich does not support prior-weighted lmer models.")
  V <- as.matrix(vcov(model))
  dimnames(V) <- list(fe, fe)
  list(V = V, type = "model-based")
}

make_ref_grid <- function(d, age_values, group_values) {
  expand.grid(
    Age_Centered = age_values,
    Group = factor(group_values, levels = c("Control", "Cancer Case")),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      index_age_z = 0,
      base_race = factor(levels(d$base_race)[1], levels = levels(d$base_race)),
      base_marital = factor(levels(d$base_marital)[1],
                            levels = levels(d$base_marital)),
      base_pckgr = factor(levels(d$base_pckgr)[1], levels = levels(d$base_pckgr))
    )
}

fit_weighted_spline <- function(d, analysis_name, cohort_name) {
  d <- d %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Cancer Case")),
      id = factor(id),
      base_race = factor(base_race),
      base_marital = factor(base_marital),
      base_pckgr = factor(base_pckgr),
      sw_ipcw = as.numeric(sw_ipcw)
    ) %>%
    filter(
      abs(Age_Centered) <= window_yrs,
      !is.na(fi_score_nocancer),
      !is.na(sw_ipcw),
      sw_ipcw > 0
    ) %>%
    droplevels()

  if (length(unique(d$Group)) < 2) {
    stop("Both Control and Cancer Case rows are required for ", analysis_name,
         " / ", cohort_name, ".")
  }

  cat("\n======== IPCW-weighted spline GLME:", analysis_name, "/", cohort_name, "========\n")
  cat("Analytic rows:", nrow(d), "\n")
  print(table(Group = d$Group))
  cat("IPCW summary in analysis window:\n")
  print(summary(d$sw_ipcw))

  boundary_knots <- range(d$Age_Centered, na.rm = TRUE)
  inner_knots <- as.numeric(quantile(
    d$Age_Centered,
    probs = seq(0, 1, length.out = spline_df + 1)[-c(1, spline_df + 1)],
    na.rm = TRUE
  ))
  spline_basis <- ns(
    d$Age_Centered,
    knots = inner_knots,
    Boundary.knots = boundary_knots
  )
  spline_terms <- paste0("S", seq_len(ncol(spline_basis)))
  d_sp <- bind_cols(d, setNames(as.data.frame(unclass(spline_basis)), spline_terms))

  fixed_rhs <- as.formula(paste0(
    "~ Group * (", paste(spline_terms, collapse = " + "), ") + ",
    "index_age_z + base_race + base_marital + base_pckgr"
  ))
  model_formula <- update(
    fixed_rhs,
    fi_score_nocancer ~ . + (1 + Age_Centered | id)
  )

  m <- lme4::lmer(
    model_formula,
    data = d_sp,
    weights = sw_ipcw,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  )
  print(summary(m), correlation = FALSE)

  beta <- fixef(m)
  Vobj <- get_vcov(m)
  V <- Vobj$V
  coef_se <- sqrt(diag(V))

  coef_tab <- data.frame(
    Analysis = analysis_name,
    Cohort = cohort_name,
    Term = names(beta),
    Estimate = as.numeric(beta),
    SE = coef_se,
    CI_low = as.numeric(beta) - 1.96 * coef_se,
    CI_high = as.numeric(beta) + 1.96 * coef_se,
    vcov_type = Vobj$type,
    row.names = NULL,
    check.names = FALSE
  )

  pred_grid <- make_ref_grid(
    d_sp,
    age_values = seq(-window_yrs, window_yrs, by = 0.1),
    group_values = c("Control", "Cancer Case")
  )
  pred_basis <- predict(spline_basis, pred_grid$Age_Centered)
  pred_grid <- bind_cols(pred_grid, setNames(as.data.frame(pred_basis), spline_terms))
  X <- model.matrix(fixed_rhs, data = pred_grid)[, names(beta), drop = FALSE]
  pred_grid$pred <- as.vector(X %*% beta)
  pred_grid$se <- sqrt(rowSums((X %*% V) * X))
  pred_grid$lwr <- pred_grid$pred - 1.96 * pred_grid$se
  pred_grid$upr <- pred_grid$pred + 1.96 * pred_grid$se
  pred_grid$Analysis <- analysis_name
  pred_grid$Cohort <- cohort_name
  pred_grid$plot_label <- paste(analysis_name, cohort_name, sep = ": ")
  pred_grid$vcov_type <- Vobj$type
  pred_out <- pred_grid %>%
    select(Analysis, Cohort, plot_label, Group, Age_Centered,
           pred, se, lwr, upr, vcov_type)

  contrast_grid <- make_ref_grid(
    d_sp,
    age_values = seq(-window_yrs, window_yrs, by = 0.1),
    group_values = "Cancer Case"
  )
  contrast_basis <- predict(spline_basis, contrast_grid$Age_Centered)
  contrast_grid <- bind_cols(
    contrast_grid,
    setNames(as.data.frame(contrast_basis), spline_terms)
  )
  contrast_case <- contrast_grid
  contrast_ctrl <- contrast_grid %>%
    mutate(Group = factor("Control", levels = c("Control", "Cancer Case")))
  X_case <- model.matrix(fixed_rhs, data = contrast_case)[, names(beta), drop = FALSE]
  X_ctrl <- model.matrix(fixed_rhs, data = contrast_ctrl)[, names(beta), drop = FALSE]
  X_diff <- X_case - X_ctrl
  diff_est <- as.vector(X_diff %*% beta)
  diff_se <- sqrt(rowSums((X_diff %*% V) * X_diff))
  diff_out <- data.frame(
    Analysis = analysis_name,
    Cohort = cohort_name,
    plot_label = paste(analysis_name, cohort_name, sep = ": "),
    Age_Centered = contrast_grid$Age_Centered,
    contrast = "Cancer Case - Control",
    estimate = diff_est,
    se = diff_se,
    CI_low = diff_est - 1.96 * diff_se,
    CI_high = diff_est + 1.96 * diff_se,
    vcov_type = Vobj$type
  )

  list(
    model = m,
    coefficients = coef_tab,
    predictions = pred_out,
    group_difference = diff_out,
    spline_knots = inner_knots,
    boundary_knots = boundary_knots,
    vcov_type = Vobj$type
  )
}

model_results <- list()

for (spec in subtype_specs) {
  if (!file.exists(spec$path)) {
    stop("IPCW-augmented subtype dataset not found at ", spec$path,
         ". Run ", spec$builder_script, " first.")
  }
  dat <- readRDS(spec$path)
  missing_cols <- setdiff(needed_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Dataset ", basename(spec$path), " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  dat <- dat %>%
    mutate(Cohort = droplevels(factor(Cohort)))

  for (cohort_name in levels(dat$Cohort)) {
    key <- paste(spec$analysis, cohort_name, sep = "__")
    model_results[[key]] <- fit_weighted_spline(
      d = dat %>% filter(Cohort == cohort_name),
      analysis_name = spec$analysis,
      cohort_name = cohort_name
    )
  }
}

coef_all <- bind_rows(lapply(model_results, `[[`, "coefficients"))
pred_all <- bind_rows(lapply(model_results, `[[`, "predictions"))
diff_all <- bind_rows(lapply(model_results, `[[`, "group_difference"))

write.csv(
  coef_all,
  file.path(results_dir, paste0(out_prefix, "_fixed_effects_CI.csv")),
  row.names = FALSE
)
write.csv(
  pred_all,
  file.path(results_dir, paste0(out_prefix, "_predicted_trajectories.csv")),
  row.names = FALSE
)
write.csv(
  diff_all,
  file.path(results_dir, paste0(out_prefix, "_group_difference.csv")),
  row.names = FALSE
)

group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")

p_traj <- ggplot(pred_all, aes(x = Age_Centered, y = pred, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.05) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  facet_wrap(~ plot_label, scales = "free_y") +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW-Weighted Frailty Trajectories by Cancer Subtype Cohort",
    subtitle = paste0("Natural-spline GLME (df = ", spline_df,
                      "); weights = sw_ipcw; model-based 95% CIs"),
    x = "Years relative to index",
    y = "Predicted frailty index (no-cancer FI)",
    color = NULL,
    fill = NULL
  ) +
  scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 8))

p_diff <- ggplot(diff_all, aes(x = Age_Centered, y = estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = "#de2d26", alpha = 0.18) +
  geom_line(color = "#de2d26", linewidth = 1.05) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  facet_wrap(~ plot_label, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW-Weighted Cancer-Control Difference by Cancer Subtype Cohort",
    subtitle = "Cancer Case - Control from natural-spline GLME; model-based 95% CIs",
    x = "Years relative to index",
    y = "Difference in predicted frailty index"
  ) +
  scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 8))

saveRDS(
  model_results,
  file.path(results_dir, paste0(out_prefix, "_models.rds"))
)

cat("\nSaved subtype IPCW-weighted spline GLME outputs with prefix '",
    out_prefix, "' to: ", results_dir, "\n", sep = "")
