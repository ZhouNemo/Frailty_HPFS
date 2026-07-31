# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.2_ipcw_glme_spline_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-30
# Date last updated: 2026-07-28 (prediction factors use participant-modal levels)
#
# Purpose:
#   Fits the IPCW-weighted natural-spline GLME for the OVERALL incident cancer
#   cohort: all cancer cases versus risk-set matched cancer-free controls. This
#   is the weighted downstream trajectory model after:
#     * 2.5_create_riskset_overall.R builds the overall matched cohort, and
#     * 7.1_ipcw_overall.R attaches stabilized post-index IPCW (`sw_ipcw`).
#
#   Model:
#     fi_score_nocancer ~ Group * ns(Age_Centered, df = 4)
#                       + index_age_z + base_race + base_marital + base_pckgr
#                       + (1 + Age_Centered | id),
#     fit with weights = sw_ipcw.
#
#   The spline basis is materialized before model fitting so the same basis can
#   be used for prediction and for cancer-minus-control contrasts. Confidence
#   intervals use model-based SEs because `clubSandwich` does not currently
#   support CR2/CR0 covariance estimation for prior-weighted `lmer` models.
#
# Inputs:
#   Data/riskset_matched_overall_long_ipcw.rds      (from 7.1)
# Outputs:
#   Results/cancer/data/7.2_ipcw_glme_spline_fixed_effects_CI.csv
#   Results/cancer/data/7.2_ipcw_glme_spline_predicted_trajectories.csv
#   Results/cancer/data/7.2_ipcw_glme_spline_group_difference.csv
#   Results/cancer/data/7.2_ipcw_glme_spline_model.rds
# =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(splines)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

matched_path <- file.path(data_dir, "riskset_matched_overall_long_ipcw.rds")
out_prefix   <- "7.2_ipcw_glme_spline"
window_yrs   <- 16
spline_df    <- 4

needed_cols <- c(
  "Cohort", "Group", "Age_Centered", "index_age_z",
  "base_race", "base_marital", "base_pckgr", "id",
  "fi_score_nocancer", "sw_ipcw"
)

if (!file.exists(matched_path)) {
  stop(
    "IPCW-augmented overall matched dataset not found at ", matched_path,
    ". Run Code/2_data_analysis/7.1_ipcw_overall.R first."
  )
}
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

matched_long <- readRDS(matched_path)
missing_cols <- setdiff(needed_cols, names(matched_long))
if (length(missing_cols) > 0) {
  stop("Matched IPCW dataset is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

analysis_data <- matched_long %>%
  mutate(
    Cohort = droplevels(factor(Cohort)),
    Group = factor(Group, levels = c("Control", "Cancer Case")),
    id = factor(id),
    base_race = factor(base_race),
    base_marital = factor(base_marital),
    base_pckgr = factor(base_pckgr),
    sw_ipcw = as.numeric(sw_ipcw)
  ) %>%
  filter(
    Cohort == "All Cancer Cohort",
    abs(Age_Centered) <= window_yrs,
    !is.na(fi_score_nocancer),
    !is.na(sw_ipcw),
    sw_ipcw > 0
  ) %>%
  droplevels()

participant_modal_level <- function(d, column) {
  factor_levels <- levels(d[[column]])
  participant_values <- d[!duplicated(d[c("id", column)]), c("id", column)][[column]]
  participant_values <- as.character(participant_values[!is.na(participant_values)])
  counts <- tabulate(match(participant_values, factor_levels),
                      nbins = length(factor_levels))
  if (!any(counts > 0L)) return(factor_levels[[1]])
  factor_levels[[which.max(counts)]]
}

cat("\nIPCW-weighted overall spline GLME analytic rows:", nrow(analysis_data), "\n")
print(table(Group = analysis_data$Group))
cat("\nIPCW summary in analysis window:\n")
print(summary(analysis_data$sw_ipcw))

if (length(unique(analysis_data$Group)) < 2) {
  stop("Both Control and Cancer Case rows are required for the overall model.")
}

# Fixed spline basis, reused for predictions and contrasts.
boundary_knots <- range(analysis_data$Age_Centered, na.rm = TRUE)
inner_knots <- as.numeric(quantile(
  analysis_data$Age_Centered,
  probs = seq(0, 1, length.out = spline_df + 1)[-c(1, spline_df + 1)],
  na.rm = TRUE
))
spline_basis <- ns(
  analysis_data$Age_Centered,
  knots = inner_knots,
  Boundary.knots = boundary_knots
)
spline_terms <- paste0("S", seq_len(ncol(spline_basis)))

analysis_data <- bind_cols(
  analysis_data,
  setNames(as.data.frame(unclass(spline_basis)), spline_terms)
)

fixed_rhs <- as.formula(paste0(
  "~ Group * (", paste(spline_terms, collapse = " + "), ") + ",
  "index_age_z + base_race + base_marital + base_pckgr"
))

model_formula <- update(
  fixed_rhs,
  fi_score_nocancer ~ . + (1 + Age_Centered | id)
)

cat("\nFitting IPCW-weighted natural-spline GLME ...\n")
m_ipcw <- lme4::lmer(
  model_formula,
  data = analysis_data,
  weights = sw_ipcw,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
)
print(summary(m_ipcw), correlation = FALSE)

get_vcov <- function(model, cluster) {
  fe <- names(lme4::fixef(model))
  has_prior_weights <- any(abs(stats::weights(model) - 1) > .Machine$double.eps^0.5)
  if (has_prior_weights) {
    message("  Prior weights detected; using model-based covariance because ",
            "clubSandwich does not support prior-weighted lmer models.")
  } else if (requireNamespace("clubSandwich", quietly = TRUE)) {
    for (ty in c("CR2", "CR0")) {
      V <- tryCatch(
        as.matrix(clubSandwich::vcovCR(model, cluster = cluster, type = ty)),
        error = function(e) {
          message("  clubSandwich ", ty, " failed: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(V)) {
        dimnames(V) <- list(fe, fe)
        return(list(V = V, type = paste0(ty, " (cluster-robust on id)")))
      }
    }
  } else {
    message("  clubSandwich not installed; using model-based covariance.")
  }
  V <- as.matrix(vcov(model))
  dimnames(V) <- list(fe, fe)
  list(V = V, type = "model-based")
}

beta <- fixef(m_ipcw)
Vobj <- get_vcov(m_ipcw, analysis_data$id)
V <- Vobj$V
se <- sqrt(diag(V))

coef_tab <- data.frame(
  Term = names(beta),
  Estimate = as.numeric(beta),
  SE = se,
  CI_low = as.numeric(beta) - 1.96 * se,
  CI_high = as.numeric(beta) + 1.96 * se,
  vcov_type = Vobj$type,
  row.names = NULL,
  check.names = FALSE
)
write.csv(
  coef_tab,
  file.path(results_dir, paste0(out_prefix, "_fixed_effects_CI.csv")),
  row.names = FALSE
)

make_ref_grid <- function(age_values, group_values) {
  race_ref <- participant_modal_level(analysis_data, "base_race")
  marital_ref <- participant_modal_level(analysis_data, "base_marital")
  pack_year_ref <- participant_modal_level(analysis_data, "base_pckgr")

  expand.grid(
    Age_Centered = age_values,
    Group = factor(group_values, levels = c("Control", "Cancer Case")),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      index_age_z = 0,
      base_race = factor(race_ref,
                         levels = levels(analysis_data$base_race)),
      base_marital = factor(marital_ref,
                            levels = levels(analysis_data$base_marital)),
      base_pckgr = factor(pack_year_ref,
                          levels = levels(analysis_data$base_pckgr))
    )
}

prediction_grid <- make_ref_grid(
  age_values = seq(-window_yrs, window_yrs, by = 0.1),
  group_values = c("Control", "Cancer Case")
)
pred_basis <- predict(spline_basis, prediction_grid$Age_Centered)
prediction_grid <- bind_cols(
  prediction_grid,
  setNames(as.data.frame(pred_basis), spline_terms)
)

X <- model.matrix(fixed_rhs, data = prediction_grid)[, names(beta), drop = FALSE]
prediction_grid$pred <- as.vector(X %*% beta)
prediction_grid$se <- sqrt(rowSums((X %*% V) * X))
prediction_grid$lwr <- prediction_grid$pred - 1.96 * prediction_grid$se
prediction_grid$upr <- prediction_grid$pred + 1.96 * prediction_grid$se
prediction_grid$Cohort <- "All Cancer Cohort"
prediction_grid$vcov_type <- Vobj$type

pred_out <- prediction_grid %>%
  select(Cohort, Group, Age_Centered, pred, se, lwr, upr, vcov_type)
write.csv(
  pred_out,
  file.path(results_dir, paste0(out_prefix, "_predicted_trajectories.csv")),
  row.names = FALSE
)

contrast_grid <- make_ref_grid(
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
  Cohort = "All Cancer Cohort",
  Age_Centered = contrast_grid$Age_Centered,
  contrast = "Cancer Case - Control",
  estimate = diff_est,
  se = diff_se,
  CI_low = diff_est - 1.96 * diff_se,
  CI_high = diff_est + 1.96 * diff_se,
  vcov_type = Vobj$type
)
write.csv(
  diff_out,
  file.path(results_dir, paste0(out_prefix, "_group_difference.csv")),
  row.names = FALSE
)

group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")

p_traj <- ggplot(pred_out, aes(x = Age_Centered, y = pred, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW-Weighted Frailty Trajectories: Overall Incident Cancer vs Control",
    subtitle = paste0(
      "Natural-spline GLME (df = ", spline_df,
      "); weights = sw_ipcw; ", Vobj$type, " 95% CIs"
    ),
    x = "Years relative to index",
    y = "Predicted frailty index (no-cancer FI)",
    color = NULL,
    fill = NULL
  ) +
  scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 4))

p_diff <- ggplot(diff_out, aes(x = Age_Centered, y = estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = "#de2d26", alpha = 0.18) +
  geom_line(color = "#de2d26", linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW-Weighted Cancer-Control Difference in Frailty",
    subtitle = paste0("Cancer Case - Control from natural-spline GLME; ",
                      Vobj$type, " 95% CIs"),
    x = "Years relative to index",
    y = "Difference in predicted frailty index"
  ) +
  scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 4))

saveRDS(
  list(
    model = m_ipcw,
    coefficients = coef_tab,
    predictions = pred_out,
    group_difference = diff_out,
    spline_knots = inner_knots,
    boundary_knots = boundary_knots,
    vcov_type = Vobj$type
  ),
  file.path(results_dir, paste0(out_prefix, "_model.rds"))
)

cat("\nSaved IPCW-weighted spline GLME outputs with prefix '", out_prefix,
    "' to: ", results_dir, "\n", sep = "")
