# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.5_compare_ipcw_vs_unweighted_spline.R
# Author:  Nemo Zhou
# Date started:      2026-06-30
# Date last updated: 2026-07-17 (visual-data outputs moved; PNG writes removed)
#
# Purpose:
#   Compares the key natural-spline GLME estimates from the same risk-set matched
#   cohorts with and without IPCW:
#     * unweighted spline models from 4.1-4.5, and
#     * IPCW-weighted spline models from 7.2 and 7.4.
#
#   The comparison includes:
#     1. Cancer-case predicted FI trajectories by cohort.
#     2. Cancer-minus-control spline contrasts by cohort.
#     3. Compact key-time summaries at -16, -12, -8, -4, 0, +4, +8, +12,
#        and +16 years relative to index.
#
#   Note: the unweighted 4.x scripts save spline predicted trajectories but do
#   not save the spline model covariance or a spline contrast table. Therefore,
#   unweighted cancer-minus-control contrast point estimates are computed from
#   saved predicted means, and their contrast SE/CI fields are left missing.
#
# Inputs:
#   Results/cancer/data/4.1_spline_predicted_trajectories.csv
#   Results/cancer/data/4.2_spline_predicted_trajectories.csv
#   Results/cancer/data/4.3_spline_predicted_trajectories.csv
#   Results/cancer/data/4.5_spline_predicted_trajectories.csv
#   Results/cancer/data/7.2_ipcw_glme_spline_predicted_trajectories.csv
#   Results/cancer/data/7.2_ipcw_glme_spline_group_difference.csv
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_predicted_trajectories.csv
#   Results/cancer/data/7.4_ipcw_glme_spline_subtypes_group_difference.csv
# Outputs:
#   Results/cancer/data/7.5_ipcw_vs_unweighted_case_all_grid.csv
#   Results/cancer/data/7.5_ipcw_vs_unweighted_contrast_all_grid.csv
#   Results/cancer/data/7.5_ipcw_vs_unweighted_case_key_times.csv
#   Results/cancer/data/7.5_ipcw_vs_unweighted_contrast_key_times.csv
#   Results/cancer/data/7.5_ipcw_vs_unweighted_case_delta_key_times.csv
#   Results/cancer/data/7.5_ipcw_vs_unweighted_contrast_delta_key_times.csv
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
results_dir <- file.path(project_dir, "Results", "cancer", "data")

key_times <- c(-16, -12, -8, -4, 0, 4, 8, 12, 16)

read_required_csv <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

read_unweighted_pred <- function(path, analysis_name) {
  read_required_csv(path) %>%
    mutate(
      Analysis = analysis_name,
      plot_label = paste(Analysis, Cohort, sep = ": "),
      Model = "Unweighted",
      se = NA_real_,
      vcov_type = "4.x saved predicted trajectory CI"
    ) %>%
    select(Analysis, Cohort, plot_label, Model, Group, Age_Centered,
           pred, se, lwr, upr, vcov_type)
}

unweighted_pred <- bind_rows(
  read_unweighted_pred(
    file.path(results_dir, "4.5_spline_predicted_trajectories.csv"),
    "Overall"
  ),
  read_unweighted_pred(
    file.path(results_dir, "4.1_spline_predicted_trajectories.csv"),
    "Cancer Burden"
  ),
  read_unweighted_pred(
    file.path(results_dir, "4.2_spline_predicted_trajectories.csv"),
    "Smoking-Related Cancer"
  ),
  read_unweighted_pred(
    file.path(results_dir, "4.3_spline_predicted_trajectories.csv"),
    "Obesity-Related Cancer"
  )
)

weighted_overall_pred <- read_required_csv(
  file.path(results_dir, "7.2_ipcw_glme_spline_predicted_trajectories.csv")
) %>%
  mutate(
    Analysis = "Overall",
    plot_label = paste(Analysis, Cohort, sep = ": "),
    Model = "IPCW"
  ) %>%
  select(Analysis, Cohort, plot_label, Model, Group, Age_Centered,
         pred, se, lwr, upr, vcov_type)

weighted_subtype_pred <- read_required_csv(
  file.path(results_dir, "7.4_ipcw_glme_spline_subtypes_predicted_trajectories.csv")
) %>%
  mutate(Model = "IPCW") %>%
  select(Analysis, Cohort, plot_label, Model, Group, Age_Centered,
         pred, se, lwr, upr, vcov_type)

weighted_pred <- bind_rows(weighted_overall_pred, weighted_subtype_pred)

case_all <- bind_rows(unweighted_pred, weighted_pred) %>%
  filter(Group == "Cancer Case") %>%
  mutate(Age_Centered = round(Age_Centered, 1)) %>%
  arrange(Analysis, Cohort, Age_Centered, Model)

make_unweighted_contrast <- function(pred) {
  case <- pred %>%
    filter(Group == "Cancer Case") %>%
    select(Analysis, Cohort, plot_label, Age_Centered, case_pred = pred)
  ctrl <- pred %>%
    filter(Group == "Control") %>%
    select(Analysis, Cohort, Age_Centered, ctrl_pred = pred)
  case %>%
    left_join(ctrl, by = c("Analysis", "Cohort", "Age_Centered")) %>%
    transmute(
      Analysis, Cohort, plot_label,
      Model = "Unweighted",
      Age_Centered = round(Age_Centered, 1),
      contrast = "Cancer Case - Control",
      estimate = case_pred - ctrl_pred,
      se = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      vcov_type = "contrast SE/CI unavailable from saved 4.x spline predictions"
    )
}

unweighted_contrast <- make_unweighted_contrast(unweighted_pred)

weighted_overall_contrast <- read_required_csv(
  file.path(results_dir, "7.2_ipcw_glme_spline_group_difference.csv")
) %>%
  mutate(
    Analysis = "Overall",
    plot_label = paste(Analysis, Cohort, sep = ": "),
    Model = "IPCW",
    Age_Centered = round(Age_Centered, 1)
  ) %>%
  select(Analysis, Cohort, plot_label, Model, Age_Centered, contrast,
         estimate, se, CI_low, CI_high, vcov_type)

weighted_subtype_contrast <- read_required_csv(
  file.path(results_dir, "7.4_ipcw_glme_spline_subtypes_group_difference.csv")
) %>%
  mutate(
    Model = "IPCW",
    Age_Centered = round(Age_Centered, 1)
  ) %>%
  select(Analysis, Cohort, plot_label, Model, Age_Centered, contrast,
         estimate, se, CI_low, CI_high, vcov_type)

contrast_all <- bind_rows(
  unweighted_contrast,
  weighted_overall_contrast,
  weighted_subtype_contrast
) %>%
  arrange(Analysis, Cohort, Age_Centered, Model)

case_key_long <- case_all %>%
  filter(Age_Centered %in% key_times)
contrast_key_long <- contrast_all %>%
  filter(Age_Centered %in% key_times)

make_case_delta <- function(x) {
  unweighted <- x %>%
    filter(Model == "Unweighted") %>%
    select(Analysis, Cohort, plot_label, Age_Centered,
           unweighted_case_pred = pred,
           unweighted_case_lwr = lwr,
           unweighted_case_upr = upr)
  ipcw <- x %>%
    filter(Model == "IPCW") %>%
    select(Analysis, Cohort, Age_Centered,
           ipcw_case_pred = pred,
           ipcw_case_lwr = lwr,
           ipcw_case_upr = upr)
  unweighted %>%
    inner_join(ipcw, by = c("Analysis", "Cohort", "Age_Centered")) %>%
    mutate(
      ipcw_minus_unweighted = ipcw_case_pred - unweighted_case_pred,
      relative_change_pct = 100 * ipcw_minus_unweighted / unweighted_case_pred
    ) %>%
    arrange(Analysis, Cohort, Age_Centered)
}

make_contrast_delta <- function(x) {
  unweighted <- x %>%
    filter(Model == "Unweighted") %>%
    select(Analysis, Cohort, plot_label, Age_Centered,
           unweighted_case_minus_control = estimate)
  ipcw <- x %>%
    filter(Model == "IPCW") %>%
    select(Analysis, Cohort, Age_Centered,
           ipcw_case_minus_control = estimate,
           ipcw_CI_low = CI_low,
           ipcw_CI_high = CI_high)
  unweighted %>%
    inner_join(ipcw, by = c("Analysis", "Cohort", "Age_Centered")) %>%
    mutate(
      ipcw_minus_unweighted = ipcw_case_minus_control - unweighted_case_minus_control,
      relative_change_pct = ifelse(
        abs(unweighted_case_minus_control) > 1e-8,
        100 * ipcw_minus_unweighted / unweighted_case_minus_control,
        NA_real_
      )
    ) %>%
    arrange(Analysis, Cohort, Age_Centered)
}

case_delta_key <- make_case_delta(case_key_long)
contrast_delta_key <- make_contrast_delta(contrast_key_long)

write.csv(
  case_all,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_case_all_grid.csv"),
  row.names = FALSE
)
write.csv(
  contrast_all,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_contrast_all_grid.csv"),
  row.names = FALSE
)
write.csv(
  case_key_long,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_case_key_times.csv"),
  row.names = FALSE
)
write.csv(
  contrast_key_long,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_contrast_key_times.csv"),
  row.names = FALSE
)
write.csv(
  case_delta_key,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_case_delta_key_times.csv"),
  row.names = FALSE
)
write.csv(
  contrast_delta_key,
  file.path(results_dir, "7.5_ipcw_vs_unweighted_contrast_delta_key_times.csv"),
  row.names = FALSE
)

p_case <- ggplot(case_delta_key,
                 aes(x = Age_Centered, y = ipcw_minus_unweighted)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_line(color = "#2b8cbe", linewidth = 0.9) +
  geom_point(color = "#2b8cbe", size = 1.8) +
  facet_wrap(~ plot_label, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW Minus Unweighted Cancer-Case Predicted FI",
    subtitle = "Natural-spline GLME estimates at selected relative times",
    x = "Years relative to index",
    y = "IPCW - unweighted predicted FI among cancer cases"
  ) +
  scale_x_continuous(breaks = key_times)

p_contrast <- ggplot(contrast_delta_key,
                     aes(x = Age_Centered, y = ipcw_minus_unweighted)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_line(color = "#de2d26", linewidth = 0.9) +
  geom_point(color = "#de2d26", size = 1.8) +
  facet_wrap(~ plot_label, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = "IPCW Minus Unweighted Cancer-Control Difference",
    subtitle = "Natural-spline GLME cancer-minus-control estimates at selected relative times",
    x = "Years relative to index",
    y = "IPCW - unweighted cancer-control difference"
  ) +
  scale_x_continuous(breaks = key_times)

cat("\nSaved IPCW-vs-unweighted spline comparison outputs to: ",
    results_dir, "\n", sep = "")
cat("\nCancer-case key-time comparison preview:\n")
print(head(case_delta_key, 12), row.names = FALSE)
cat("\nCancer-minus-control key-time comparison preview:\n")
print(head(contrast_delta_key, 12), row.names = FALSE)
