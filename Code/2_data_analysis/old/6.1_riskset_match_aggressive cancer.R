# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  6.1_riskset_match_aggressive cancer.R
# Author:  Nemo Zhou
# Date started:      2026-06-28
# Date last updated: 2026-06-29
#
# Purpose:
#   Fits the piecewise mixed-effects trajectory model for the high-burden versus
#   low/moderate-burden risk-set matched cohorts. Risk-set matching is no longer
#   created in this script; run Code/2_data_analysis/2.1_create_riskset_high_low_burden.R
#   first to overwrite Data/riskset_matched_analysis_long.rds.
#
# Input:
#   Data/riskset_matched_analysis_long.rds
#
# Outputs:
#   Results/cancer/6.1_riskset_trajectories_by_cohort.png
#   Results/cancer/6.1_riskset_index_age_density.png
# =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(gtsummary)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer")
matched_path <- file.path(data_dir, "riskset_matched_analysis_long.rds")

WINDOW_YRS <- 20
needed_cols <- c("Cohort", "Group", "Age_Centered", "Age_Post", "Post",
                 "index_age", "index_age_z", "base_race", "base_marital",
                 "base_living", "base_pckgr", "id", "match_set",
                 "fi_score_nocancer")

if (!file.exists(matched_path)) {
  stop("Matched dataset not found. Run Code/2_data_analysis/2.1_create_riskset_high_low_burden.R first.")
}

matched_long <- readRDS(matched_path)
missing_cols <- setdiff(needed_cols, names(matched_long))
if (length(missing_cols) > 0) {
  stop("Matched dataset is missing required columns: ", paste(missing_cols, collapse = ", "))
}

matched_long <- matched_long %>%
  mutate(
    Cohort = droplevels(factor(Cohort)),
    Group = factor(Group, levels = c("Control", "Cancer Case")),
    id = factor(id),
    match_set = factor(match_set)
  )

cat("\nMatched analytic rows:", nrow(matched_long), "\n")
print(table(Cohort = matched_long$Cohort, Group = matched_long$Group))

record_summary <- matched_long %>%
  group_by(Cohort, match_set, id, Group) %>%
  summarize(
    n_pre = sum(Post == 0L),
    n_post = sum(Post == 1L),
    index_age = first(index_age),
    fi_baseline = first(fi_score_nocancer[order(Age_Centered)]),
    base_race = first(base_race),
    base_marital = first(base_marital),
    base_living = first(base_living),
    base_pckgr = first(base_pckgr),
    .groups = "drop"
  )

record_summary %>%
  group_by(Cohort, Group) %>%
  summarize(
    n_records = n(),
    n_persons = n_distinct(id),
    mean_pre_obs = round(mean(n_pre), 2),
    mean_post_obs = round(mean(n_post), 2),
    mean_index_age = round(mean(index_age), 1),
    .groups = "drop"
  ) %>%
  print()

for (ch in levels(matched_long$Cohort)) {
  cat("\n==== Table 1:", ch, "====\n")
  t1 <- record_summary %>%
    filter(Cohort == ch) %>%
    select(Group, index_age, fi_baseline, base_race, base_marital, base_living, base_pckgr) %>%
    tbl_summary(
      by = Group,
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)"),
      label = list(index_age ~ "Index Age (years)",
                   fi_baseline ~ "Baseline Frailty Index (distal)"),
      missing_text = "(Missing)"
    ) %>%
    add_p() %>%
    add_overall() %>%
    bold_labels()
  print(t1)
}

group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")

p_index <- ggplot(record_summary, aes(x = index_age, fill = Group, color = Group)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  facet_wrap(~ Cohort) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 13) +
  labs(title = "Index-Age Distribution by Group (Separate Risk-Set Matching)",
       x = "Index age (years)", y = "Density", fill = NULL, color = NULL) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p_index)

lme_formula <- fi_score_nocancer ~
  Group * Age_Centered + Group * Age_Post +
  index_age_z + base_race + base_marital + base_pckgr +
  (1 + Age_Centered | id) + (1 | match_set)

models <- list()
for (ch in levels(matched_long$Cohort)) {
  d <- matched_long %>% filter(Cohort == ch) %>% droplevels()
  cat("\n======== Model:", ch, "========\n")
  models[[ch]] <- lmer(
    lme_formula,
    data = d,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )
  print(summary(models[[ch]]))
}

pred_all <- bind_rows(lapply(levels(matched_long$Cohort), function(ch) {
  d <- matched_long %>% filter(Cohort == ch) %>% droplevels()
  grid <- expand.grid(
    Age_Centered = seq(-WINDOW_YRS, WINDOW_YRS, by = 0.1),
    Group = factor(c("Control", "Cancer Case"), levels = c("Control", "Cancer Case"))
  ) %>%
    mutate(
      Age_Post = if_else(Age_Centered >= 0, Age_Centered, 0),
      index_age_z = 0,
      base_race = factor(levels(d$base_race)[1], levels = levels(d$base_race)),
      base_marital = factor(levels(d$base_marital)[1], levels = levels(d$base_marital)),
      base_pckgr = factor(levels(d$base_pckgr)[1], levels = levels(d$base_pckgr)),
      Cohort = ch
    )
  grid$pred <- predict(models[[ch]], newdata = grid, re.form = NA)
  grid
})) %>%
  mutate(Cohort = factor(Cohort, levels = levels(matched_long$Cohort)))

p_traj <- ggplot(pred_all, aes(x = Age_Centered, y = pred, color = Group)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  facet_wrap(~ Cohort) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Frailty Trajectories: Separate Risk-Set Matching by Burden Cohort",
    subtitle = "Risk-set cohorts generated by 2.1; centered on own attained age at index",
    x = "Years relative to index",
    y = "Predicted frailty index (no-cancer FI)",
    color = NULL
  ) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(-WINDOW_YRS, WINDOW_YRS, by = 4))

print(p_traj)

ggsave(file.path(results_dir, "6.1_riskset_trajectories_by_cohort.png"),
       p_traj, width = 10, height = 5, dpi = 300)
ggsave(file.path(results_dir, "6.1_riskset_index_age_density.png"),
       p_index, width = 10, height = 5, dpi = 300)

cat("\nSaved figures to Results/cancer/\n")
