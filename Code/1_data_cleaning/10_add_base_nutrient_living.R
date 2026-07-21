# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health
#          Professionals Follow-up Study
# Script: 10_add_base_nutrient_living.R
# Author: Nemo Zhou
# Date started: 2026-07-05
# Date last updated: 2026-07-05
# Purpose:
#   Adds baseline cycle-86 nutritional covariates and living arrangement to the
#   canonical cancer analytic longitudinal dataset. The script derives one
#   participant-level baseline row from cycle 86, preserves missing nutritional
#   values as NA, recodes missing baseline living arrangement as "Missing", and
#   joins the resulting base_* covariates back to every participant-cycle row.
# ==============================================================================

library(dplyr)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Data")
input_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")

baseline_vars <- c(
  "base_calor",
  "base_sat",
  "base_diet_chol",
  "base_nblnk",
  "base_alco",
  "base_living"
)

required_source_vars <- c(
  "id",
  "cycle",
  "calor",
  "sat",
  "diet_chol",
  "nblnk",
  "alco",
  "living_arr"
)

print("Loading canonical cancer analytic dataset...")
fi_long <- readRDS(input_path)

missing_source_vars <- setdiff(required_source_vars, names(fi_long))
if (length(missing_source_vars) > 0) {
  stop(
    "Missing required source variables in input dataset: ",
    paste(missing_source_vars, collapse = ", ")
  )
}

baseline_86 <- fi_long %>%
  filter(as.character(cycle) == "86")

duplicate_baseline_ids <- baseline_86 %>%
  count(id, name = "baseline_rows") %>%
  filter(baseline_rows > 1)

if (nrow(duplicate_baseline_ids) > 0) {
  stop(
    "Cycle 86 contains duplicate participant rows. Example duplicated ids: ",
    paste(head(duplicate_baseline_ids$id, 10), collapse = ", ")
  )
}

base_covariates <- baseline_86 %>%
  transmute(
    id = as.character(id),
    base_calor = as.numeric(as.character(calor)),
    base_sat = as.numeric(as.character(sat)),
    base_diet_chol = as.numeric(as.character(diet_chol)),
    base_nblnk = as.numeric(as.character(nblnk)),
    base_alco = as.numeric(as.character(alco)),
    base_living = factor(if_else(is.na(living_arr), "Missing", as.character(living_arr)))
  )

fi_long_with_base <- fi_long %>%
  mutate(id = as.character(id)) %>%
  select(-any_of(baseline_vars)) %>%
  left_join(base_covariates, by = "id")

print("Baseline covariate verification:")
cat("Total rows:", nrow(fi_long_with_base), "\n")
cat("Distinct participants:", n_distinct(fi_long_with_base$id), "\n")
cat("Baseline cycle 86 rows:", nrow(baseline_86), "\n")
cat("Baseline cycle 86 distinct participants:", n_distinct(baseline_86$id), "\n")
cat("Duplicate baseline IDs:", nrow(duplicate_baseline_ids), "\n")

missing_summary <- fi_long_with_base %>%
  filter(as.character(cycle) == "86") %>%
  summarize(
    across(
      all_of(baseline_vars),
      list(
        missing_n = ~ sum(is.na(.x)),
        missing_pct = ~ round(mean(is.na(.x)) * 100, 2)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_pattern = "(.+)_(missing_n|missing_pct)$"
  ) %>%
  arrange(match(variable, baseline_vars))

print(missing_summary)

base_living_missing_category <- sum(
  fi_long_with_base$cycle == "86" &
    !is.na(fi_long_with_base$base_living) &
    fi_long_with_base$base_living == "Missing"
)
cat(
  "Baseline base_living values recoded to \"Missing\":",
  base_living_missing_category,
  "\n"
)

saveRDS(fi_long_with_base, input_path)
cat("\nSaved baseline covariates to:", input_path, "\n")
