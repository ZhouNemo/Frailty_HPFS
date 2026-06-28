# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 5_logical imputation.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Applies logical imputation rules to longitudinal frailty deficit items, especially irreversible disease indicators, while distinguishing missingness from unasked or inapplicable questionnaire items.
# ==============================================================================

library(dplyr)
library(tidyr)

# ==============================================================================
# 0. LOAD DATA & SETUP ITEMS
# ==============================================================================
target_dir <- "/n/home06/xyzhou/Frailty"

print("Loading long-format dataset...")
fi_long_final <- readRDS(file.path(target_dir, "FI_longitudinal_1986_2020.rds"))

# The full list of 28 FI items
fi_items <- c(
  "fi_cancer", "fi_chol", "fi_hbp", "fi_clau", "fi_pe", "fi_visual", "fi_cad", "fi_cvd", 
  "fi_gout", "fi_hfra", "fi_ost", "fi_arth", "fi_dvrt", "fi_ulcer", "fi_ucol", 
  "fi_park", "fi_db", "fi_asthma", 
  "fi_antidepressant", "fi_tranq", "fi_polypharmacy", 
  "fi_stairs", "fi_balance", "fi_pace", 
  "fi_act_score", "fi_tv_score", "fi_abnormal_bmi_score", "fi_wtloss_score"
)

# Diseases that are irreversible (Excluding 'fi_cancer' per instructions)
disease_vars <- c(
  "fi_chol", "fi_hbp", "fi_clau", "fi_pe", "fi_visual", "fi_cad", "fi_cvd", 
  "fi_gout", "fi_hfra", "fi_ost", "fi_arth", "fi_dvrt", "fi_ulcer", "fi_ucol", 
  "fi_park", "fi_db", "fi_asthma"
)

# ==============================================================================
# 1. FIX CHRONOLOGICAL SORTING & APPLY IMPUTATIONS
# ==============================================================================
print("Applying logical imputations...")

fi_imputed <- fi_long_final %>%
  
  # --- NEW FIX: Force participated = 1 for the 1986 baseline ---
  mutate(
    participated = if_else(cycle == "86", 1, participated)
  ) %>%
  
  # --- FIX: Create a true numeric year for chronological sorting ---
  mutate(
    year_num = case_when(
      as.numeric(cycle) >= 86 ~ 1900 + as.numeric(cycle),
      as.numeric(cycle) <= 20 ~ 2000 + as.numeric(cycle),
      TRUE ~ NA_real_
    )
  ) %>%
  
  # CRITICAL: Ensure the data is strictly ordered by person, then by true year
  arrange(id, year_num) %>%
  group_by(id) %>%
  
  # ----------------------------------------------------------------------------
# RULE 1: Within-cycle Physical Activity Imputation
# If act_score is missing, but they are explicitly limited in balance, pace, or stairs -> 1.0
# ----------------------------------------------------------------------------
mutate(
  fi_act_score = case_when(
    is.na(fi_act_score) & (fi_balance %in% 1 | fi_pace %in% 1 | fi_stairs %in% 1) ~ 1.0,
    TRUE ~ fi_act_score
  )
) %>%
  
  # ----------------------------------------------------------------------------
# RULE 2: 1-Cycle (2-Year) Last Observation Carried Forward (LOCF) for ALL variables
# 'lag(.)' looks back exactly one cycle. It will NOT cascade indefinitely.
# ----------------------------------------------------------------------------
mutate(across(all_of(fi_items), ~ ifelse(is.na(.) & !is.na(lag(.)), lag(.), .))) %>%
  
  # ----------------------------------------------------------------------------
# RULE 3: Once a disease, always a disease (Strict participation logic)
# This forces a 1.0 forward indefinitely for chronic diseases, but ONLY if the 
# participant actually returned the survey (participated == 1) that year.
# ----------------------------------------------------------------------------
mutate(across(all_of(disease_vars), ~ {
  
  # Calculate if they have EVER had the disease (treating NAs as 0 so cummax doesn't break)
  ever_had <- cummax(ifelse(is.na(.), 0, .))
  
  # Force to 1.0 ONLY IF they've had it before AND they actually participated this cycle.
  # Otherwise, leave it as it was (which preserves NAs for missed survey years).
  ifelse(ever_had == 1 & participated == 1, 1.0, .)
  
})) %>%
  
  ungroup() %>%
  
  # Optional: Remove the temporary year_num column if you don't need it later
  select(-year_num)

print("Imputations complete!")

# ==============================================================================
# 2. SAVE THE IMPUTED DATASET
# ==============================================================================
save_path_rds <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds")
save_path_csv <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.csv")

saveRDS(fi_imputed, file = save_path_rds)
# write.csv(fi_imputed, file = save_path_csv, row.names = FALSE)

cat("Imputed dataset successfully saved to:", save_path_rds, "\n")
