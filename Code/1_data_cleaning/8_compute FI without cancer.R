# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 8_compute FI without cancer.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Computes a frailty index that excludes cancer-related deficit items to avoid circularity when cancer is the exposure of interest.
# ==============================================================================

library(dplyr)
library(tidyr)

# ==============================================================================
# 1. LOAD YOUR EXISTING DATA (which already has your standard FI)
# ==============================================================================
target_dir <- "/n/home06/xyzhou/Frailty"
fi_data <- readRDS(file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds"))

# Define the 27 FI items (EXCLUDING CANCER)
fi_items_nocancer <- c(
  "fi_chol", "fi_hbp", "fi_clau", "fi_pe", "fi_visual", "fi_cad", "fi_cvd", 
  "fi_gout", "fi_hfra", "fi_ost", "fi_arth", "fi_dvrt", "fi_ulcer", "fi_ucol", 
  "fi_park", "fi_db", "fi_asthma", 
  "fi_antidepressant", "fi_tranq", "fi_polypharmacy", 
  "fi_stairs", "fi_balance", "fi_pace", 
  "fi_act_score", "fi_tv_score", "fi_abnormal_bmi_score", "fi_wtloss_score"
)

# ==============================================================================
# 2. GET DENOMINATOR FOR NON-CANCER ITEMS ONLY
# ==============================================================================
cycle_denominators_nocancer <- fi_data %>%
  group_by(cycle) %>%
  summarize(across(all_of(fi_items_nocancer), ~ if_else(all(is.na(.)), 0, 1)), .groups = "drop") %>%
  rowwise() %>%
  mutate(items_asked_nocancer = sum(c_across(all_of(fi_items_nocancer)))) %>%
  select(cycle, items_asked_nocancer) %>%
  ungroup()

# ==============================================================================
# 3. APPEND NON-CANCER FI TO EXISTING DATA
# ==============================================================================
fi_final <- fi_data %>%
  left_join(cycle_denominators_nocancer, by = "cycle") %>%
  mutate(
    # Only calculate the non-cancer metrics
    sum_deficits_nocancer = rowSums(select(., all_of(fi_items_nocancer)), na.rm = TRUE),
    n_answered_nocancer = rowSums(!is.na(select(., all_of(fi_items_nocancer)))),
    n_missing_nocancer = items_asked_nocancer - n_answered_nocancer,
    
    # Compute the no-cancer index
    fi_score_nocancer = if_else(n_answered_nocancer == 0, NA_real_, sum_deficits_nocancer / n_answered_nocancer),
    
    # Categorize the no-cancer index
    frailty_cat_nocancer = case_when(
      is.na(fi_score_nocancer) ~ NA_character_,
      fi_score_nocancer < 0.1 ~ "Non-frail",
      fi_score_nocancer >= 0.1 & fi_score_nocancer < 0.2 ~ "Pre-frail",
      fi_score_nocancer >= 0.2 & fi_score_nocancer < 0.3 ~ "Mildly frail",
      fi_score_nocancer >= 0.3 ~ "Mod/Sev frail"
    )
  )

# ==============================================================================
# 4. INSPECT AND SAVE
# ==============================================================================
print(
  fi_final %>% 
    filter(cycle == "12") %>% 
    select(id, cycle, fi_score, fi_score_nocancer) %>% 
    head()
)

# Overwrite the file with the newly added columns
saveRDS(fi_final, file.path("/n/home06/xyzhou/Frailty/2_data analysis/before and after cancer/data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"))
cat("\nAppended fi_score_nocancer and saved successfully.\n")
