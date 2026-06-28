# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 6_compute FI.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Computes the longitudinal frailty index from imputed frailty deficit items using cycle-specific denominators and produces the scored frailty-index dataset.
# ==============================================================================

library(dplyr)
library(tidyr)

# ==============================================================================
# 1. LOAD THE IMPUTED LONGITUDINAL DATA
# ==============================================================================
target_dir <- "/n/home06/xyzhou/Frailty"
print("Loading imputed dataset...")
fi_imputed <- readRDS(file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds"))

# Define the 28 FI items (with 'fi_' prefix from our renaming step)
fi_items <- c(
  "fi_cancer", "fi_chol", "fi_hbp", "fi_clau", "fi_pe", "fi_visual", "fi_cad", "fi_cvd", 
  "fi_gout", "fi_hfra", "fi_ost", "fi_arth", "fi_dvrt", "fi_ulcer", "fi_ucol", 
  "fi_park", "fi_db", "fi_asthma", 
  "fi_antidepressant", "fi_tranq", "fi_polypharmacy", 
  "fi_stairs", "fi_balance", "fi_pace", 
  "fi_act_score", "fi_tv_score", "fi_abnormal_bmi_score", "fi_wtloss_score"
)

# ==============================================================================
# 2. DYNAMIC DENOMINATOR (ADAPTIVE TO MISSING CYCLES)
# ==============================================================================
# Determine how many items were actually measurable in each cycle. 
# (If an item is 100% NA, it wasn't asked, so it shouldn't penalize the participant)
cycle_denominators <- fi_imputed %>%
  group_by(cycle) %>%
  summarize(across(all_of(fi_items), ~ if_else(all(is.na(.)), 0, 1)), .groups = "drop") %>%
  rowwise() %>%
  mutate(items_asked = sum(c_across(all_of(fi_items)))) %>%
  select(cycle, items_asked) %>%
  ungroup()

# ==============================================================================
# 3. COMPUTE FI SCORE & CATEGORIES
# ==============================================================================
print("Computing Frailty Index & Categories...")

fi_final <- fi_imputed %>%
  # Only compute FI for participants who actually returned the survey that year
  filter(participated == 1) %>%
  left_join(cycle_denominators, by = "cycle") %>%
  mutate(
    # A. Sum of deficits present (Numerator)
    sum_deficits = rowSums(select(., all_of(fi_items)), na.rm = TRUE),
    
    # B. Count how many items this specific person actually answered
    n_answered = rowSums(!is.na(select(., all_of(fi_items)))),
    
    # C. Count how many items they skipped (relative to what was asked that year)
    n_missing = items_asked - n_answered,
    
    # D. Compute Index FOR EVERYONE 
    # (Only returns NA if they answered exactly 0 items to prevent division by zero)
    fi_score = if_else(n_answered == 0, NA_real_, sum_deficits / n_answered),
    
    # E. Assign Frailty Categories based on the unrestricted score
    frailty_cat = case_when(
      is.na(fi_score) ~ NA_character_,
      fi_score < 0.1 ~ "Non-frail",
      fi_score >= 0.1 & fi_score < 0.2 ~ "Pre-frail",
      fi_score >= 0.2 & fi_score < 0.3 ~ "Mildly frail",
      fi_score >= 0.3 ~ "Mod/Sev frail"
    )
  ) 
# NOTE: I removed the select() line that dropped items_asked and n_answered.
# Now, these variables will remain in fi_final so you can use them later!

# ==============================================================================
# 4. INSPECT & SAVE
# ==============================================================================
# Check the first few rows for 2012 to ensure it worked
cat("\n--- Sample Output (Cycle 2012) ---\n")
print(
  fi_final %>% 
    filter(cycle == "12") %>% 
    select(id, cycle, sum_deficits, items_asked, n_answered, n_missing, fi_score) %>% 
    head()
)

# Save updated, scored dataset
save_path_rds <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds")
save_path_csv <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.csv")

saveRDS(fi_final, file = save_path_rds)
# write.csv(fi_final, file = save_path_csv, row.names = FALSE)

cat("\nFrailty Index computed successfully and saved to:", save_path_rds, "\n")