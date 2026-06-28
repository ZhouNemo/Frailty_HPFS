# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 7.2_nutrient vars.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Merges nutrient variables into the imputed longitudinal frailty dataset, including calories, saturated fat, dietary cholesterol, blank FFQ items, and alcohol intake, with basic quality checks.
# ==============================================================================

# Load the dplyr package
library(dplyr)

# 1. Load your main imputed dataset
FI_longitudinal_1986_2020_IMPUTED <- readRDS("~/Frailty/FI_longitudinal_1986_2020_IMPUTED.rds")

# 2. Generate the nutrition dataset using your newly fixed function
# (Ensure the get_hpfs_nutrients() function is run in your environment first)
#nutrients_data <- get_hpfs_nutrients()
# 3. Merge the datasets
fi_long_merged <- FI_longitudinal_1986_2020_IMPUTED %>%
  left_join(
    # Convert 'id' to character, RENAME 'chol', then select
    nutrients_data %>% 
      mutate(id = as.character(id)) %>% 
      rename(diet_chol = chol) %>%  # <--- Renaming step here
      select(id, cycle, calor, sat, diet_chol, nblnk, alco), 
    by = c("id", "cycle") 
  )

# ---------------------------------------------------------
# QUALITY CHECKS (Updated with diet_chol)
# ---------------------------------------------------------

# Check 1: Check if the new columns exist
c("calor", "sat", "diet_chol", "nblnk", "alco") %in% names(fi_long_merged)

# Check 2: Check the summary
summary(fi_long_merged[c("calor", "sat", "diet_chol", "nblnk", "alco")])

# Check 3: Spot check a few rows 
fi_long_merged %>%
  select(id, cycle, calor, sat, diet_chol, nblnk, alco) %>%
  filter(!is.na(calor) | !is.na(alco)) %>%
  head(40)

# ---------------------------------------------------------
# SAVE DATA
# ---------------------------------------------------------

# Make sure to define target_dir if it isn't already in your environment
target_dir <- "~/Frailty"

# Save updated dataset with a distinct suffix so you don't overwrite the original
save_path_rds <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds")
# save_path_csv <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.csv")

saveRDS(fi_long_merged, file = save_path_rds)
# write.csv(fi_long_merged, file = save_path_csv, row.names = FALSE)

cat("\nDataset merged and saved to:", save_path_rds, "\n")