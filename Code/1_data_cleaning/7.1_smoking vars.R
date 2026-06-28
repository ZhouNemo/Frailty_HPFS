# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 7.1_smoking vars.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Merges smoking variables, including smoking status and pack-year measures, into the imputed longitudinal frailty dataset and checks merge quality.
# ==============================================================================

# Load the dplyr package
library(dplyr)

FI_longitudinal_1988_2020 <- readRDS("~/Frailty/FI_longitudinal_1988_2020.rds")
FI_longitudinal_1986_2020_IMPUTED <- readRDS("~/Frailty/FI_longitudinal_1986_2020_IMPUTED.rds")

# Merge the datasets
FI_longitudinal_1986_2020_IMPUTED <- FI_longitudinal_1986_2020_IMPUTED %>%
  left_join(
    # Select only the keys and the 3 variables you want from the source dataset
    # This prevents duplicating other columns (which creates .x and .y columns)
    FI_longitudinal_1988_2020 %>% select(id, cycle, smoke, pckyr, pckgr), 
    
    # Specify the columns that link the two datasets together
    by = c("id", "cycle") 
  )


# 1. Check if the new columns exist in the dataset
# This will print TRUE for each variable if it was successfully added
c("smoke.y", "pckyr.y", "pckgr.y") %in% names(FI_longitudinal_1986_2020_IMPUTED)

# 2. Check the summary of the newly added columns
# This shows you the min, max, mean, and how many NAs are in these columns. 
# If they are 100% NA, the merge keys (e.g., ID or Year) didn't match properly.
summary(FI_longitudinal_1986_2020_IMPUTED[c("smoke.y", "pckyr.y", "pckgr.y")])

# 3. Check for row duplication
# The number of rows in the imputed dataset SHOULD NOT change after a left_join.
# If nrow() is much larger than you expect, you have duplicate ID/Year combinations 
# in the FI_longitudinal_1988_2020 dataset.
nrow(FI_longitudinal_1986_2020_IMPUTED)

# 4. Spot check a few rows 
# Look at a few rows where we know the person has data in the new columns
FI_longitudinal_1986_2020_IMPUTED %>%
  select(id, cycle, smoke.y, pckyr.y, pckgr.y) %>%
  filter(!is.na(smoke.y) | !is.na(pckyr.y)) %>%
  head(40)

# Save updated dataset with a distinct '_COVARIATES' suffix
save_path_rds <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds")
save_path_csv <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.csv")

saveRDS(fi_long_merged, file = save_path_rds)
# write.csv(fi_long_merged, file = save_path_csv, row.names = FALSE)

cat("\nDataset merged and saved to:", save_path_rds, "\n")
