# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 3_rename vars.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Standardizes frailty-item variable names across cycles, applies the fi_ prefix convention, and prepares a consistent long-format item structure for imputation and scoring.
# ==============================================================================

library(dplyr)

fi <- readRDS("/n/home06/xyzhou/Frailty/FI_longitudinal_1986_2020.rds")

# 1. Define the 28 Frailty Deficits
# Note: 'cancer' covers bladder and kidney as merged in the previous extraction step.
fi_items <- c(
  "cancer", "chol", "hbp", "clau", "pe", "visual", "cad", "cvd", 
  "gout", "hfra", "ost", "arth", "dvrt", "ulcer", "ucol", 
  "park", "db", "asthma", 
  "antidepressant", "tranq", "polypharmacy", 
  "stairs", "balance", "pace", 
  "act_score", "tv_score", "abnormal_bmi_score", "wtloss_score"
)

# 2. UPDATE: All 18 Biennial Cycles
cycles <- c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")

# 3. Generate the expected names WITH underscores
expected_cols <- as.vector(outer(fi_items, cycles, paste, sep = "_"))

# 4. Find the exact overlap between our expected names and the actual dataset
cols_that_exist <- intersect(expected_cols, names(fi))

# 5. Only attempt the rename if R actually found matching columns!
if (length(cols_that_exist) > 0) {
  
  fi <- fi %>%
    rename_with(~ paste0("fi_", .x), all_of(cols_that_exist))
  
  print(paste("Success! Renamed", length(cols_that_exist), "variables."))
  
  # Print the first few to verify they look correct
  print(head(grep("^fi_", names(fi), value = TRUE)))
  
} else {
  
  print("WARNING: R found 0 matching columns. The items in the dataset might not have underscores.")
  print("Here are some of the columns that currently end in 88 so we can check the spelling:")
  print(head(grep("88$", names(fi), value = TRUE), 10))
  
}

# ==============================================================================
# Final Cleanup: Rename all '_XX' variables to 'XX'
# ==============================================================================

fi <- fi %>%
  rename_with(
    # .fn tells it HOW to rename: 
    # gsub looks for an underscore followed by 2 digits at the end of the string ("$")
    # and replaces it with just the 2 digits ("\\1")
    .fn = ~ gsub("_(\\d{2})$", "\\1", .x), 
    
    # .cols tells it WHICH columns to apply this to:
    # Only columns that match our "_XX" pattern at the very end
    .cols = matches("_\\d{2}$")
  )

# Quick check to verify the names changed correctly!
print("Column name check (should show fi_hbp86, fi_tv_score88, etc.):")
head(grep("\\d{2}$", names(fi), value = TRUE), 15)