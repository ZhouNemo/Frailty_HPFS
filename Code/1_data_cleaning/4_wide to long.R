# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 4_wide to long.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Converts frailty and questionnaire variables from wide cycle-specific format into long person-cycle format and merges follow-up participation information for longitudinal analysis.
# ==============================================================================

library(dplyr)
library(tidyr)
library(stringr)

# Load Baseline and Follow-up
baseline_df <- load_hpfs_baseline()
followup_df <- load_hpfs_followup()

# ==============================================================================
# 0. DEFINE THE TARGET LONGITUDINAL CYCLES
# ==============================================================================
# UPDATE: Expanded to all 18 biennial cycles
cycles_target <- c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")

# ==============================================================================
# 1. PRE-PROCESS FOLLOW-UP DATA
# ==============================================================================
print("Processing Follow-up data...")

followup_clean <- followup_df %>%
  mutate(
    # Standardize ID to Character to match main 'fi' dataframe
    id = as.character(id),
    
    # Standardize Cycle to 2-digit character string (e.g., forces '0' to '00')
    cycle = str_pad(as.character(cycle), 2, pad = "0"),
    
    # Define Binary Participation Indicator (rtmnyr > 0 means they returned it)
    participated = if_else(rtmnyr > 0, 1, 0)
  ) %>%
  # Filter only to the chronological cycles we care about
  filter(cycle %in% cycles_target) %>%
  select(id, cycle, worked_rtmnyr = rtmnyr, participated) %>%
  
  # CRITICAL: Deduplicate. If a person sent two forms for one cycle, we just need one participation flag.
  distinct(id, cycle, .keep_all = TRUE)

# ==============================================================================
# 2. PRE-PROCESS BASELINE DATA
# ==============================================================================
print("Processing Baseline data...")

baseline_clean <- baseline_df %>%
  mutate(id = as.character(id)) %>%
  # Deduplicate just to be safe
  distinct(id, .keep_all = TRUE)

# ==============================================================================
# 3. PIVOT FRAILTY DATA TO LONG FORMAT
# ==============================================================================
print("Pivoting main 'fi' dataset (Wide to Long)...")

# Crucial step: Make IDs character before the big pivot
fi_unique <- fi %>%
  mutate(id = as.character(id)) %>%
  distinct(id, .keep_all = TRUE) 

fi_long <- fi_unique %>%
  pivot_longer(
    # Select cols ending in target years (handles 'fi_hbp88', 'fi_tv_score00', etc.)
    cols = ends_with(cycles_target),
    
    # Split: Part 1 goes to separate cols (.value), Part 2 (year) goes to 'cycle'
    names_to = c(".value", "cycle"),
    
    # Regex: (.*) = any chars before year, (\\d{2}) = last 2 digits.
    names_pattern = "(.*)(\\d{2})$"
  )

# ==============================================================================
# 4. MEMORY MANAGEMENT (Prevent Crash)
# ==============================================================================
# Remove the large original dataframes to free up RAM before merging
print("Cleaning up memory for the big merge...")
rm(fi_unique, followup_df, baseline_df)
gc()

# ==============================================================================
# 5. FINAL MERGE: Combine Everything
# ==============================================================================
print("Performing final longitudinal merge...")

# Use fi_long as the base to duplicate static baseline variables appropriately
fi_long_final <- fi_long %>%
  # Join Baseline (Static info joins by ID only)
  left_join(baseline_clean, by = "id") %>%
  # Join Follow-up Status (Time-varying info joins by both ID and Cycle)
  left_join(followup_clean, by = c("id", "cycle")) %>%
  
  # Organize: Sort by person, then timeline
  arrange(id, cycle)

# Optional memory cleanup
rm(fi_long, baseline_clean, followup_clean)
gc()

# ==============================================================================
# 6. VERIFY SUCCESS
# ==============================================================================
print(paste("Data unified! Total Rows:", nrow(fi_long_final), 
            "| Total Columns:", ncol(fi_long_final)))

# Look at the structure to confirm it merged correctly
head(fi_long_final)


# ==============================================================================
# SAVE OUTPUT
# ==============================================================================
target_dir <- "/n/home06/xyzhou/Frailty"

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

write.csv(fi_long_final, file = file.path(target_dir, "FI_longitudinal_1986_2020.csv"), row.names = FALSE)
saveRDS(fi_long_final, file = file.path(target_dir, "FI_longitudinal_1986_2020.rds"))

cat("Files successfully saved to:", target_dir)