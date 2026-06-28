# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 7_merge covariate.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Merges external covariates, death data, race, smoking, and other participant-level or cycle-level variables into the scored longitudinal frailty-index dataset.
# ==============================================================================

library(dplyr)
library(tidyr)
library(hpfs)
library(stringr)

# Load the scored longitudinal Frailty Index dataset
target_dir <- "/n/home06/xyzhou/Frailty"
fi_long <- readRDS(file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds"))

print("Loading covariate datasets...")
deaths_df <- load_hpfs_deaths()
race_df <- load_hpfs_race()
smoking_df <- load_hpfs_smoking()
# nutrients_df <- load_hpfs_nutrients()

# ==============================================================================
# 1. MERGE EXTERNAL COVARIATES AND DEATH DATA
# ==============================================================================

print("Merging all covariates and death data...")

fi_long_merged <- fi_long %>%
  
  # A. ADD RACE (Convert ID to character on the fly)
  left_join(race_df %>% mutate(id = as.character(id)), by = "id") %>%
  
  # B. ADD SMOKING (Convert ID and Cycle to character just in case)
  left_join(smoking_df %>% mutate(id = as.character(id), 
                                  cycle = str_pad(as.character(cycle), 2, pad = "0")), 
            by = c("id", "cycle")) %>%
  
  # C. ADD NUTRIENTS (Convert ID and Cycle)
  # left_join(nutrients_df %>% mutate(id = as.character(id), 
  #                                   cycle = str_pad(as.character(cycle), 2, pad = "0")), 
  #           by = c("id", "cycle")) %>%
  
  # D. ADD DEATH DATA (Convert ID)
  left_join(deaths_df %>% mutate(id = as.character(id)), by = "id")

# ==============================================================================
# 2. HARMONIZE AND EXTRACT CROSS-SECTIONAL COVARIATES
# ==============================================================================
print("Extracting and harmonizing cross-sectional covariates...")

# Expanded to all 18 biennial cycles
cycles <- c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")
extracted_list <- list()

for (cyc in cycles) {
  df_name <- paste0("hp", cyc)
  
  if (df_name %in% names(hp_data_list)) {
    df <- hp_data_list[[df_name]]
    cols <- names(df)
    
    # Identify exact column names for this cycle
    mar_col  <- intersect(paste0(c("mar", "marital"), cyc), cols)[1]
    work_col <- intersect(paste0(c("work", "wrkst"), cyc), cols)[1]
    liv_col  <- intersect(paste0(c("livng", "livar", "living"), cyc), cols)[1]
    
    c_wife   <- paste0("lwife", cyc)
    c_alone  <- paste0("lalon", cyc)
    c_nurs   <- paste0("lnurs", cyc)
    c_aslv   <- paste0("laslv", cyc)
    c_senior <- paste0("lsenior", cyc)
    c_fam    <- paste0("lfam", cyc)
    
    # Build a standardized temporary dataframe. 
    # If a variable doesn't exist in this cycle, it safely becomes NA.
    temp_df <- data.frame(
      id = as.character(df$id),
      cycle = as.character(cyc),
      raw_mar  = if(!is.na(mar_col)) df[[mar_col]] else NA,
      raw_work = if(!is.na(work_col)) df[[work_col]] else NA,
      raw_liv  = if(!is.na(liv_col)) df[[liv_col]] else NA,
      raw_wife = if(c_wife %in% cols) df[[c_wife]] else NA,
      raw_alone = if(c_alone %in% cols) df[[c_alone]] else NA,
      raw_nurs = if(c_nurs %in% cols) df[[c_nurs]] else NA,
      raw_aslv = if(c_aslv %in% cols) df[[c_aslv]] else NA,
      raw_senior = if(c_senior %in% cols) df[[c_senior]] else NA,
      raw_fam  = if(c_fam %in% cols) df[[c_fam]] else NA,
      stringsAsFactors = FALSE
    )
    
    # Safely harmonize the values using the standardized columns
    temp_df <- temp_df %>%
      mutate(
        marital_status = case_when(
          is.na(raw_mar) ~ NA_character_,
          # 1988 Coding (Assuming 86 follows the same if present)
          cycle %in% c("86", "88") & raw_mar == 1 ~ "Married",
          cycle %in% c("86", "88") & raw_mar == 2 ~ "Widowed",
          cycle %in% c("86", "88") & raw_mar == 3 ~ "Never Married",
          cycle %in% c("86", "88") & raw_mar == 4 ~ "Divorced/Separated",
          # 1990 - 2020 Coding
          !cycle %in% c("86", "88") & raw_mar == 1 ~ "Married",
          !cycle %in% c("86", "88") & raw_mar == 2 ~ "Divorced/Separated",
          !cycle %in% c("86", "88") & raw_mar == 3 ~ "Widowed",
          !cycle %in% c("86", "88") & raw_mar == 4 ~ "Never Married",
          TRUE ~ NA_character_
        ),
        
        work_status = case_when(
          is.na(raw_work) ~ NA_character_,
          raw_work == 1 ~ "Full-time",
          raw_work == 2 ~ "Part-time",
          raw_work == 3 ~ "Retired",
          raw_work == 4 ~ "Disabled",
          # Kept original logic, expanded safely for missing unemployed tags
          raw_work == 5 ~ "Unemployed", 
          TRUE ~ NA_character_
        ),
        
        living_arr = case_when(
          # Pre-2004: 1986/1988 Mapping
          !is.na(raw_liv) & cycle %in% c("86", "88") & raw_liv == 1 ~ "With spouse/partner",
          !is.na(raw_liv) & cycle %in% c("86", "88") & raw_liv == 2 ~ "Alone",
          !is.na(raw_liv) & cycle %in% c("86", "88") & raw_liv == 3 ~ "With other family",
          !is.na(raw_liv) & cycle %in% c("86", "88") & raw_liv == 4 ~ "Other",
          
          # Pre-2004: 1990-2002 Mapping
          !is.na(raw_liv) & !cycle %in% c("86", "88") & raw_liv == 1 ~ "Alone",
          !is.na(raw_liv) & !cycle %in% c("86", "88") & raw_liv == 2 ~ "With spouse/partner",
          !is.na(raw_liv) & !cycle %in% c("86", "88") & raw_liv == 3 ~ "With other family",
          !is.na(raw_liv) & !cycle %in% c("86", "88") & raw_liv == 4 ~ "Nursing home/Assisted Living",
          !is.na(raw_liv) & !cycle %in% c("86", "88") & raw_liv == 5 ~ "Other",
          
          # 2004+: Checkbox Mapping (Using hierarchy if multiple are checked)
          is.na(raw_liv) & raw_nurs == 1 ~ "Nursing home/Assisted Living",
          is.na(raw_liv) & raw_aslv == 1 ~ "Nursing home/Assisted Living",
          is.na(raw_liv) & raw_senior == 1 ~ "Nursing home/Assisted Living",
          is.na(raw_liv) & raw_wife == 1 ~ "With spouse/partner",
          is.na(raw_liv) & raw_fam == 1 ~ "With other family",
          is.na(raw_liv) & raw_alone == 1 ~ "Alone",
          
          TRUE ~ NA_character_
        )
      ) %>%
      select(id, cycle, marital_status, work_status, living_arr)
    
    extracted_list[[cyc]] <- temp_df
  }
}

# Combine all harmonized years
harmonized_covariates <- bind_rows(extracted_list)

# ==============================================================================
# 3. FINAL MERGE & SAVE
# ==============================================================================
print("Merging cleaned text covariates into the main dataset...")

fi_long_merged <- fi_long_merged %>%
  left_join(harmonized_covariates, by = c("id", "cycle"))

# Verify a few rows
print(
  fi_long_merged %>%
    filter(cycle %in% c("86", "88", "96", "04", "16")) %>%
    select(id, cycle, marital_status, living_arr, work_status) %>%
    drop_na(marital_status) %>%
    group_by(cycle) %>%
    slice(1)
)

# Save updated dataset with a distinct '_COVARIATES' suffix
save_path_rds <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.rds")
save_path_csv <- file.path(target_dir, "FI_longitudinal_1986_2020_IMPUTED.csv")

saveRDS(fi_long_merged, file = save_path_rds)
# write.csv(fi_long_merged, file = save_path_csv, row.names = FALSE)

cat("\nDataset merged and saved to:", save_path_rds, "\n")