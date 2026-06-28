# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 9_time variables generation.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Generates demographic, questionnaire-cycle, death, and cancer-timing variables used to align frailty measurements before and after cancer diagnosis.
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)

# Load the merged cancer analytic dataset from the project Data folder.
project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Data")
input_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
fi_long_merged <- readRDS(input_path)

# ==============================================================================
# 1. CALCULATE DEMOGRAPHIC & CLINICAL TIMELINES (FULL COHORT)
# ==============================================================================
print("Calculating demographic and cancer timelines for the full cohort...")

fi_time_prepped <- fi_long_merged %>%
  mutate(
    # 0. FIX DATA TYPES: Convert time variables from Character to Numeric
    # (Note: R might throw a warning about "NAs introduced by coercion" 
    # if it turns blank spaces into NAs. That is exactly what we want!)
    cancer_dateca = as.numeric(as.character(cancer_dateca)),
    dtdth         = as.numeric(as.character(dtdth)),
    dbmy09        = as.numeric(as.character(dbmy09)),
    worked_rtmnyr = as.numeric(as.character(worked_rtmnyr)), # Included just in case!
    participated  = as.numeric(as.character(participated)),
    fi_score_nocancer = as.numeric(as.character(fi_score_nocancer)),
    
    # 1. Create readable Dates (using 15th as mid-month)
    dob = make_date(year = 1900 + floor((dbmy09 - 1) / 12), 
                    month = ((dbmy09 - 1) %% 12) + 1, 
                    day = 15),
    
    dod = make_date(year = 1900 + floor((dtdth - 1) / 12), 
                    month = ((dtdth - 1) %% 12) + 1, 
                    day = 15),
    
    cancer_date = make_date(year = 1900 + floor((cancer_dateca - 1) / 12), 
                            month = ((cancer_dateca - 1) %% 12) + 1, 
                            day = 15),
    
    # 2. Calculate Ages at key events (will be NA if event hasn't occurred)
    age_at_death = (dtdth - dbmy09) / 12,
    age_at_cancer = (cancer_dateca - dbmy09) / 12,
    
    # 3. Status Indicators 
    death_status = ifelse(!is.na(dtdth), 1, 0)
  )
# ==============================================================================
# 2. CREATE CYCLE-SPECIFIC TIME VARIABLES
# ==============================================================================
print("Calculating trajectory timescales (Age at cycle, Time to Event)...")

fi_trajectory <- fi_time_prepped %>%
  # Filter to rows where the participant returned the survey and the cancer-free
  # FI outcome used by active cancer trajectory analyses is observed.
  filter(participated == 1, !is.na(fi_score_nocancer)) %>%
  
  mutate(
    # 1. Readable Date of Return
    date_of_return = make_date(year = 1900 + floor((worked_rtmnyr - 1) / 12), 
                               month = ((worked_rtmnyr - 1) %% 12) + 1, 
                               day = 15),
    
    # 2. Exact Age at this specific cycle
    age_at_cycle = (worked_rtmnyr - dbmy09) / 12,
    
    # 3. "Time to Death" (Retrospective timeline)
    # Negative means the survey was taken BEFORE death.
    years_to_death = (worked_rtmnyr - dtdth) / 12,
    
    # 4. "Time to Cancer" (Retrospective/Prospective timeline)
    # Negative = survey was BEFORE cancer diagnosis. Positive = survey was AFTER.
    years_to_cancer = (worked_rtmnyr - cancer_dateca) / 12
  ) %>%
  
  # Group by ID to calculate "Time Since Baseline" 
  arrange(id, worked_rtmnyr) %>%
  group_by(id) %>%
  mutate(
    # How many years have passed since this participant's first survey?
    years_since_baseline = (worked_rtmnyr - first(worked_rtmnyr)) / 12
  ) %>%
  ungroup()

# ==============================================================================
# 3. CATEGORIZE VARIABLES
# ==============================================================================
fi_trajectory <- fi_trajectory %>%
  mutate(
    # Detailed Age at Death Category (Only populates for decedents)
    age_at_death_cat5 = case_when(
      age_at_death <= 70 ~ "<=70",
      age_at_death > 70 & age_at_death < 80 ~ "70-79",
      age_at_death >= 80 & age_at_death < 90 ~ "80-89",
      age_at_death >= 90 & age_at_death < 100 ~ "90-99",
      age_at_death >= 100 ~ "100+",
      TRUE ~ NA_character_
    ),
    
    # Binary Age at Death Category (<=80 vs >80)
    age_at_death_cat2 = case_when(
      age_at_death <= 80 ~ "<=80",
      age_at_death > 80 ~ ">80",
      TRUE ~ NA_character_
    )
  ) %>%
  # Convert them to ordered factors so they plot correctly in ggplot
  mutate(
    age_at_death_cat5 = factor(age_at_death_cat5, 
                               levels = c("<=70", "70-79", "80-89", "90-99", "100+")),
    age_at_death_cat2 = factor(age_at_death_cat2, 
                               levels = c("<=80", ">80"))
  )

# ==============================================================================
# 4. VERIFICATION & SAVE
# ==============================================================================
print("Sanity Check for participants who got cancer:")
fi_trajectory %>%
  filter(!is.na(cancer_date)) %>% # Look at someone with a cancer diagnosis
  filter(id == first(id)) %>%     # Just look at their specific trajectory
  select(id, cycle, date_of_return, cancer_date, age_at_cycle, age_at_cancer, years_to_cancer) %>%
  print()

# Save the trajectory-ready dataset as a distinct derived file. Do not overwrite
# the broader cancer analytic dataset used as the input to this step.
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

save_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer_timevars.rds")
saveRDS(fi_trajectory, save_path)

cat("\nDataset with timelines calculated and saved to:", save_path, "\n")
