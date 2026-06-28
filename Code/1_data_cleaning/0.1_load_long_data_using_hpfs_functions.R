# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 0.1_load_long_data_using_hpfs_functions.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Loads core HPFS longitudinal data using the internal hpfs package, including baseline, follow-up, deaths, race, anthropometrics, smoking, physical activity, and other questionnaire-derived datasets needed for frailty-index construction.
# ==============================================================================

#install.packages("labelled")
#nstall.packages("/n/hpnh/HPFS/R/hpfs_0.2.0.tar.gz", repos = NULL)

library(dplyr)
library(hpfs)

# --- Main Cohort Data ---
print("Loading Baseline...")
baseline_df <- load_hpfs_baseline()

print("Loading Follow-up...")
followup_df <- load_hpfs_followup()

print("Loading Deaths...")
deaths_df <- load_hpfs_deaths()

print("Loading Race...")
race_df <- load_hpfs_race()


# --- Lifestyle & Measurements ---
print("Loading Anthropometrics...")
anthropometrics_df <- load_hpfs_anthropometrics()

print("Loading Smoking...")
smoking_df <- load_hpfs_smoking()

print("Loading Physical Activity...")
physact_df <- load_hpfs_physact()

print("Loading Nutrients...")
nutrients_df <- load_hpfs_nutrients()


# --- Clinical & Biospecimen Data ---
print("Loading Blood Data...")
blood_df <- load_hpfs_blood()

print("Loading Blood Pressure...")
bloodpressure_df <- load_hpfs_bloodpressure()

print("Loading Cheek Cell Data...")
cheek_df <- load_hpfs_cheek()

print("Loading Cardio...")
cardio_df <- load_hpfs_cardio()


# --- Prostate Specific Data ---
print("Loading Family History (Prostate)...")
famhx_prostate_df <- load_hpfs_famhx_prostate()

print("Loading Prostate Casefile...")
prostate_casefile_df <- load_hpfs_prostate_casefile()

print("Loading PSA Screening...")
psascreening_df <- load_hpfs_psascreening()







###############################################
# --- look at the age at baseline and 2010 ---

baseline_df <- baseline_df %>%
  mutate(age_years = age_enrollment / 12)

# 2. Summary Statistics
print("Summary of Age at Enrollment (Years):")
summary(baseline_df$age_years)

# 3. Create Histogram
# Using a binwidth of 1 year to see the distribution clearly
ggplot(baseline_df, aes(x = age_years)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Age at Enrollment (1986)",
    subtitle = "Health Professionals Follow-Up Study",
    x = "Age (Years)",
    y = "Number of Participants"
  )



# 1. Identify participants who responded in Cycle 10 (2010)
#    We filter the follow-up data specifically for this cycle.
participants_2010 <- followup_df %>% 
  filter(cycle == "10") %>% 
  select(id)

# 2. Merge with Baseline Age and Calculate 2010 Age
#    Formula: (Age at Enrollment in Months / 12) + 24 Years
df_age_2010 <- participants_2010 %>% 
  inner_join(baseline_df, by = "id") %>% 
  mutate(
    age_2010 = (age_enrollment / 12) + 24
  )

# 3. Summary Statistics
print("Summary of Age in 2010 (Derived from Baseline):")
summary(df_age_2010$age_2010)

# 4. Histogram
ggplot(df_age_2010, aes(x = age_2010)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Age Distribution of Participants in 2010",
    subtitle = "Calculated as: Baseline Age + 24 Years",
    x = "Age (Years)",
    y = "Count"
  )


# --- Optional: PHS Data ---
# (Included in the package, uncomment if needed)
# print("Loading PHS Prostate Casefile...")
# phs_prostate_casefile_df <- load_phs_prostate_casefile()

print("All datasets loaded successfully.")