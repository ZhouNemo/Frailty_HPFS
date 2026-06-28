# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 7.3_cancer records.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Reads fixed-width HPFS cancer record files, assigns cancer-related variable names and positions, and prepares cancer diagnosis and cancer-type information for linkage to frailty data.
# ==============================================================================

library(readr)
library(dplyr)


# ---------------------------------------------------------
# READ IN CANCER FILE
# ---------------------------------------------------------
# 1. Define start positions for 2020
starts_20 <- c(
  1, 8, 12, 16, 18, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 
  36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 52, 56, 58, 62, 66, 70, 
  74, 77, 80, 83, 86, 90, 94, 98, 102, 105, 108, 111, 114, 118, 122, 126
)

# 2. Define end positions for 2020
ends_20 <- c(
  6, 11, 15, 16, 18, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 
  36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 54, 56, 60, 64, 68, 72, 
  75, 78, 81, 84, 88, 92, 96, 100, 103, 106, 109, 112, 116, 120, 124, 128
)

# 3. Define column names for 2020
cols_20 <- c(
  "hpfsid", "dateca", "dod", "cancer", "mel", "lal", "other", "colorect", 
  "lung", "carcin", "newag", "cancerp", "carcinp", "prosnoa1", "acan", 
  "bladder", "pancreas", "kidney", "stomach", "liver", "brain", 
  "sarcoma", "oralc", "status", "esoph", "leuk", "lymp_nh", "lymp_h", 
  "myeloma", "pharyn", "oroph", "conf", "cancerdth", "icdx", "dupnum", 
  "firsticd", "secondicd", "thirdicd", "fourthicd", "firstconf", "secondconf", 
  "thirdconf", "fourthconf", "firstyodx", "secondyodx", "thirdyodx", "fourthyodx", 
  "firstmm", "secondmm", "thirdmm", "fourthmm", "firstyy", "secondyy", 
  "thirdyy", "fourthyy"
)

# 4. Read the 2020 file
cancer2020 <- read_fwf(
  file = "/n/hpnh/HPFS/Endpoint/Cancer/AllCancer/2020/Cancer2020_March2024_withduplicate.dat",
  col_positions = fwf_positions(start = starts_20, end = ends_20, col_names = cols_20)
)

# View the structure of your loaded dataset
str(cancer2020)


# ---------------------------------------------------------
# MERGE TO THE MAIN DATASET 
# ---------------------------------------------------------


# 1. Load your main imputed dataset
FI_longitudinal_1986_2020_IMPUTED <- readRDS("~/Frailty/FI_longitudinal_1986_2020_IMPUTED.rds")

# Assuming 'cancer2020' is already loaded in your environment from the previous step

# 2. Process and merge the cancer dataset
fi_long_merged <- FI_longitudinal_1986_2020_IMPUTED %>%
  left_join(
    cancer2020 %>%
      # Rename hpfsid to id to match the main dataset
      rename(id = hpfsid) %>%
      # Ensure id is a character (matching your previous merge logic)
      mutate(id = as.character(id)) %>%
      # Add "cancer_" prefix to all columns EXCEPT 'id'
      rename_with(~ paste0("cancer_", .), -id),
    by = "id" # Merging only by id since cancer is a person-level endpoint
  )

# ---------------------------------------------------------
# QUALITY CHECKS 
# ---------------------------------------------------------

# Check 1: Check if some of the new columns exist
c("cancer_dateca", "cancer_cancer", "cancer_prostate") %in% names(fi_long_merged)

# Check 2: Check the summary of a few key cancer variables
summary(fi_long_merged[c("cancer_dateca", "cancer_cancer", "cancer_cancerdth")])

# Check 3: Spot check a few rows (showing id and some newly prefixed cancer columns)
fi_long_merged %>%
  select(id, cycle, starts_with("cancer_")) %>%
  filter(!is.na(cancer_dateca)) %>%
  head(15)

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
