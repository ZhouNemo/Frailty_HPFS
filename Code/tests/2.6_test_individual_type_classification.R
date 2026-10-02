# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  2.6_test_individual_type_classification.R
# Author:  Nemo Zhou
# Date started:      2026-09-29
# Date last updated: 2026-09-29
# Purpose: Fast classification and CLI contract checks for the six site cohorts.
#          This file does not read project data, match risk sets, or fit models.
# =============================================================================

source("Code/2_data_analysis/2.6_create_riskset_individual_types.R")

person <- data.frame(
  id = c("lung", "colon", "rectal_anal", "prostate", "bladder", "pancreas",
         "kidney", "tie", "later_site", "none", "boundary", "undated"),
  cancer_index_dateca = c(rep(1100, 9), NA, 1100, NA),
  cancer_index_icds = c("162", "153", "154", "185", "188", "157", "189",
                        "153;157", "185", NA, "1880", "162"),
  stringsAsFactors = FALSE
)
membership <- site_membership_table(person)
has <- function(id, site) membership[[site]][match(id, membership$id)]
stopifnot(has("lung", "lung"), !has("lung", "bladder"))
stopifnot(has("colon", "colorectal"), has("rectal_anal", "colorectal"))
stopifnot(has("prostate", "prostate"), has("bladder", "bladder"))
stopifnot(has("pancreas", "pancreas"), has("kidney", "kidney"))
stopifnot(has("tie", "colorectal"), has("tie", "pancreas"))
stopifnot(!has("boundary", "bladder"), !any(unlist(membership[membership$id == "none", -1])))
stopifnot(!has("undated", "lung"))
# The source has only the first-date ICD list: a later lung event cannot turn
# this prostate index case into a lung case.
stopifnot(has("later_site", "prostate"), !has("later_site", "lung"))
stopifnot(identical(parse_index_icds("153; 154;153"), c(153L, 154L)))
stopifnot(identical(requested_sites(character()), names(site_icd_codes)))
stopifnot(identical(requested_sites("--site=kidney"), "kidney"))
stopifnot(inherits(try(requested_sites("--site=unknown"), silent = TRUE), "try-error"))
stopifnot(inherits(try(parse_index_icds("153;abc"), silent = TRUE), "try-error"))
cat("Individual cancer-site classification contracts passed.\n")
