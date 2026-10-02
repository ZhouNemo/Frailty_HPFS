# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  2.6_create_riskset_individual_types.R
# Author:  Nemo Zhou
# Date started:      2026-09-29
# Date last updated: 2026-09-29
#
# Purpose:
#   Independently match six first-cancer site cohorts to cancer-free risk-set
#   controls using the active 2.0 incidence-density design. The canonical
#   earliest-diagnosis ICD list, not lifetime site flags, defines case type.
#   Same-date different-site cancers enter each applicable cohort. One script
#   builds all six cohorts by default; --site=<slug> selects one. It never runs
#   GLME models or changes the existing overall/grouped/prostate artifacts.
#
# Inputs:  Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
# Outputs: Data/riskset_matched_site_<slug>_{long,assignments}.rds;
#          Results/cancer/data/matching_diagnostics/riskset_matched_site_*
#          and 2.6_individual_type_run_status.csv.
# Run after 7.4_cancer_subtypes.R; Nemo runs this long matching script.
# =============================================================================

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
source(file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R"))

site_icd_codes <- list(
  lung = 162L,
  colorectal = c(153L, 154L),
  prostate = 185L,
  bladder = 188L,
  pancreas = 157L,
  kidney = 189L
)
site_labels <- c(
  lung = "Lung/Trachea/Bronchus Cancer Cohort",
  colorectal = "Colorectal/Anal Cancer Cohort",
  prostate = "Prostate Cancer Cohort",
  bladder = "Bladder Cancer Cohort",
  pancreas = "Pancreatic Cancer Cohort",
  kidney = "Kidney/Other Urinary Cancer Cohort"
)

parse_index_icds <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) return(integer())
  tokens <- trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
  if (any(!grepl("^[0-9]+$", tokens)))
    stop("Malformed cancer_index_icds: ", as.character(x), call. = FALSE)
  codes <- suppressWarnings(as.integer(tokens))
  if (anyNA(codes)) stop("Out-of-range cancer_index_icds: ", as.character(x), call. = FALSE)
  unique(codes)
}

site_membership_table <- function(person) {
  required <- c("id", "cancer_index_dateca", "cancer_index_icds")
  if (!all(required %in% names(person))) stop("Missing canonical index fields.", call. = FALSE)
  ids <- as.character(person$id)
  if (anyNA(ids) || anyDuplicated(ids)) stop("Expected one nonmissing row per participant.", call. = FALSE)
  dates <- suppressWarnings(as.numeric(as.character(person$cancer_index_dateca)))
  codes <- lapply(person$cancer_index_icds, parse_index_icds)
  out <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (site in names(site_icd_codes)) {
    wanted <- site_icd_codes[[site]]
    out[[site]] <- is.finite(dates) & vapply(codes, function(x) any(x %in% wanted), logical(1))
  }
  out
}

requested_sites <- function(args) {
  if (!length(args)) return(names(site_icd_codes))
  if (length(args) != 1L || !grepl("^--site=", args[[1]]))
    stop("Usage: Rscript Code/2_data_analysis/2.6_create_riskset_individual_types.R [--site=<slug>]", call. = FALSE)
  site <- sub("^--site=", "", args[[1]])
  if (!site %in% names(site_icd_codes))
    stop("Unknown site: ", site, "; choose ", paste(names(site_icd_codes), collapse = ", "), call. = FALSE)
  site
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  sites <- requested_sites(args)
  data_dir <- file.path(project_dir, "Data")
  input_path <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
  diagnostics_dir <- file.path(project_dir, "Results", "cancer", "data", "matching_diagnostics")
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(diagnostics_dir, "2.6_individual_type_run_status.csv")
  builder_path <- file.path(project_dir, "Code", "2_data_analysis", "2.6_create_riskset_individual_types.R")
  engine_path <- file.path(project_dir, "Code", "2_data_analysis", "2.0_riskset_matching_functions.R")
  input_md5 <- unname(tools::md5sum(input_path))
  builder_md5 <- unname(tools::md5sum(builder_path))
  engine_md5 <- unname(tools::md5sum(engine_path))

  fi <- readRDS(input_path)
  needed <- c("id", "cancer_index_dateca", "cancer_index_icds")
  if (!all(needed %in% names(fi))) stop("Run 7.4_cancer_subtypes.R before 2.6.", call. = FALSE)
  person <- dplyr::distinct(dplyr::select(fi, dplyr::all_of(needed)))
  if (anyDuplicated(as.character(person$id)))
    stop("Canonical index fields vary within participant.", call. = FALSE)
  membership <- site_membership_table(person)
  rm(fi, person)

  write_status <- function(site, status, output_md5 = NA_character_, message = "") {
    previous <- if (file.exists(status_path)) read.csv(status_path, stringsAsFactors = FALSE) else data.frame()
    if (nrow(previous)) previous <- previous[previous$site != site, , drop = FALSE]
    row <- data.frame(
      site = site, status = status, input_md5 = input_md5,
      builder_md5 = builder_md5, matching_engine_md5 = engine_md5,
      output_md5 = output_md5, message = message,
      updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      stringsAsFactors = FALSE
    )
    write.csv(rbind(previous, row), status_path, row.names = FALSE)
  }
  for (site in sites) write_status(site, "pending")
  failures <- character()
  for (site in sites) {
    write_status(site, "running")
    tryCatch({
      case_ids <- membership$id[membership[[site]]]
      label <- unname(site_labels[[site]])
      matched <- build_riskset_matched_long(
        input_path = input_path,
        classification_vars = character(0),
        classify_fn = function(pl) ifelse(as.character(pl$id) %in% case_ids, label, NA_character_),
        cohort_levels = label,
        target_cycles = c("88", "92", "96", "00", "04", "08", "12", "16", "20"),
        match_ratio = 5L, age_caliper = 2, seed = 20260703L
      )
      matched$run_metadata$site <- site
      matched$run_metadata$site_icds <- as.integer(site_icd_codes[[site]])
      matched$run_metadata$site_label <- label
      matched$run_metadata$case_definition <- "canonical earliest dated cancer ICDs only"
      matched$run_metadata$tie_policy <- "include in each applicable site cohort"
      matched$run_metadata$builder_md5 <- builder_md5
      matched$run_metadata$matching_engine_md5 <- engine_md5
      output_path <- file.path(data_dir, paste0("riskset_matched_site_", site, "_long.rds"))
      save_riskset_match(matched, output_path, paste(label, "risk-set cohort"))
      write_status(site, "matched", unname(tools::md5sum(output_path)))
    }, error = function(e) {
      failures <<- c(failures, site)
      write_status(site, "failed", message = conditionMessage(e))
      message("Site ", site, " failed: ", conditionMessage(e))
    })
  }
  if (length(failures)) stop("Matching failed for: ", paste(failures, collapse = ", "),
                             ". See ", status_path, call. = FALSE)
  invisible(status_path)
}

if (sys.nframe() == 0L) main()
