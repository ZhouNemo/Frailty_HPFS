# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.0.2_GLME_descriptive_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-09-21
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
#
# Purpose: Build model-free, observed frailty-index trajectory summaries for
# the active 4.x risk-set-matched cohorts. The script uses the exact M1
# analytic population and the same diagnosis-relative 2-year support bins as
# 4.0, but does not fit a model or calculate inferential uncertainty.
#
# Input:  Data/riskset_matched_{analysis,smoking,obesity,overall}_long.rds
# Output: Results/cancer/data/4.0.2_descriptive_*.csv and metadata RDS.
# Render: Code/2_data_analysis/4.0.2_GLME_descriptive_trajectories.Rmd
# Usage:  Rscript Code/2_data_analysis/4.0.2_GLME_descriptive_trajectories.R
# =============================================================================

source("Code/2_data_analysis/4.0_GLME_spline_functions.R")

.descriptive_cohort_inputs <- data.frame(
  prefix = c("4.1", "4.2", "4.3", "4.5"),
  matched_file = c(
    "riskset_matched_analysis_long.rds",
    "riskset_matched_smoking_long.rds",
    "riskset_matched_obesity_long.rds",
    "riskset_matched_overall_long.rds"
  ),
  stringsAsFactors = FALSE
)

.descriptive_cohort_order <- c(
  "Low/Moderate Burden Cohort",
  "High Burden Cohort",
  "Smoking-Related Cancer Cohort",
  "Obesity-Related Cancer Cohort",
  "All Cancer Cohort"
)

.descriptive_bin_midpoint <- function(bin) {
  bin <- as.character(bin)
  midpoint <- rep(NA_real_, length(bin))
  inside <- !(bin %in% c("<= -20 years", "> +20"))
  bounds <- strsplit(bin[inside], " to ", fixed = TRUE)
  midpoint[inside] <- vapply(bounds, function(x) mean(as.numeric(gsub("\\+", "", x))), numeric(1))
  midpoint
}

.prepare_descriptive_m1_data <- function(data) {
  required <- c(
    "Cohort", "Group", "Age_Centered", "id", "match_set", "role", "cycle",
    "fi_score_nocancer", "post_own_cancer", .primary_covars
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Descriptive input is missing: ", paste(missing, collapse = ", "))
  if (!"trajectory_id" %in% names(data)) {
    data$trajectory_id <- paste(data$Cohort, data$match_set, data$id, data$role, sep = "__")
  }
  if (anyNA(data[c("id", "Group", "trajectory_id", "cycle", "post_own_cancer")])) {
    stop("Descriptive input has missing identifiers, Group, cycle, or own-cancer censoring flag")
  }
  if (!all(as.character(data$Group) %in% c("Control", "Cancer Case"))) stop("Unexpected Group values")
  if (any(!is.na(data$fi_score_nocancer) &
          (!is.finite(data$fi_score_nocancer) | data$fi_score_nocancer < 0 | data$fi_score_nocancer > 1))) {
    stop("FI must be finite and in [0, 1]")
  }
  if (any(!is.na(data$Age_Centered) & !is.finite(data$Age_Centered))) stop("Relative time must be finite")

  out <- data %>%
    mutate(
      Cohort = as.character(Cohort),
      Group = factor(as.character(Group), levels = c("Control", "Cancer Case")),
      id = as.character(id),
      trajectory_id = as.character(trajectory_id),
      cycle = as.character(cycle)
    ) %>%
    filter(
      !post_own_cancer,
      !is.na(fi_score_nocancer),
      !is.na(Age_Centered),
      if_all(all_of(.primary_covars), ~ !is.na(.x))
    ) %>%
    add_relative_time_bin()

  if (!nrow(out) || n_distinct(out$Group) != 2L) stop("M1 descriptive population lacks two-arm data")
  if (anyDuplicated(paste(out$trajectory_id, out$cycle, sep = "__"))) {
    stop("Duplicated trajectory_id x cycle rows remain after descriptive M1 filtering")
  }
  out
}

.summarize_descriptive_cohort <- function(data, prefix, provenance,
                                           min_case_bin = 50L, min_ctrl_bin = 250L,
                                           window_yrs = 20) {
  observed <- data %>%
    group_by(Cohort, Group, rel_time_bin, .drop = FALSE) %>%
    summarize(
      n_observations = n(),
      n_assignments = n_distinct(trajectory_id),
      n_participants = n_distinct(id),
      mean_fi = mean(fi_score_nocancer),
      .groups = "drop"
    ) %>%
    mutate(prefix = prefix, relative_time_mid = .descriptive_bin_midpoint(rel_time_bin))

  support <- data %>%
    group_by(Cohort, rel_time_bin, .drop = FALSE) %>%
    summarize(
      n_case_observations = sum(Group == "Cancer Case"),
      n_control_observations = sum(Group == "Control"),
      n_case_assignments = n_distinct(trajectory_id[Group == "Cancer Case"]),
      n_control_assignments = n_distinct(trajectory_id[Group == "Control"]),
      n_case_participants = n_distinct(id[Group == "Cancer Case"]),
      n_control_participants = n_distinct(id[Group == "Control"]),
      .groups = "drop"
    ) %>%
    mutate(
      prefix = prefix,
      relative_time_mid = .descriptive_bin_midpoint(rel_time_bin),
      support_ok = n_case_participants >= min_case_bin & n_control_participants >= min_ctrl_bin
    )

  window <- continuous_support_window(support %>% select(rel_time_bin, support_ok), window_yrs)
  support <- support %>% mutate(
    supported_window_min = window[[1]],
    supported_window_max = window[[2]],
    display_supported = is.finite(relative_time_mid) & is.finite(window[[1]]) &
      relative_time_mid >= window[[1]] & relative_time_mid <= window[[2]] & support_ok
  )

  observed <- observed %>%
    left_join(support %>% select(Cohort, rel_time_bin, display_supported), by = c("Cohort", "rel_time_bin"))
  difference <- observed %>%
    select(Cohort, prefix, rel_time_bin, relative_time_mid, Group, mean_fi, display_supported) %>%
    tidyr::pivot_wider(names_from = Group, values_from = mean_fi) %>%
    transmute(
      Cohort, prefix, rel_time_bin, relative_time_mid, display_supported,
      mean_fi_control = Control,
      mean_fi_cancer_case = `Cancer Case`,
      observed_difference = `Cancer Case` - Control
    )
  eligibility <- data.frame(
    prefix = prefix,
    Cohort = unique(data$Cohort),
    n_observations = nrow(data),
    n_assignments = dplyr::n_distinct(data$trajectory_id),
    n_participants = dplyr::n_distinct(data$id),
    input_md5 = provenance$input_md5,
    stringsAsFactors = FALSE
  )
  list(observed = observed, difference = difference, support = support, eligibility = eligibility,
       provenance = provenance)
}

run_descriptive_trajectories <- function(project_dir = getwd(), window_yrs = 20L,
                                         min_case_bin = 50L, min_ctrl_bin = 250L) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  data_dir <- file.path(project_dir, "Data")
  results_dir <- file.path(project_dir, "Results", "cancer", "data")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(identical(window_yrs, 20L), identical(min_case_bin, 50L), identical(min_ctrl_bin, 250L))

  population_flow <- list()
  pieces <- lapply(seq_len(nrow(.descriptive_cohort_inputs)), function(i) {
    info <- .descriptive_cohort_inputs[i, ]
    path <- file.path(data_dir, info$matched_file)
    provenance <- validate_matching_provenance(path)
    raw <- readRDS(path)
    prepared <- .prepare_descriptive_m1_data(raw)
    ledger <- read_matching_assignments(path)
    population_flow[[i]] <<- ledger %>% group_by(Cohort, Group) %>% summarize(
      n_matched_assignments = n(),
      n_fi_contributing_assignments = sum(trajectory_id %in% raw$trajectory_id[!raw$post_own_cancer]),
      n_m1_complete_assignments = sum(trajectory_id %in% prepared$trajectory_id),
      n_without_usable_fi = n_matched_assignments - n_fi_contributing_assignments,
      n_excluded_m1_completeness = n_fi_contributing_assignments - n_m1_complete_assignments,
      .groups = "drop") %>% mutate(prefix = info$prefix,
        input_md5 = provenance$input_md5, assignment_md5 = provenance$assignment_md5)
    lapply(split(prepared, prepared$Cohort), function(one_cohort) {
      .summarize_descriptive_cohort(
        one_cohort, info$prefix, provenance,
        min_case_bin = min_case_bin, min_ctrl_bin = min_ctrl_bin, window_yrs = window_yrs
      )
    })
  })
  pieces <- unlist(pieces, recursive = FALSE)

  write.csv(bind_rows(population_flow), file.path(results_dir, "4.0.2_descriptive_assignment_flow.csv"), row.names = FALSE)
  observed <- bind_rows(lapply(pieces, `[[`, "observed"))
  difference <- bind_rows(lapply(pieces, `[[`, "difference"))
  support <- bind_rows(lapply(pieces, `[[`, "support"))
  eligibility <- bind_rows(lapply(pieces, `[[`, "eligibility"))
  order_cohort <- function(x) transform(x, Cohort = factor(Cohort, levels = .descriptive_cohort_order)) %>%
    arrange(Cohort, rel_time_bin)
  observed <- order_cohort(observed)
  difference <- order_cohort(difference)
  support <- order_cohort(support)
  eligibility <- eligibility %>% mutate(Cohort = factor(Cohort, levels = .descriptive_cohort_order)) %>% arrange(Cohort)

  write.csv(observed, file.path(results_dir, "4.0.2_descriptive_observed_means.csv"), row.names = FALSE)
  write.csv(difference, file.path(results_dir, "4.0.2_descriptive_group_difference.csv"), row.names = FALSE)
  write.csv(support, file.path(results_dir, "4.0.2_descriptive_support.csv"), row.names = FALSE)
  write.csv(eligibility, file.path(results_dir, "4.0.2_descriptive_m1_eligibility.csv"), row.names = FALSE)
  metadata <- list(
    schema_version = 1L,
    analysis = "model-free observed FI trajectories",
    outcome = "fi_score_nocancer",
    filters = c("M1 primary covariate complete case", "own-cancer censoring", "observed FI and relative time"),
    time_scale = "diagnosis-relative Age_Centered",
    bin_definition = "4.0 relative-time bins: 2 years, right-closed",
    window_years = window_yrs,
    support_thresholds = c(case_participants = min_case_bin, control_participants = min_ctrl_bin),
    uncertainty_intervals = FALSE,
    assignment_weighted_means = TRUE,
    input_provenance = lapply(pieces, `[[`, "provenance"),
    script_md5 = unname(tools::md5sum(file.path(project_dir, "Code", "2_data_analysis", "4.0.2_GLME_descriptive_trajectories.R"))),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  saveRDS(metadata, file.path(results_dir, "4.0.2_descriptive_run_metadata.rds"))
  message("Saved descriptive trajectory data to: ", results_dir)
  invisible(list(observed = observed, difference = difference, support = support,
                 eligibility = eligibility, metadata = metadata))
}

render_descriptive_trajectories <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  report_rmd <- file.path(project_dir, "Code", "2_data_analysis", "4.0.2_GLME_descriptive_trajectories.Rmd")
  visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")
  output_html <- file.path(visuals_dir, "4.0.2_GLME_descriptive_trajectories.html")
  if (!file.exists(report_rmd)) stop("Descriptive report is missing: ", report_rmd)
  dir.create(visuals_dir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    rmarkdown::render(report_rmd, output_file = basename(output_html), output_dir = visuals_dir,
                      quiet = TRUE, envir = new.env(parent = globalenv()))
  } else {
    if (!requireNamespace("knitr", quietly = TRUE) || !requireNamespace("markdown", quietly = TRUE)) {
      stop("Rendering without Pandoc requires both knitr and markdown")
    }
    knitted <- tempfile(fileext = ".md")
    knitr::knit(report_rmd, output = knitted, quiet = TRUE, envir = new.env(parent = globalenv()))
    markdown::markdownToHTML(
      knitted, output = output_html,
      options = c("+embed_resources", "+toc", "+table", "+auto_identifiers"),
      title = "Observed Frailty Trajectories in Risk-Set Matched Cancer Cohorts"
    )
  }
  if (!file.exists(output_html) || file.info(output_html)$size <= 0) {
    stop("Descriptive report renderer did not create a nonempty HTML file")
  }
  message("Saved descriptive HTML report to: ", output_html)
  invisible(output_html)
}

if (sys.nframe() == 0L) run_descriptive_trajectories()
