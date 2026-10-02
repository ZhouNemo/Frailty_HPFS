# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-29 (M2 alignment; explicit unresolved date-policy block)
# Purpose: Pure preparation, contract checks, supported nlme/JMbayes2 fitting,
# posterior extraction, and fixed-effect comparisons for the overall mortality
# sensitivity. Sourcing this file never fits a model or writes an output.
# One matched assignment per participant is selected before outcome exclusions.
# Pre-index FI is history; survival starts at index, conditional on entry alive.
# Canonical M2 spline and attained-age random designs update at hazard times;
# dietary covariates use saved M2 transformations exactly once.
# All profiles and direct fitting helpers currently block on unresolved dates.
# =============================================================================

source("/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS/Code/2_data_analysis/2.0_matching_provenance.R")

jm_alignment_model <- "M2_full_spline"
jm_primary_covars <- c("index_age_z", "base_race", "base_marital", "base_living")
jm_diet_covars <- c("base_calor", "base_sat", "base_diet_chol", "base_alco")
jm_covars <- c(jm_primary_covars, "base_pckgr", jm_diet_covars)
jm_factors <- c("base_race", "base_marital", "base_living", "base_pckgr")
jm_hazard_control <- list(basis="bs", Bsplines_degree=2L,
  base_hazard_segments=9L, timescale_base_hazard="identity", diff=2L)
# Deliberately no environment/CLI override: a later reviewed date-policy revision
# must resolve this gate. Synthetic tests replace this helper only in test scope.
jm_date_policy <- function() list(status="unresolved", date_precision="month",
  decision="Keep fitting blocked; preserve same-month cases and source dates",
  exact_day_availability="Not verified on the inaccessible HPFS source server")
jm_require_resolved_dates <- function() {
  if (jm_date_policy()$status != "resolved") stop(structure(list(
    message="Fitting blocked: month-level diagnosis/death and administrative-end date policy remains unresolved. Review saved date audits; no dates or cases were changed.",
    call=NULL), class=c("jm_date_policy_error", "error", "condition")))
  invisible(TRUE)
}
jm_admin_date <- (2020 - 1900) * 12 + 6

jm_assert <- function(ok, message, audit = NULL, audit_name = "validation_failures") {
  if (!isTRUE(ok)) stop(structure(list(message = message, call = NULL,
    audit = audit, audit_name = audit_name), class = c("jm_validation_error", "error", "condition")))
}
jm_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
jm_profile <- function(x = Sys.getenv("JM_PROFILE", "prepare")) {
  jm_assert(length(x) == 1L && x %in% c("prepare", "pilot", "full"),
            "JM_PROFILE must be exactly prepare, pilot, or full.")
  x
}
jm_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
jm_save <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp)
  jm_assert(file.rename(tmp, path), paste("Could not commit checkpoint:", path))
}
jm_new_run <- function(root, profile) {
  profile <- jm_profile(profile)
  parent <- file.path(root, profile)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  run <- tempfile(paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "_"), parent)
  jm_assert(dir.create(run), "Could not create an isolated JM run directory.")
  normalizePath(run)
}
jm_need <- function(data, columns, label) {
  miss <- setdiff(columns, names(data))
  jm_assert(!length(miss), paste(label, "is missing:", paste(miss, collapse = ", ")))
}

# Strict per-participant/assignment consistency, allowing repeated missing dates.
jm_consistent <- function(data, keys, columns, label) {
  bad <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(columns),
                                  ~ dplyr::n_distinct(.x, na.rm = TRUE)), .groups = "drop")
  invalid <- rowSums(as.matrix(bad[columns]) > 1L) > 0
  jm_assert(!any(invalid), paste("Conflicting", label, "within identifier."), bad[invalid, ])
}

jm_validate_inputs <- function(project_dir) {
  data_dir <- file.path(project_dir, "Data")
  result_dir <- file.path(project_dir, "Results", "cancer", "data")
  stem <- "riskset_matched_overall_long"
  paths <- c(matched = file.path(data_dir, paste0(stem, ".rds")),
             assignments = file.path(data_dir, "riskset_matched_overall_assignments.rds"),
             panel = file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds"),
             gate = file.path(result_dir, "matching_diagnostics", paste0(stem, "_gate_g4.csv")),
             matching = file.path(result_dir, "matching_diagnostics", paste0(stem, "_run_metadata.rds")),
             age_scaling = file.path(result_dir, "matching_diagnostics", paste0(stem, "_scaling_metadata.rds")),
             spline = file.path(result_dir, "glme_spline_metadata.rds"))
  jm_assert(all(file.exists(paths)), paste("Missing prerequisite:", paste(paths[!file.exists(paths)], collapse = ", ")))
  hashes <- stats::setNames(as.character(tools::md5sum(paths)), names(paths))
  gate <- utils::read.csv(paths[["gate"]])
  jm_need(gate, "gate_pass", "Gate G4")
  jm_assert(nrow(gate) > 0 && all(as.logical(gate$gate_pass) %in% TRUE), "Matching Gate G4 failed.")
  provenance <- readRDS(paths[["matching"]])
  assert_matching_metadata(provenance, paths[["matched"]])
  validate_assignment_provenance(provenance, paths[["matched"]])
  jm_assert(identical(unname(hashes[["matched"]]), unname(provenance$output_md5)),
            "Matched input hash differs from Gate G4 provenance. Nemo must rebuild 2.5.")
  jm_assert(identical(unname(hashes[["panel"]]), unname(provenance$input_md5)),
            "Canonical panel differs from matching source hash. Nemo must reconcile/rebuild 2.5.",
            data.frame(source = paths[["panel"]], expected_md5 = provenance$input_md5,
                       observed_md5 = hashes[["panel"]]), "source_hash_mismatch")
  matched <- readRDS(paths[["matched"]])
  assignments <- readRDS(paths[["assignments"]])
  assert_assignment_rows(assignments, "JM assignment ledger")
  panel <- readRDS(paths[["panel"]])
  metadata <- readRDS(paths[["spline"]])
  age_scaling <- readRDS(paths[["age_scaling"]])
  jm_need(matched, c("id", "Cohort", "Group", "role", "match_set", "cycle", "index_date",
                    "index_age", "age_at_cycle", "worked_rtmnyr", "Age_Centered", "fi_score_nocancer",
                    "cancer_dateca", "post_own_cancer", jm_covars), "Matched data")
  jm_need(panel, c("id", "dtdth"), "Panel")
  expected_cycles <- c("88", "92", "96", "00", "04", "08", "12", "16", "20")
  jm_assert(setequal(as.character(provenance$target_cycles), expected_cycles),
            "Matching provenance does not specify the active four-year analytic cycles.")
  jm_assert(all(as.character(matched$cycle) %in% expected_cycles), "Unexpected analytic cycle in matched data.")
  jm_need(age_scaling, c("index_age_mean", "index_age_sd"), "Index-age scaling")
  jm_assert(nrow(age_scaling) == 1L && is.finite(age_scaling$index_age_sd) && age_scaling$index_age_sd > 0,
            "Invalid index-age scaling artifact.")
  error <- abs(jm_numeric(matched$index_age_z) -
                 (jm_numeric(matched$index_age) - age_scaling$index_age_mean) / age_scaling$index_age_sd)
  jm_assert(all(is.finite(error)) && max(error) < 1e-8, "Upstream index-age scaling is inconsistent.")
  jm_assert(identical(metadata$input_md5, unname(hashes[["matched"]])) &&
              identical(metadata$assignment_md5, unname(hashes[["assignments"]])),
            "M2 alignment metadata has stale ledger/outcome hashes; Nemo must rerun Code/2_data_analysis/4.5_GLME_overall.R.")
  jm_validate_metadata(metadata, matched)
  list(matched = matched, assignments = assignments, panel = panel, metadata = metadata,
       provenance = list(paths = paths, hashes = hashes, matching = provenance,
                         age_scaling = age_scaling,
                         spline_validation = "M1 reference-frame knots and M2 complete-case scaling verified against current matched input"))
}

jm_validate_metadata <- function(m, matched = NULL) {
  required <- c("schema_version", "cohort", "outcome", "primary_model_id", "spline_df",
    "knots", "boundary_knots", "primary_covariates", "model_scaling", "model_specification",
    "spline_reference_model_id", "random_clock", "random_time_center", "random_time_scale")
  jm_assert(all(required %in% names(m)), "Canonical M2 alignment metadata is incomplete; regenerate 4.5 metadata.")
  spec <- m$model_specification[[jm_alignment_model]]
  jm_assert(m$schema_version == 3L && m$cohort == "All Cancer Cohort" &&
    m$outcome == "fi_score_nocancer" && m$primary_model_id == "M1_primary_spline" &&
    m$spline_reference_model_id == "M1_primary_spline" && m$spline_df == 3L &&
    identical(m$primary_covariates, jm_primary_covars) &&
    identical(spec$covars, jm_covars) && identical(spec$basis_key, "adjusted_spline") &&
    identical(spec$time_structure, "natural spline") && identical(spec$spline_df, 3L) &&
    identical(spec$random, "(1 + .random_time | id)") && identical(spec$matching_set_random, FALSE) &&
    identical(m$random_clock, "attained_age") && m$random_time_center == 60 && m$random_time_scale == 4,
    "Canonical metadata does not describe the current overall M2 df-3/attained-age specification.")
  jm_assert(length(m$knots) == 2L && length(m$boundary_knots) == 2L &&
    all(is.finite(c(m$knots, m$boundary_knots))) &&
    all(diff(c(m$boundary_knots[1], m$knots, m$boundary_knots[2])) > 0), "Invalid shared M1-reference spline knots/boundaries.")
  sc <- m$model_scaling[[jm_alignment_model]]
  jm_need(sc, c("column", "center", "scale"), "M2 scaling")
  jm_assert(identical(sc$column, c(paste0("S", 1:3), jm_diet_covars)) &&
    all(is.finite(sc$center)) && all(is.finite(sc$scale) & sc$scale > 0), "Invalid M2 spline/dietary scaling.")
  if (is.null(matched)) return(invisible(TRUE))
  jm_need(matched, c(jm_covars, "age_at_cycle", "index_age"), "M2 matched frame")
  jm_assert(is.logical(matched$post_own_cancer) && !anyNA(matched$post_own_cancer), "post_own_cancer must be a complete logical flag.")
  primary <- matched[matched$Cohort == "All Cancer Cohort" & !matched$post_own_cancer &
    stats::complete.cases(matched[c("Age_Centered", "fi_score_nocancer", jm_primary_covars)]), ]
  t <- jm_numeric(primary$Age_Centered)
  check <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance=1e-8))
  jm_assert(length(t)>0 && all(is.finite(t)) && check(m$knots, stats::quantile(t,c(1/3,2/3))) &&
    check(m$boundary_knots, range(t)) && identical(as.integer(m$n_obs), as.integer(nrow(primary))) &&
    identical(as.integer(m$n_id), as.integer(dplyr::n_distinct(primary$id))),
    "Canonical metadata does not match the current primary spline reference frame. Nemo must rerun 4.5.")
  d <- primary[stats::complete.cases(primary[jm_covars]), ]
  jm_assert(nrow(d)>1 && all(vapply(d[jm_diet_covars], is.numeric, logical(1))), "Invalid M2 complete-case frame/dietary types.")
  values <- cbind(as.data.frame(unclass(splines::ns(jm_numeric(d$Age_Centered),
    knots=m$knots, Boundary.knots=m$boundary_knots))), d[jm_diet_covars])
  jm_assert(all(is.finite(as.matrix(values))) && check(sc$center,colMeans(values)) &&
    check(sc$scale,apply(values,2,stats::sd)), "Stored M2 scaling differs from the current M2 complete-case frame.")
  jm_assert(all(is.finite(d$age_at_cycle)) && all(is.finite(d$index_age)) &&
    max(abs(d$age_at_cycle-d$index_age-d$Age_Centered)) < 1e-8, "Attained-age and index-relative clocks disagree.")
  invisible(TRUE)
}

jm_assignment_ledger <- function(d, assignments) {
  d <- as.data.frame(d)
  for (v in c("id", "Cohort", "role", "Group", "match_set", "cycle")) d[[v]] <- as.character(d[[v]])
  jm_assert(!anyNA(d[c("id", "Cohort", "role", "Group", "match_set", "cycle")]), "Missing assignment identifiers.")
  jm_assert(all(d$Cohort == "All Cancer Cohort"), "Unexpected cohort in overall matched input.")
  jm_assert(all(d$role %in% c("Case", "Control")) &&
              all(d$Group == ifelse(d$role == "Case", "Cancer Case", "Control")), "Group/role mismatch.")
  d$trajectory_id <- paste(d$Cohort, d$match_set, d$id, d$role, sep = "__")
  jm_assert(!anyDuplicated(paste(d$trajectory_id, d$cycle)), "Duplicate assignment-cycle observations.")
  numeric_cols <- intersect(c("index_date", "index_age", "worked_rtmnyr", "Age_Centered",
                              "fi_score_nocancer", "index_age_z", "age_at_cycle", "cancer_dateca", jm_diet_covars), names(d))
  for (v in numeric_cols) {
    value <- jm_numeric(d[[v]])
    jm_assert(!any(!is.na(d[[v]]) & is.na(value)), paste("Unparseable numeric/date field:",v))
    d[[v]] <- value
  }
  jm_assert(all(is.finite(d$index_date)), "Missing or invalid index date.")
  jm_consistent(d, "id", "cancer_dateca", "cancer dates")
  jm_consistent(d, "trajectory_id", c("index_date", "index_age", jm_covars), "assignment fields")
  # Include NA in the distinct count for fields that must be invariant, so a
  # partly missing baseline covariate cannot create multiple survival records.
  fields <- c("trajectory_id", "id", "Cohort", "Group", "role", "match_set",
              "index_date", "index_age", "cancer_dateca", jm_covars)
  observed_assignments <- unique(d[fields])
  jm_assert(!anyDuplicated(observed_assignments$trajectory_id),
            "Assignment has partially missing/inconsistent baseline fields.")
  a <- as.data.frame(assignments)
  jm_need(a, fields, "Complete assignment ledger")
  for (v in c("trajectory_id", "id", "Cohort", "Group", "role", "match_set")) a[[v]] <- as.character(a[[v]])
  jm_assert(!anyNA(a[c("trajectory_id", "id", "Cohort", "Group", "role", "match_set")]) &&
              !anyDuplicated(a$trajectory_id), "Invalid complete assignment ledger keys.")
  jm_assert(all(a$Cohort == "All Cancer Cohort") &&
              all(a$role %in% c("Case", "Control")) &&
              all(a$Group == ifelse(a$role == "Case", "Cancer Case", "Control")),
            "Invalid assignment ledger cohort or role.")
  for (v in c("index_date", "index_age", "index_age_z", "cancer_dateca", jm_diet_covars)) {
    value <- jm_numeric(a[[v]])
    jm_assert(!any(!is.na(a[[v]]) & is.na(value)), paste("Unparseable ledger field:", v))
    a[[v]] <- value
  }
  jm_assert(all(is.finite(a$index_date)) && all(is.finite(a$index_age)), "Invalid ledger index date/age.")
  jm_consistent(a, "id", "cancer_dateca", "ledger cancer dates")
  a <- a[unique(c(fields, grep("^base_", names(a), value = TRUE)))]
  jm_assert(nrow(dplyr::anti_join(observed_assignments, a, by = fields)) == 0,
            "Observed assignments disagree with complete assignment ledger.")
  a$case_priority <- a$role == "Case"
  a$match_set_num <- jm_numeric(a$match_set)
  a <- a[order(a$id, -as.integer(a$case_priority), a$index_date,
               a$match_set_num, a$match_set, na.last = TRUE), ]
  retained <- a[!duplicated(a$id), ]
  a$selected <- a$trajectory_id %in% retained$trajectory_id
  list(long = d, all = a, retained = retained)
}

# Audit all selected assignments before FI/M2 exclusions. Conflicts are retained
# as flags rather than resolved by selecting an arbitrary death date.
jm_date_audit <- function(retained, panel, admin_date=jm_admin_date) {
  jm_need(panel,c("id","dtdth"),"Death-source panel")
  raw <- data.frame(id=as.character(panel$id), original=as.character(panel$dtdth),
                    value=jm_numeric(panel$dtdth))
  raw$invalid <- !is.na(raw$original) & (!is.finite(raw$value) | raw$value<=0)
  death <- raw |> dplyr::group_by(id) |> dplyr::summarise(
    source_present=TRUE, invalid_death_source=any(invalid),
    n_death_dates=dplyr::n_distinct(value[is.finite(value) & value>0]),
    death_values=paste(sort(unique(value[is.finite(value) & value>0])),collapse=";"),
    dtdth=if(n_death_dates==1) unique(value[is.finite(value) & value>0])[1] else NA_real_, .groups="drop")
  s <- as.data.frame(dplyr::left_join(retained,death,by="id"))
  s$missing_death_source <- is.na(s$source_present)
  s$invalid_death_source[is.na(s$invalid_death_source)] <- FALSE
  s$conflicting_death_dates <- !is.na(s$n_death_dates) & s$n_death_dates>1
  s$same_month_death <- !is.na(s$dtdth) & s$dtdth==s$index_date
  s$death_before_index <- !is.na(s$dtdth) & s$dtdth<s$index_date
  s$index_on_after_admin <- s$index_date>=admin_date
  cancer_end <- ifelse(s$role=="Control" & !is.na(s$cancer_dateca),s$cancer_dateca,Inf)
  s$nonpositive_followup <- pmin(admin_date,cancer_end,ifelse(is.na(s$dtdth),Inf,s$dtdth))<=s$index_date
  s$control_cancer_at_before_index <- s$role=="Control" & !is.na(s$cancer_dateca) & s$cancer_dateca<=s$index_date
  s$case_index_mismatch <- s$role=="Case" & (is.na(s$cancer_dateca) | s$cancer_dateca!=s$index_date)
  s$m2_complete_covariates <- stats::complete.cases(s[jm_covars]) &
    apply(s[c("index_age_z",jm_diet_covars)],1,function(z) all(is.finite(z)))
  flags <- c("same_month_death","death_before_index","conflicting_death_dates",
    "invalid_death_source","missing_death_source","index_on_after_admin",
    "nonpositive_followup","control_cancer_at_before_index","case_index_mismatch")
  s$n_date_issues <- rowSums(s[flags]); s$any_date_issue <- s$n_date_issues>0
  summary <- dplyr::bind_rows(lapply(c(flags,"any_date_issue"),function(flag)
    s |> dplyr::group_by(Group) |> dplyr::summarise(issue=flag,
      affected=sum(.data[[flag]]),denominator=dplyr::n(),
      affected_with_complete_M2=sum(.data[[flag]] & m2_complete_covariates),.groups="drop")))
  availability <- dplyr::bind_rows(lapply(jm_covars,function(v)
    s |> dplyr::group_by(Group) |> dplyr::summarise(variable=v,
      missing_or_nonfinite=sum(if(is.numeric(.data[[v]])) !is.finite(.data[[v]]) else is.na(.data[[v]])),
      denominator=dplyr::n(),.groups="drop")))
  overlap <- s |> dplyr::group_by(Group,same_month_death,index_on_after_admin,m2_complete_covariates) |>
    dplyr::summarise(participants=dplyr::n(),.groups="drop")
  list(records=s,summary=summary,overlap=overlap,m2_availability=availability,flags=flags)
}

jm_deaths <- function(panel) {
  d <- data.frame(id = as.character(panel$id), dtdth = jm_numeric(panel$dtdth))
  jm_assert(!any(!is.na(panel$dtdth) & is.na(d$dtdth)), "Unparseable death date.")
  jm_assert(!anyNA(d$id) && all(is.na(d$dtdth) | (is.finite(d$dtdth) & d$dtdth > 0)), "Invalid death date/participant identifier.")
  jm_consistent(d, "id", "dtdth", "death dates")
  d |>
    dplyr::group_by(id) |>
    dplyr::summarise(dtdth = if (all(is.na(dtdth))) NA_real_ else unique(dtdth[!is.na(dtdth)])[1], .groups = "drop")
}

jm_survival <- function(retained, deaths, admin_date = jm_admin_date) {
  jm_assert(all(retained$id %in% deaths$id), "Some retained participants are absent from the death-source panel.")
  s <- dplyr::left_join(retained, deaths, by = "id")
  invalid_death <- !is.na(s$dtdth) & s$dtdth <= s$index_date
  jm_assert(!any(invalid_death), "Death on/before index: same-month values require a date-precision policy; earlier values require source review.",
            s[invalid_death, ], "invalid_death_at_or_before_index")
  jm_assert(!any(s$role == "Control" & !is.na(s$cancer_dateca) & s$cancer_dateca <= s$index_date),
            "Retained control was not cancer-free at index.")
  jm_assert(all(s$role != "Case" | (!is.na(s$cancer_dateca) & s$index_date == s$cancer_dateca)),
            "Case index differs from own diagnosis date.")
  cancer_end <- ifelse(s$role == "Control" & !is.na(s$cancer_dateca), s$cancer_dateca, Inf)
  s$censor_date <- pmin(admin_date, cancer_end)
  s$cancer_death_tie <- is.finite(cancer_end) & !is.na(s$dtdth) &
    s$dtdth == cancer_end & cancer_end <= admin_date
  # Death on the administrative endpoint is an event; cancer-date ties censor.
  s$death_event <- as.integer(!is.na(s$dtdth) & s$dtdth <= admin_date & s$dtdth < cancer_end)
  s$surv_end_date <- ifelse(s$death_event == 1L, s$dtdth, s$censor_date)
  s$endpoint_reason <- ifelse(s$death_event == 1L, "death",
                              ifelse(cancer_end <= admin_date, "own_cancer", "administrative"))
  s$surv_time <- (s$surv_end_date - s$index_date) / 12
  invalid_duration <- !is.finite(s$surv_time) | s$surv_time <= 0
  jm_assert(!any(invalid_duration), "Nonpositive index-to-endpoint duration.",
            s[invalid_duration, ], "invalid_followup_duration")
  jm_assert(!anyDuplicated(s$id), "Duplicate survival contribution.")
  as.data.frame(s)
}

jm_prepare <- function(matched, panel, assignments) {
  ledger <- jm_assignment_ledger(matched, assignments)
  survival <- jm_survival(ledger$retained, jm_deaths(panel))
  d <- dplyr::inner_join(ledger$long, survival[c("trajectory_id", "surv_end_date", "endpoint_reason", "surv_time")], by = "trajectory_id")
  log <- list()
  record <- function(stage, x) {
    z <- x |>
      dplyr::group_by(Group) |>
      dplyr::summarise(rows = dplyr::n(), participants = dplyr::n_distinct(id),
                       assignments = dplyr::n_distinct(trajectory_id), .groups = "drop")
    z$stage <- stage
    log[[length(log) + 1L]] <<- z
  }
  record("full_matched_assignments", ledger$all)
  record("selected_assignments_before_FI", ledger$retained)
  record("full_matched_outcome_rows", ledger$long)
  record("selected_assignment", d)
  jm_assert(is.logical(d$post_own_cancer) && !anyNA(d$post_own_cancer), "Invalid post_own_cancer flag.")
  dated <- is.finite(d$worked_rtmnyr) & is.finite(d$Age_Centered)
  jm_assert(all(abs(d$Age_Centered[dated] - (d$worked_rtmnyr[dated] - d$index_date[dated])/12) < 1e-8),
            "Longitudinal and survival clocks do not share index date/month units.")
  expected_flag <- d$role == "Control" & !is.na(d$cancer_dateca) &
    d$cancer_dateca > d$index_date & !is.na(d$worked_rtmnyr) & d$worked_rtmnyr >= d$cancer_dateca
  jm_assert(identical(d$post_own_cancer, expected_flag), "Own-cancer flag disagrees with dates.")
  # Deliberate endpoint exclusions are logged; they are not passed to the model.
  keep <- dated & !d$post_own_cancer & d$worked_rtmnyr <= d$surv_end_date &
    !(d$endpoint_reason == "own_cancer" & d$worked_rtmnyr >= d$surv_end_date)
  excluded_rows <- d[!keep, c("id", "trajectory_id", "cycle", "worked_rtmnyr", "surv_end_date", "post_own_cancer")]
  d <- d[keep, ]
  record("within_followup", d)
  d <- d[is.finite(d$fi_score_nocancer) & stats::complete.cases(d[jm_covars]) &
    apply(d[c("index_age_z",jm_diet_covars)],1,function(z) all(is.finite(z))), ]
  jm_assert(all(d$fi_score_nocancer >= 0 & d$fi_score_nocancer <= 1), "FI outside [0,1].")
  record("M2_complete_case", d)
  excluded_people <- survival[!survival$id %in% d$id, ]
  survival <- survival[survival$id %in% d$id, ]
  jm_assert(nrow(survival) > 0 && length(unique(survival$Group)) == 2, "No two-arm JM cohort after exclusions.")
  jm_assert(all(is.finite(d$age_at_cycle)) &&
    max(abs(d$age_at_cycle-d$index_age-d$Age_Centered)) < 1e-8, "Attained-age clock mismatch.")
  d <- d[order(d$id, d$Age_Centered), ]
  survival <- survival[order(survival$id), ]
  for (v in c("Group", jm_factors)) {
    lev <- if (v == "Group") c("Control", "Cancer Case") else if(is.factor(d[[v]])) levels(droplevels(d[[v]])) else sort(unique(as.character(d[[v]])))
    jm_assert(length(lev) >= 2, paste("Insufficient factor levels for M2:", v))
    d[[v]] <- factor(d[[v]], levels = lev)
    survival[[v]] <- factor(survival[[v]], levels = lev)
  }
  id_levels <- survival$id
  d$id <- factor(d$id, levels = id_levels)
  survival$id <- factor(survival$id, levels = id_levels)
  jm_assert(identical(as.character(unique(d$id)), as.character(survival$id)), "Longitudinal/survival ordering mismatch.")
  jm_assert(all(d$Age_Centered <= d$surv_time + 1e-8), "FI after survival endpoint.")
  list(long = as.data.frame(d), survival = as.data.frame(survival),
       assignments = ledger$all, retained = ledger$retained,
       filter_log = dplyr::bind_rows(log), excluded_rows = excluded_rows,
       excluded_people = excluded_people)
}

jm_support <- function(d, min_cases = 50L, min_controls = 250L) {
  breaks <- c(-Inf, seq(-20, 20, 2), Inf)
  bins <- cut(d$Age_Centered, breaks, labels = FALSE, right = TRUE)
  tab <- data.frame(bin = seq_len(length(breaks)-1L), lo = head(breaks,-1), hi = tail(breaks,-1))
  for (g in c("Control", "Cancer Case")) {
    tab[[g]] <- vapply(tab$bin, function(k) length(unique(d$id[bins == k & d$Group == g])), integer(1))
  }
  tab$support_ok <- tab[["Control"]] >= min_controls & tab[["Cancer Case"]] >= min_cases
  ref <- which(tab$lo == -2 & tab$hi == 0)
  window <- c(NA_real_, NA_real_)
  if (tab$support_ok[ref]) {
    left <- right <- ref
    while (left > 1 && tab$support_ok[left-1]) left <- left-1
    while (right < nrow(tab) && tab$support_ok[right+1]) right <- right+1
    window <- c(max(-20, tab$lo[left]), min(20, tab$hi[right]))
  }
  list(table = tab, window = window,
       theta_supported = all(is.finite(window)) && window[1] <= -8 && window[2] >= 8)
}

jm_extrapolation <- function(p, metadata) {
  ranges <- p$long |>
    dplyr::group_by(id) |>
    dplyr::summarise(first_FI = min(Age_Centered), last_FI = max(Age_Centered),
                     n_FI = dplyr::n(), n_pre = sum(Age_Centered < 0), n_post = sum(Age_Centered >= 0), .groups = "drop")
  z <- dplyr::left_join(p$survival[c("id", "Group", "surv_time", "death_event")], ranges, by = "id")
  z$years_after_last_FI <- pmax(0, z$surv_time - pmax(0, z$last_FI))
  z$years_before_first_FI <- pmin(z$surv_time, pmax(0, z$first_FI))
  z$beyond_canonical_boundary <- pmax(0, z$surv_time - metadata$boundary_knots[2])
  z$beyond_cohort_FI <- pmax(0, z$surv_time - max(p$long$Age_Centered))
  z
}

jm_balance <- function(matched, prepared) {
  # Each original assignment is a unit in the original matched population;
  # each retained participant is one unit in the selected populations.
  a <- jm_assignment_ledger(matched, prepared$assignments)
  visits <- a$long |>
    dplyr::group_by(trajectory_id) |>
    dplyr::summarise(n_observed_FI = sum(is.finite(fi_score_nocancer)),
      preindex_FI = { ii <- which(is.finite(fi_score_nocancer) & Age_Centered < 0)
                     if (length(ii)) fi_score_nocancer[ii[which.max(Age_Centered[ii])]] else NA_real_ }, .groups = "drop")
  base <- dplyr::left_join(a$all, visits, by = "trajectory_id")
  base$n_observed_FI[is.na(base$n_observed_FI)] <- 0L
  base$index_year <- 1900 + base$index_date / 12
  if (!"base_pckgr" %in% names(base) && "base_pckgr" %in% names(a$long)) {
    smoking <- unique(a$long[c("trajectory_id", "base_pckgr")])
    jm_assert(!anyDuplicated(smoking$trajectory_id), "Inconsistent baseline pack-years.")
    base <- dplyr::left_join(base, smoking, by = "trajectory_id")
  }
  populations <- list(matched_assignments = base,
                      selected_assignments = base[base$selected, ],
                      JM_complete_case = base[base$trajectory_id %in% prepared$survival$trajectory_id, ])
  rows <- list()
  for (pop in names(populations)) {
    d <- populations[[pop]]
    vars <- intersect(c("index_age", "index_year", jm_factors, jm_diet_covars, "preindex_FI", "n_observed_FI"), names(d))
    for (v in vars) {
      categorical <- v %in% c(jm_factors, "base_pckgr")
      levels_v <- if (categorical) sort(unique(as.character(base[[v]][!is.na(base[[v]])]))) else "numeric"
      for (lev in levels_v) {
        x <- if (categorical) as.numeric(as.character(d[[v]]) == lev) else jm_numeric(d[[v]])
        case <- x[d$Group == "Cancer Case"]; ctrl <- x[d$Group == "Control"]
        m1 <- mean(case,na.rm=TRUE); m0 <- mean(ctrl,na.rm=TRUE)
        denom <- if (categorical) sqrt((m1*(1-m1)+m0*(1-m0))/2) else sqrt((stats::var(case,na.rm=TRUE)+stats::var(ctrl,na.rm=TRUE))/2)
        smd <- if (is.finite(denom) && denom > 0) (m1-m0)/denom else if (isTRUE(m1 == m0)) 0 else NA_real_
        rows[[length(rows)+1L]] <- data.frame(population=pop, variable=v, level=lev,
          case_mean_or_proportion=m1, control_mean_or_proportion=m0, smd=smd,
          n_case=sum(!is.na(case)), n_control=sum(!is.na(ctrl)),
          missing_case=sum(is.na(case)), missing_control=sum(is.na(ctrl)))
      }
    }
  }
  dplyr::bind_rows(rows)
}

# A formula-local closure is serialized with the model. It contains only the
# fixed basis constants, not analytic data, and reevaluates splines at NEW times.
jm_formula <- function(metadata) {
  jm_validate_metadata(metadata)
  env <- new.env(parent = baseenv())
  env$knots <- metadata$knots
  env$boundary <- metadata$boundary_knots
  env$centers <- metadata$model_scaling[[jm_alignment_model]]$center[1:3]
  env$scales <- metadata$model_scaling[[jm_alignment_model]]$scale[1:3]
  basis <- function(t) {
    b <- splines::ns(t, knots = knots, Boundary.knots = boundary)
    b <- matrix(as.numeric(b), nrow = length(t), ncol = 3L)
    b <- sweep(sweep(b, 2, centers, "-"), 2, scales, "/")
    colnames(b) <- paste0("S", seq_len(ncol(b)))
    b
  }
  environment(basis) <- env
  env$jm_time_basis <- basis
  stats::as.formula(paste("fi_score_nocancer ~ Group * jm_time_basis(Age_Centered) +",paste(jm_covars,collapse=" + ")), env=env)
}

jm_random_formula <- function() ~ 1 + I((index_age + Age_Centered - 60)/4)

# Keep raw dietary columns in explicit companions; reject any second scaling.
jm_scale_data <- function(p, metadata) {
  jm_validate_metadata(metadata)
  sc <- metadata$model_scaling[[jm_alignment_model]]
  transform <- function(d) {
    jm_assert(!any(paste0(jm_diet_covars,"_raw") %in% names(d)), "M2 dietary covariates already scaled.")
    for(v in jm_diet_covars) {
      j <- match(v,sc$column); d[[paste0(v,"_raw")]] <- d[[v]]
      d[[v]] <- (d[[v]]-sc$center[j])/sc$scale[j]
    }
    d
  }
  p$long <- transform(p$long); p$survival <- transform(p$survival)
  p$scaling <- sc
  p
}

jm_check_clock_basis <- function(metadata) {
  jm_assert(requireNamespace("JMbayes2", quietly = TRUE), "Install JMbayes2 before preparing/fitting this workflow.")
  ns <- asNamespace("JMbayes2")
  fun <- get("design_matrices_functional_forms", ns)
  form <- jm_formula(metadata)
  # Two histories with negative times; all covariates remain fixed within id.
  d <- data.frame(id = rep(c("A","B"), each=4), Age_Centered=rep(c(-8,-4,1,5),2),
                  Group=factor(rep(c("Control","Cancer Case"),each=4),levels=c("Control","Cancer Case")),
                  index_age=c(rep(65,4),rep(80,4)), index_age_z=0, base_race=factor(rep(c("r1","r2"),4)),
                  base_marital=factor(rep(c("m1","m2"),4)), base_living=factor(rep(c("l1","l2"),4)))
  # Use constant factors within each subject so carrying baseline fields is valid.
  for (v in jm_factors) d[[v]] <- factor(rep(c("a","b"), each=4))
  for (v in jm_diet_covars) d[[v]] <- 0
  rhs <- stats::delete.response(stats::terms(form))
  times <- matrix(c(0.25,2.5,6,0.5,3.5,7),nrow=2,byrow=TRUE)
  args <- list(time=times, terms=list(rhs), data=d, timeVar="Age_Centered", idVar="id",
               idT=c("A","B"), Fun_Forms=list(c("value","slope")), Xbar=NULL,
               eps=list(slope=0.001), direction=list("both"), zero_ind=NULL,
               time_window=list(NULL), standardise=list(NULL), IE_time=list(NULL))
  # API arrays are per outcome, not per association.
  args$eps <- list(0.001)
  matrices <- do.call(fun,args)
  direct <- d[rep(c(1,5),each=3), ]
  direct$Age_Centered <- c(t(times))
  X <- stats::model.matrix(rhs,direct)
  plus <- minus <- direct
  plus$Age_Centered <- direct$Age_Centered+0.001
  minus$Age_Centered <- direct$Age_Centered-0.001
  DX <- (stats::model.matrix(rhs,plus)-stats::model.matrix(rhs,minus))/0.002
  # Package output is one list per outcome, containing named functional forms.
  value <- matrices[[1]][["value"]]
  slope <- matrices[[1]][["slope"]]
  jm_assert(is.matrix(value) && is.matrix(slope), "Unsupported JMbayes2 functional-matrix schema.")
  ve <- max(abs(value-X)); se <- max(abs(slope-DX))
  jm_assert(is.finite(ve) && is.finite(se) && ve < 1e-8 && se < 1e-7,
            "JMbayes2 does not correctly reevaluate the index-origin spline/history.")
  args$terms <- list(stats::terms(jm_random_formula()))
  random <- do.call(fun,args)[[1]]
  Z <- stats::model.matrix(jm_random_formula(),direct)
  DZ <- cbind(0,rep(1/4,nrow(direct)))
  ze <- max(abs(random$value-Z)); dze <- max(abs(random$slope-DZ))
  jm_assert(ze<1e-8 && dze<1e-7, "JMbayes2 attained-age random value/slope contract failed.")
  code <- paste(deparse(get("jm.default",ns)),collapse=" ")
  jm_assert(grepl('Time_left <- Time_start <- trunc_Time <- rep(0',code,fixed=TRUE),
            "JMbayes2 right-censoring origin changed; review its survival integration before fitting.")
  data.frame(package_version=as.character(utils::packageVersion("JMbayes2")),
             value_max_error=ve, slope_max_error=se,
             random_value_max_error=ze, random_slope_max_error=dze,
             negative_history_verified=TRUE, survival_origin=0)
}

jm_pilot_subset <- function(p, seed = 20260703L, fraction = 0.25) {
  set.seed(seed)
  controls <- as.character(p$survival$id[p$survival$Group == "Control"])
  cases <- as.character(p$survival$id[p$survival$Group == "Cancer Case"])
  keep <- c(cases, sample(controls, ceiling(length(controls)*fraction)))
  p$long <- droplevels(p$long[as.character(p$long$id) %in% keep, ])
  p$survival <- droplevels(p$survival[as.character(p$survival$id) %in% keep, ])
  jm_assert(identical(as.character(unique(p$long$id)),as.character(p$survival$id)), "Pilot ID order mismatch.")
  for (v in jm_factors) jm_assert(nlevels(p$long[[v]]) >= 2, paste("Pilot has insufficient levels:",v))
  p
}

jm_fit_components <- function(p, formula, run_dir) {
  jm_require_resolved_dates()
  jm_assert(requireNamespace("nlme",quietly=TRUE) && requireNamespace("survival",quietly=TRUE), "nlme and survival are required.")
  d <- p$long; s <- p$survival
  X <- stats::model.matrix(formula, d)
  jm_assert(qr(X)$rank == ncol(X), "M2 longitudinal fixed-effect matrix is rank deficient.")
  attempts <- list(); accepted <- NULL; chosen <- NULL
  for (rung in c("correlated_intercept_slope", "uncorrelated_intercept_slope")) {
    warnings <- character(); error <- NULL
    random <- if (rung == "correlated_intercept_slope") list(id=nlme::pdSymm(jm_random_formula())) else list(id=nlme::pdDiag(jm_random_formula()))
    start <- proc.time()[[3]]
    fit <- tryCatch(withCallingHandlers(
      nlme::lme(fixed=formula, random=random, data=d, method="REML", na.action=stats::na.fail,
                control=nlme::lmeControl(opt="optim", maxIter=200, msMaxIter=400,
                                        niterEM=50, returnObject=FALSE)),
      warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")}),
      error=function(e) {error <<- conditionMessage(e); NULL})
    valid <- FALSE; eigen_ratio <- NA_real_
    if (!is.null(fit)) {
      D <- as.matrix(nlme::getVarCov(fit, type="random.effects"))
      eig <- eigen(D, symmetric=TRUE, only.values=TRUE)$values
      eigen_ratio <- min(eig)/max(eig)
      valid <- !length(warnings) && all(is.finite(nlme::fixef(fit))) &&
        all(is.finite(stats::vcov(fit))) && is.matrix(fit$apVar) &&
        all(is.finite(fit$apVar)) && all(is.finite(eig)) &&
        min(eig) > 1e-10 && eigen_ratio > 1e-8
      # Reject unsupported random-effect dimension rather than falling back.
      valid <- valid && nrow(D) == 2L
      jm_save(fit,file.path(run_dir,paste0("longitudinal_",rung,".rds")))
    }
    attempts[[length(attempts)+1]] <- data.frame(rung=rung, accepted=valid,
      eigen_ratio=eigen_ratio, elapsed_seconds=proc.time()[[3]]-start,
      warnings=paste(warnings,collapse=" | "), error=if(is.null(error)) "" else error)
    jm_write_csv(dplyr::bind_rows(attempts),file.path(run_dir,"component_attempts.csv"))
    if (valid) {accepted <- fit; chosen <- rung; break}
  }
  jm_assert(!is.null(accepted), "No supported random-slope LME passed convergence/variance checks. Components saved; linked fitting stopped.")
  # Component coefficients use canonical M2-scaled spline/dietary columns. JM priors are
  # the installed package defaults on this parameterization and are persisted.
  cox_formula <- stats::as.formula(paste("survival::Surv(surv_time, death_event) ~ Group +",paste(jm_covars,collapse=" + ")))
  XS <- stats::model.matrix(cox_formula,s)
  jm_assert(qr(XS)$rank==ncol(XS), "M2 mortality design is rank deficient.")
  jm_assert(sum(s$death_event)>0, "No mortality events in this run's cohort.")
  warning_cox <- character()
  cox <- withCallingHandlers(survival::coxph(cox_formula,data=s,x=TRUE,model=TRUE,
    ties="efron",na.action=stats::na.fail), warning=function(w) {
      warning_cox <<- c(warning_cox,conditionMessage(w)); invokeRestart("muffleWarning")})
  jm_save(cox,file.path(run_dir,"cox_model.rds"))
  jm_assert(!length(warning_cox) && all(is.finite(stats::coef(cox))) && all(is.finite(stats::vcov(cox))),
            paste("Cox component failed checks:",paste(warning_cox,collapse=" | ")))
  zph <- tryCatch(survival::cox.zph(cox),error=function(e) list(error=conditionMessage(e)))
  km <- survival::survfit(survival::Surv(surv_time,death_event)~Group,data=s)
  components <- list(longitudinal=accepted, survival=cox, zph=zph, km=km, random_structure=chosen,
                      fixed_formula=formula, covariance="nlme model-based; selected cohort retains possible matched-set dependence")
  jm_save(components,file.path(run_dir,"components.rds"))
  capture.output(summary(accepted),file=file.path(run_dir,"longitudinal_summary.txt"))
  capture.output(summary(cox),file=file.path(run_dir,"cox_summary.txt"))
  capture.output(print(zph),file=file.path(run_dir,"cox_ph_diagnostics.txt"))
  components
}

jm_design <- function(p, formula, support) {
  d <- p$long
  reference <- data.frame(index_age_z=0)
  for (v in jm_factors) {
    per_id <- d[!duplicated(d$id),v]
    counts <- table(per_id)
    reference[[v]] <- factor(names(counts)[which.max(counts)],levels=levels(d[[v]]))
  }
  for(v in jm_diet_covars) {
    jm_assert(paste0(v,"_raw") %in% names(d), "Prediction requires raw and scaled M2 dietary fields.")
    reference[[v]] <- mean(d[[paste0(v,"_raw")]])
  }
  raw_reference <- reference
  for(v in jm_diet_covars) {
    j <- match(v,p$scaling$column)
    reference[[v]] <- (reference[[v]]-p$scaling$center[j])/p$scaling$scale[j]
  }
  rhs <- stats::delete.response(stats::terms(formula))
  grid <- function(t,g) {
    nd <- reference[rep(1,length(t)),,drop=FALSE]
    nd$Age_Centered <- t
    nd$Group <- factor(g,levels=c("Control","Cancer Case"))
    nd
  }
  matrix_at <- function(t,g) stats::model.matrix(rhs,grid(t,g))
  difference <- function(t) matrix_at(t,"Cancer Case")-matrix_at(t,"Control")
  derivative <- function(t) (difference(t+1e-4)-difference(t-1e-4))/2e-4
  names_beta <- colnames(matrix_at(0,"Control"))
  times <- if (all(is.finite(support$window))) seq(support$window[1],support$window[2],0.25) else numeric()
  empty <- matrix(numeric(),0,length(names_beta),dimnames=list(NULL,names_beta))
  theta <- as.numeric((difference(8)-difference(0))/8-(difference(0)-difference(-8))/8)
  names(theta) <- names_beta
  numeric_theta <- colMeans(derivative(seq(0.025,7.975,0.05)))-colMeans(derivative(seq(-7.975,-0.025,0.05)))
  # Independent basis-level check, followed by estimate-level check at inference.
  jm_assert(max(abs(theta-numeric_theta)) < 1e-5, "Endpoint/derivative theta basis mismatch.")
  list(times=times, control=if(length(times)) matrix_at(times,"Control") else empty,
       case=if(length(times)) matrix_at(times,"Cancer Case") else empty,
       difference=if(length(times)) difference(times) else empty,
       derivative=if(length(times)) derivative(times) else empty,
       theta=theta, numeric_theta=numeric_theta, theta_supported=support$theta_supported,
       reference=reference, raw_reference=raw_reference)
}

jm_component_summaries <- function(model, design) {
  b <- nlme::fixef(model); V <- stats::vcov(model)
  jm_assert(setequal(names(b),names(design$theta)), "Component design coefficient mismatch.")
  summarize <- function(X) {
    X <- X[,names(b),drop=FALSE]
    est <- as.numeric(X%*%b); se <- sqrt(pmax(0,rowSums((X%*%V)*X)))
    data.frame(estimate=est, se=se, lower=est-1.96*se, upper=est+1.96*se)
  }
  curves <- dplyr::bind_rows(lapply(c("control","case","difference","derivative"),function(q) {
    X <- design[[q]]
    if(!nrow(X)) return(NULL)
    cbind(data.frame(quantity=q,time=design$times),summarize(X))
  }))
  theta <- data.frame(quantity="theta", supported=design$theta_supported,
                      estimate=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_)
  if (design$theta_supported) {
    crosscheck <- sum((design$theta-design$numeric_theta)*b[names(design$theta)])
    jm_assert(abs(crosscheck)<1e-6,"Endpoint/derivative component theta differs by >1e-6.")
    theta[1,c("estimate","se","lower","upper")] <- summarize(matrix(design$theta,nrow=1,dimnames=list(NULL,names(design$theta))))
  }
  list(curves=curves,theta=theta,covariance="nlme model-based normal approximation")
}

jm_chain_array <- function(chains, required_names, expected_chains, label) {
  jm_assert(is.list(chains) && length(chains)==expected_chains, paste("Missing/wrong chain count:",label))
  matrices <- lapply(chains,as.matrix)
  jm_assert(length(required_names)>0 && !anyDuplicated(required_names),paste("Invalid required coefficients:",label))
  for (m in matrices) {
    jm_assert(is.numeric(m) && nrow(m)>=4 && !is.null(colnames(m)) && !anyDuplicated(colnames(m)) &&
                setequal(colnames(m),required_names) && all(is.finite(m)),paste("Incomplete/nonfinite posterior coefficients:",label))
  }
  sizes <- vapply(matrices,nrow,integer(1))
  jm_assert(length(unique(sizes))==1L,paste("Unequal posterior chain lengths:",label))
  out <- array(NA_real_,c(sizes[1],expected_chains,length(required_names)),
               dimnames=list(iteration=NULL,chain=as.character(seq_len(expected_chains)),parameter=required_names))
  for(i in seq_along(matrices)) out[,i,] <- matrices[[i]][,required_names,drop=FALSE]
  out
}

jm_posterior <- function(fit, beta_names, gamma_names, structure, expected_chains) {
  jm_assert(is.list(fit$mcmc), "Joint model has no posterior draws.")
  alpha_names <- unlist(lapply(fit$model_data$U_H,colnames),use.names=FALSE)
  expected_alpha_count <- if(structure=="value") 1L else 2L
  jm_assert(length(alpha_names)==expected_alpha_count && all(nzchar(alpha_names)) &&
              any(grepl("value",alpha_names,fixed=TRUE)) &&
              (structure=="value" || any(grepl("slope",alpha_names,fixed=TRUE))),
            "Missing or unexpected current-value/slope association parameters.")
  jm_assert(is.matrix(fit$initial_values$D) && identical(dim(fit$initial_values$D),c(2L,2L)),
    "Expected the supported two-dimensional participant covariance.")
  covariance_names <- c("D[1, 1]","D[2, 1]","D[2, 2]")
  D <- jm_chain_array(fit$mcmc$D,covariance_names,expected_chains,"D")
  zeros <- fit$model_data$ind_zero_D
  constrained <- character()
  if(!is.null(zeros) && length(zeros)) {
    jm_assert(is.matrix(zeros) && ncol(zeros)==2L && nrow(zeros)==1L &&
      identical(as.integer(zeros[1,]),c(1L,2L)), "Unsupported covariance zero constraints.")
    constrained <- "D[2, 1]"
    jm_assert(all(abs(D[,,constrained])<1e-12), "Constrained covariance contains nonzero draws.")
  }
  jm_assert(all(D[,,"D[1, 1]"]>0) && all(D[,,"D[2, 2]"]>0) &&
    all(D[,,"D[1, 1]"]*D[,,"D[2, 2]"]-D[,,"D[2, 1]"]^2>0), "Nonpositive posterior covariance.")
  nb <- ncol(fit$model_data$W0_H)
  jm_assert(length(nb)==1L && is.finite(nb) && nb>0, "Missing baseline-hazard design.")
  out <- list(betas=jm_chain_array(fit$mcmc$betas1,beta_names,expected_chains,"betas1"),
    gammas=jm_chain_array(fit$mcmc$gammas,gamma_names,expected_chains,"gammas"),
    alphas=jm_chain_array(fit$mcmc$alphas,alpha_names,expected_chains,"alphas"),
    sigmas=jm_chain_array(fit$mcmc$sigmas,"sigmas_1",expected_chains,"sigmas"),
    D=D[,,setdiff(covariance_names,constrained),drop=FALSE],
    bs_gammas=jm_chain_array(fit$mcmc$bs_gammas,paste0("bs_gammas_",seq_len(nb)),expected_chains,"bs_gammas"),
    tau_bs_gammas=jm_chain_array(fit$mcmc$tau_bs_gammas,"tau_bs_gammas",expected_chains,"tau_bs_gammas"))
  jm_assert(all(out$sigmas>0) && all(out$tau_bs_gammas>0), "Nonpositive scale/smoothing posterior draws.")
  attr(out,"constrained_parameters") <- data.frame(block=rep("D",length(constrained)),parameter=constrained,
    fixed_value=rep(0,length(constrained)),status=rep("verified_fixed_zero_not_sampled",length(constrained)))
  out
}

jm_diagnostics <- function(draws, profile) {
  jm_assert(requireNamespace("posterior",quietly=TRUE), "posterior package is required for rank-normalized diagnostics.")
  rows <- list()
  for (block in names(draws)) {
    a <- draws[[block]]
    for(j in seq_len(dim(a)[3])) {
      x <- matrix(a[,,j],nrow=dim(a)[1],ncol=dim(a)[2])
      rh <- if(ncol(x)>=2) posterior::rhat(x) else NA_real_
      bulk <- posterior::ess_bulk(x); tail <- posterior::ess_tail(x)
      rh_ok <- is.finite(rh) && rh < 1.05
      ess_ok <- is.finite(bulk) && is.finite(tail) && bulk > 400 && tail > 400
      rows[[length(rows)+1]] <- data.frame(block=block,parameter=dimnames(a)[[3]][j],
          rhat=rh,ess_bulk=bulk,ess_tail=tail,chains=ncol(x),
          pass=profile=="full" && rh_ok && ess_ok)
    }
  }
  dplyr::bind_rows(rows)
}

jm_flatten <- function(a) {
  # Iterations within chain are stacked only AFTER chain-preserving diagnostics.
  matrix(a,nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3],dimnames=list(NULL,dimnames(a)[[3]]))
}
jm_draw_summary <- function(x) {
  c(estimate=mean(x),sd=stats::sd(x),lower=unname(stats::quantile(x,0.025)),upper=unname(stats::quantile(x,0.975)))
}
jm_joint_summaries <- function(draws,design,comparator) {
  b <- jm_flatten(draws$betas)
  summarize <- function(X) {
    X <- X[,colnames(b),drop=FALSE]
    if(!nrow(X)) return(data.frame())
    # Rowwise evaluation avoids a large grid x draw temporary for full runs.
    as.data.frame(t(vapply(seq_len(nrow(X)),function(i) jm_draw_summary(as.numeric(b%*%X[i,])),numeric(4))))
  }
  curves <- dplyr::bind_rows(lapply(c("control","case","difference","derivative"),function(q) {
    if(!nrow(design[[q]])) return(NULL)
    cbind(data.frame(quantity=q,time=design$times),summarize(design[[q]]))
  }))
  theta <- data.frame(quantity="theta",supported=design$theta_supported,estimate=NA_real_,sd=NA_real_,
                      lower=NA_real_,upper=NA_real_,sensitivity_difference_vs_longitudinal=NA_real_)
  if(design$theta_supported) {
    td <- as.numeric(b%*%design$theta[colnames(b)])
    check <- as.numeric(b%*%(design$theta-design$numeric_theta)[colnames(b)])
    jm_assert(max(abs(check))<1e-6,"Posterior endpoint/derivative theta differs by >1e-6.")
    theta[1,c("estimate","sd","lower","upper")] <- as.list(jm_draw_summary(td))
    theta$sensitivity_difference_vs_longitudinal <- mean(td)-comparator$theta$estimate
  }
  a <- jm_flatten(draws$alphas)
  associations <- dplyr::bind_rows(lapply(seq_len(ncol(a)),function(j) {
    cbind(data.frame(parameter=colnames(a)[j],scale=if(grepl("slope",colnames(a)[j],fixed=TRUE)) "HR per 0.1 FI/year" else "HR per 0.1 FI"),
          as.data.frame(as.list(setNames(jm_draw_summary(a[,j]),paste0("alpha_",names(jm_draw_summary(a[,j])))))),
          as.data.frame(as.list(setNames(jm_draw_summary(exp(0.1*a[,j])),paste0("hr_",names(jm_draw_summary(a[,j])))))))
  }))
  g <- jm_flatten(draws$gammas)
  gi <- grep("^Group",colnames(g))
  jm_assert(length(gi)==1,"Expected one Group mortality coefficient.")
  group_hr <- as.data.frame(as.list(jm_draw_summary(exp(g[,gi]))))
  list(curves=curves,theta=theta,associations=associations,group_mortality_hr=group_hr)
}

jm_trace_data <- function(draws) {
  # Persist all monitored draws; no thinning beyond the declared MCMC setting.
  dplyr::bind_rows(lapply(names(draws),function(block) {
    a <- draws[[block]]
    z <- expand.grid(iteration=seq_len(dim(a)[1]),chain=seq_len(dim(a)[2]),parameter=dimnames(a)[[3]],stringsAsFactors=FALSE)
    z$value <- as.numeric(a); z$block <- block
    z
  }))
}

# Separate adapter allows orchestration/checkpoint tests with synthetic model
# objects, without invoking a statistical fitter or changing package namespaces.
jm_configuration <- function(metadata, profile, provenance, project_dir) {
  list(schema_version=2L, alignment_model_id=jm_alignment_model,profile=profile,seed=20260703L,
    selected_cohort="one matched case assignment if available, otherwise earliest control assignment",
    time_origin="index; survival conditional on alive at index; pre-index FI retained as history",
    admin_date=jm_admin_date,own_cancer_censoring=TRUE,date_policy=jm_date_policy(),
    fixed_formula=jm_formula(metadata),random_formula=jm_random_formula(),
    covariates=jm_covars,knots=metadata$knots,boundary_knots=metadata$boundary_knots,
    spline_reference_model_id=metadata$spline_reference_model_id,
    canonical_scaling=metadata$model_scaling[[jm_alignment_model]],
    baseline_hazard=jm_hazard_control,theta_times=c(-8,0,8),
    priors="Installed JMbayes2 defaults; actual priors and internal standardization saved per fit. No FI x 29 conversion.",
    mcmc=if(profile=="full") list(n_chains=3L,n_iter=15000L,n_burnin=5000L,n_thin=5L) else
      list(n_chains=1L,n_iter=1500L,n_burnin=500L,n_thin=2L),
    input_hashes=provenance$hashes,
    script_hashes=tools::md5sum(file.path(project_dir,"Code","2_data_analysis",
      c("8.0_joint_model_functions.R","8.1_joint_model_overall.R"))),
    software_versions=vapply(c("JMbayes2","nlme","survival","posterior"),
      function(pkg) as.character(utils::packageVersion(pkg)),character(1)))
}

jm_fit_details <- function(fit) {
  for(nm in names(jm_hazard_control)) jm_assert(isTRUE(all.equal(fit$control[[nm]],jm_hazard_control[[nm]])),
    paste("Realized baseline hazard differs from declared control:",nm))
  jm_assert(!is.null(fit$control$knots),"Joint model did not retain realized hazard knots.")
  scaling <- fit[c("Wlong_bar","Wlong_sds","Wlong_std","W_bar","W_sds","W_std")]
  jm_assert(all(vapply(scaling,function(z) !is.null(z),logical(1))) && length(scaling)==6L,
    "Joint model internal standardization was not retained.")
  list(alignment_model_id=jm_alignment_model,control=fit$control,
    priors=fit$priors,internal_standardization=scaling,
    functional_design_columns=lapply(fit$model_data$U_H,colnames),
    covariance_zero_constraints=fit$model_data$ind_zero_D,
    package_version=as.character(utils::packageVersion("JMbayes2")))
}

jm_fit_criteria <- function(fit, inference_status) {
  x <- fit$fit_stats$marginal
  fields <- c("DIC","WAIC","LPML")
  jm_assert(all(fields %in% names(x)) && all(vapply(x[fields],
    function(v) is.numeric(v) && length(v)==1L && is.finite(v),logical(1))), "Missing/nonfinite marginal fit statistics.")
  data.frame(criterion=fields,value=as.numeric(unlist(x[fields])),likelihood="marginal",
    inference_status=inference_status,automatic_model_selection=FALSE)
}

jm_fit_linked <- function(components, p, config, functional_form) {
  jm_require_resolved_dates()
  JMbayes2::jm(components$survival, components$longitudinal,
    time_var="Age_Centered", data_Surv=p$survival, id_var="id",
    functional_forms=functional_form, control=config$baseline_hazard, n_chains=config$mcmc$n_chains,
    n_iter=config$mcmc$n_iter, n_burnin=config$mcmc$n_burnin,
    n_thin=config$mcmc$n_thin, cores=1L, save_random_effects=FALSE, seed=config$seed)
}

jm_run <- function(project_dir, profile=jm_profile()) {
  profile <- jm_profile(profile)
  run_dir <- jm_new_run(file.path(project_dir,"Results","cancer","data","8.1_joint_model_overall"),profile)
  message("JM output directory: ",run_dir)
  started <- Sys.time()
  status <- list(profile=profile,alignment_model_id=jm_alignment_model,status="preparing",started=started,run_dir=run_dir,
                 inference="Not available; preparation or unreviewed model output")
  save_status <- function() jm_save(status,file.path(run_dir,"run_status.rds"))
  save_status()
  tryCatch({
    x <- jm_validate_inputs(project_dir)
    config <- jm_configuration(x$metadata,profile,x$provenance,project_dir)
    jm_save(config,file.path(run_dir,"configuration.rds"))
    jm_save(x$provenance,file.path(run_dir,"provenance.rds"))
    jm_save(x$metadata,file.path(run_dir,"canonical_M2_alignment_metadata.rds"))
    capture.output(utils::sessionInfo(),file=file.path(run_dir,"sessionInfo.txt"))
    ledger <- jm_assignment_ledger(x$matched,x$assignments)
    audit <- jm_date_audit(ledger$retained,x$panel)
    jm_write_csv(ledger$all,file.path(run_dir,"assignment_selection.csv"))
    jm_write_csv(audit$records,file.path(run_dir,"date_audit.csv"))
    jm_write_csv(audit$summary,file.path(run_dir,"date_audit_summary.csv"))
    jm_write_csv(audit$overlap,file.path(run_dir,"date_issue_overlap.csv"))
    jm_write_csv(audit$m2_availability,file.path(run_dir,"M2_covariate_availability.csv"))
    jm_write_csv(jm_check_clock_basis(x$metadata),file.path(run_dir,"clock_basis_contract.csv"))
    status$selected_participants <- nrow(ledger$retained)
    status$participants_with_date_issues <- sum(audit$records$any_date_issue)
    status$date_policy <- jm_date_policy(); save_status()
    # Every profile records audits first, then stops. There is no runtime switch
    # to exclude early deaths, add an epsilon, or change the administrative date.
    jm_require_resolved_dates()
    jm_assert(!any(audit$records$any_date_issue),"Unresolved date validation issues remain; fitting blocked.")
    p <- jm_prepare(x$matched,x$panel,x$assignments)
    jm_write_csv(p$filter_log,file.path(run_dir,"filter_log.csv"))
    jm_write_csv(p$excluded_rows,file.path(run_dir,"excluded_longitudinal_rows.csv"))
    jm_write_csv(p$excluded_people,file.path(run_dir,"excluded_participants.csv"))
    balance <- jm_balance(x$matched,p)
    jm_write_csv(balance,file.path(run_dir,"cohort_balance.csv"))
    if(profile=="pilot") p <- jm_pilot_subset(p)
    support <- jm_support(p$long)
    formula <- jm_formula(x$metadata)
    config$prediction_window <- support$window
    p <- jm_scale_data(p,x$metadata)
    jm_save(config,file.path(run_dir,"configuration.rds"))
    jm_save(p,file.path(run_dir,"prepared_data.rds"))
    jm_write_csv(p$long,file.path(run_dir,"longitudinal_dataset.csv"))
    jm_write_csv(p$survival,file.path(run_dir,"survival_dataset.csv"))
    jm_write_csv(support$table,file.path(run_dir,"support_by_time.csv"))
    jm_write_csv(jm_extrapolation(p,x$metadata),file.path(run_dir,"mortality_extrapolation.csv"))
    summary <- p$survival |>
      dplyr::group_by(Group) |>
      dplyr::summarise(participants=dplyr::n(),deaths=sum(death_event),
        cancer_censored=sum(endpoint_reason=="own_cancer"),cancer_death_ties=sum(cancer_death_tie),
        median_followup=stats::median(surv_time),.groups="drop")
    jm_write_csv(summary,file.path(run_dir,"survival_summary.csv"))
    # Save actual selected/pilot counts separately from full-cohort selection logs.
    status$participants <- nrow(p$survival); status$observations <- nrow(p$long)
    rm(x)
    if(profile=="prepare") {
      status$status <- "prepared_no_models_fitted"
      status$completed <- Sys.time(); save_status()
      message("Preparation complete; no models fitted. ",run_dir)
      return(invisible(run_dir))
    }
    status$status <- "fitting_components"; save_status()
    components <- jm_fit_components(p,formula,run_dir)
    design <- jm_design(p,formula,support)
    comparator <- jm_component_summaries(components$longitudinal,design)
    jm_save(design,file.path(run_dir,"comparison_design.rds"))
    jm_save(comparator,file.path(run_dir,"longitudinal_only_summaries.rds"))
    jm_write_csv(comparator$curves,file.path(run_dir,"longitudinal_only_curves.csv"))
    jm_write_csv(comparator$theta,file.path(run_dir,"longitudinal_only_theta.csv"))
    status$status <- "components_checkpointed"; save_status()
    structures <- list(value=~value(fi_score_nocancer),value_slope=~value(fi_score_nocancer)+slope(fi_score_nocancer))
    results <- list()
    for(nm in names(structures)) {
      subdir <- file.path(run_dir,nm); dir.create(subdir)
      elapsed <- proc.time()[[3]]
      results[[nm]] <- tryCatch({
        fit <- jm_fit_linked(components,p,config,structures[[nm]])
        jm_save(fit,file.path(subdir,"joint_model.rds"))
        jm_save(fit$priors,file.path(subdir,"priors.rds"))
        jm_save(jm_fit_details(fit),file.path(subdir,"realized_specification.rds"))
        draws <- jm_posterior(fit,names(nlme::fixef(components$longitudinal)),
                              names(stats::coef(components$survival)),nm,config$mcmc$n_chains)
        diag <- jm_diagnostics(draws,profile)
        jm_write_csv(diag,file.path(subdir,"mcmc_diagnostics.csv"))
        jm_write_csv(attr(draws,"constrained_parameters"),file.path(subdir,"constrained_parameters.csv"))
        jm_save(draws,file.path(subdir,"monitored_draws.rds"))
        jm_save(jm_trace_data(draws),file.path(subdir,"trace_data.rds"))
        passed <- profile=="full" && all(diag$pass)
        # Numerical failures cannot produce inferential-looking tables. Pilots
        # may export explicitly non-inferential feasibility summaries.
        if(profile=="pilot" || passed) {
          sums <- jm_joint_summaries(draws,design,comparator)
          sums$inference_status <- if(passed) "numerical_diagnostics_passed_pending_trace_review" else "NON_INFERENTIAL_PILOT"
          criteria <- jm_fit_criteria(fit,sums$inference_status)
          jm_save(sums,file.path(subdir,"posterior_summaries.rds"))
          jm_write_csv(criteria,file.path(subdir,"marginal_fit_statistics.csv"))
          for(q in c("curves","theta","associations","group_mortality_hr")) {
            z <- sums[[q]]; z$inference_status <- rep(sums$inference_status,nrow(z))
            jm_write_csv(z,file.path(subdir,paste0("joint_",q,".csv")))
          }
        }
        capture.output(summary(fit),file=file.path(subdir,"joint_summary_UNREVIEWED.txt"))
        list(status=if(profile=="pilot") "NON_INFERENTIAL_PILOT" else if(passed) "numerical_pass_pending_trace_review" else "failed_convergence_no_inference",
             numerical_pass=passed,elapsed_seconds=proc.time()[[3]]-elapsed)
      },error=function(e) list(status="failed",error=conditionMessage(e),numerical_pass=FALSE,
                               elapsed_seconds=proc.time()[[3]]-elapsed))
      jm_save(results[[nm]],file.path(subdir,"status.rds"))
      jm_save(results,file.path(run_dir,"association_status.rds"))
    }
    status$associations <- results
    status$status <- if(all(vapply(results,function(z) z$status!="failed" && z$status!="failed_convergence_no_inference",logical(1))))
      if(profile=="pilot") "pilot_complete_noninferential" else "numerical_pass_pending_trace_review" else "incomplete_or_failed_no_final_inference"
    status$completed <- Sys.time(); status$elapsed_seconds <- as.numeric(difftime(Sys.time(),started,units="secs")); save_status()
    if(status$status=="incomplete_or_failed_no_final_inference") stop("One or more joint models failed fitting, extraction, or diagnostics. Check association_status.rds and the report.",call.=FALSE)
    message("Run saved: ",run_dir,". Review trace plots before interpretation.")
    invisible(run_dir)
  },error=function(e) {
    if (!is.null(e$audit)) jm_write_csv(e$audit, file.path(run_dir,paste0(e$audit_name,".csv")))
    status$status <- if(inherits(e,"jm_date_policy_error")) "blocked_date_policy" else "failed_no_final_inference"
    status$error <- conditionMessage(e)
    status$completed <- Sys.time()
    # This handler has its own frame; save the updated local status directly.
    jm_save(status,file.path(run_dir,"run_status.rds"))
    writeLines(conditionMessage(e),file.path(run_dir,"failure.txt"))
    stop(conditionMessage(e),"\nSaved diagnostics: ",run_dir,call.=FALSE)
  })
}
