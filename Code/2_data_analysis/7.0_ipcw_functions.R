# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  7.0_ipcw_functions.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-20 (canonical earliest-cancer index date)
#
# Purpose:
#   Shared inverse-probability-of-censoring-weighting (IPCW) utilities used by
#   7.1+. Builds STABILIZED weights for the observation (questionnaire-response)
#   process so the weighted trajectory models (4.x event-study, 6.x quadratic)
#   reduce bias from informative loss to follow-up and questionnaire nonresponse
#   after the index.
#
#   IMPORTANT scope of the weights:
#     * The weight model is for REMAINING OBSERVED (alive AND responding), modeled
#       only among person-cycles still at risk (alive and observed at the prior
#       analytic cycle). DEATH is NOT weighted away: a dead person has no frailty
#       to recover, so death is left as a truncating event handled by the estimand
#       (per project guidance), not as recoverable censoring. These weights target
#       loss to follow-up / nonresponse among survivors.
#
#   Censoring-model predictors (all measured at the PRIOR observed cycle, i.e.
#   lagged, so we never condition on the unobserved current questionnaire):
#     * prior frailty: most recent observed fi_score_nocancer (p_fi) + baseline FI
#     * cancer status/timing: time-varying cancer_now + years since diagnosis (tsc)
#     * aging / calendar: natural spline in attained age + questionnaire cycle
#     * response history: cumulative count of prior observed cycles
#     * Section-10 covariates (lagged): race, marital_status, living_arr, smoke,
#       pckyr, bmi, act, alco
#   Time-varying comorbidity / disease items are deliberately EXCLUDED (prior FI
#   already summarizes health); marital_status and living_arr use the raw
#   time-varying columns, not the baseline factors.
#
#   Stabilized weight = cumulative product over at-risk cycles of
#     P(observed | numerator model) / P(observed | denominator model),
#   numerator = baseline + exposure + age + cycle only; denominator adds the
#   lagged time-varying predictors and response history. Weights are then made
#   index-conditional (post-index only) and truncated.
#
# Units:
#   Dates are months since 1900; age = (date - dbmy09)/12. Nominal cycle dates are
#   used so non-responding person-cycles still have an age and alive status.
# =============================================================================

library(dplyr)
library(tidyr)
library(splines)

# Default censoring-model covariates, all from section 10 of
# Documents/Data Dictionary/Analysis_Required_Variables.md.
IPCW_CAT_COVS <- c("race", "marital_status", "living_arr", "smoke")
IPCW_NUM_COVS <- c("bmi", "act", "alco", "pckyr")

# Nominal mid-year date (months since 1900) for a 2-digit cycle code.
ipcw_cycle_to_date <- function(cycle) {
  yr2  <- as.integer(cycle)
  year <- ifelse(yr2 >= 50, 1900L + yr2, 2000L + yr2)
  (year - 1900) * 12 + 6
}

# ---- 1) Build stabilized id-cycle IPCW from the longitudinal panel ----------
build_ipcw_weights <- function(panel_path,
                               target_cycles,
                               cat_covs = IPCW_CAT_COVS,
                               num_covs = IPCW_NUM_COVS,
                               age_df = 4) {

  if (!file.exists(panel_path)) stop("Panel not found at ", panel_path)
  panel <- readRDS(panel_path)
  if (!"cancer_index_dateca" %in% names(panel)) {
    stop(
      "Panel is missing cancer_index_dateca; run 7.4_cancer_subtypes.R first.",
      call. = FALSE
    )
  }

  present_cat <- intersect(cat_covs, names(panel))
  present_num <- intersect(num_covs, names(panel))
  miss <- setdiff(c(cat_covs, num_covs), c(present_cat, present_num))
  if (length(miss) > 0) message("IPCW: covariates not in panel, skipped: ",
                                paste(miss, collapse = ", "))

  pan <- panel %>%
    transmute(
      id    = as.character(id),
      cycle = as.character(cycle),
      participated = as.numeric(as.character(participated)),
      fi    = as.numeric(as.character(fi_score_nocancer)),
      dbmy09 = as.numeric(as.character(dbmy09)),
      dtdth  = as.numeric(as.character(dtdth)),
      cancer_dateca = as.numeric(as.character(cancer_index_dateca)),
      across(all_of(present_cat), ~ as.character(.x)),
      across(all_of(present_num), ~ as.numeric(as.character(.x)))
    ) %>%
    filter(cycle %in% target_cycles)

  person <- pan %>%
    group_by(id) %>%
    summarize(
      dbmy09        = dbmy09[which(!is.na(dbmy09))[1]],
      dtdth         = suppressWarnings(max(dtdth, na.rm = TRUE)),
      cancer_dateca = suppressWarnings(min(cancer_dateca, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      dtdth = ifelse(is.infinite(dtdth), NA_real_, dtdth),
      cancer_dateca = ifelse(is.infinite(cancer_dateca), NA_real_, cancer_dateca)
    ) %>%
    filter(!is.na(dbmy09))

  # one record per observed id-cycle with covariate snapshot
  obs_cov <- pan %>%
    group_by(id, cycle) %>%
    summarize(
      observed_raw = as.integer(any(participated == 1 & !is.na(fi))),
      fi_obs = fi[which(participated == 1 & !is.na(fi))[1]],
      across(all_of(present_cat), ~ .x[1]),
      across(all_of(present_num), ~ .x[1]),
      .groups = "drop"
    )

  cyc_tab <- tibble(
    cycle = target_cycles,
    cycle_order = seq_along(target_cycles),
    cycle_date = ipcw_cycle_to_date(target_cycles)
  )

  skel <- tidyr::expand_grid(id = unique(person$id), cycle = target_cycles) %>%
    left_join(cyc_tab,  by = "cycle") %>%
    left_join(person,   by = "id") %>%
    left_join(obs_cov,  by = c("id", "cycle")) %>%
    mutate(
      observed   = ifelse(is.na(observed_raw), 0L, observed_raw),
      age_nom    = (cycle_date - dbmy09) / 12,
      alive      = as.integer(is.na(dtdth) | dtdth >= cycle_date),
      cancer_now = as.integer(!is.na(cancer_dateca) & cancer_dateca <= cycle_date),
      tsc        = ifelse(!is.na(cancer_dateca) & cancer_dateca <= cycle_date,
                          (cycle_date - cancer_dateca) / 12, 0)
    ) %>%
    arrange(id, cycle_order)

  # observed-only covariate snapshots, carried forward, then lagged one cycle
  skel <- skel %>%
    mutate(
      across(all_of(present_cat), ~ if_else(observed == 1, as.character(.x), NA_character_), .names = "{.col}__m"),
      across(all_of(present_num), ~ if_else(observed == 1, as.numeric(.x),   NA_real_),      .names = "{.col}__m"),
      fi__m = if_else(observed == 1, fi_obs, NA_real_)
    ) %>%
    group_by(id) %>%
    tidyr::fill(ends_with("__m"), .direction = "down") %>%
    mutate(
      across(ends_with("__m"), ~ dplyr::lag(.x), .names = "prior_{.col}"),
      entered_prior      = dplyr::lag(cumsum(observed) > 0),
      prior_observed     = dplyr::lag(observed),
      cum_observed_prior = dplyr::lag(cumsum(observed)),
      baseline_fi        = { v <- na.omit(fi__m); if (length(v)) v[1] else NA_real_ }
    ) %>%
    ungroup() %>%
    mutate(
      entered_prior      = ifelse(is.na(entered_prior), FALSE, entered_prior),
      prior_observed     = ifelse(is.na(prior_observed), 0L, prior_observed),
      cum_observed_prior = ifelse(is.na(cum_observed_prior), 0L, cum_observed_prior),
      at_risk = as.integer(alive == 1 & entered_prior & prior_observed == 1 & !is.na(age_nom))
    )

  # tidy predictor names: prior_<cov>__m -> p_<cov>; prior_fi__m -> p_fi
  skel <- skel %>%
    rename_with(~ gsub("__m$", "", gsub("^prior_", "p_", .x)),
                .cols = starts_with("prior_") & ends_with("__m"))

  p_cat <- paste0("p_", present_cat)
  p_num <- paste0("p_", present_num)

  # ---- model frame on at-risk person-cycles; drop degenerate predictors ----
  atr <- skel %>% filter(at_risk == 1)
  num_med <- sapply(c(p_num, "p_fi", "baseline_fi"),
                    function(v) stats::median(atr[[v]], na.rm = TRUE))

  build_md <- function(d) {
    d <- d %>%
      mutate(
        cyc_f = droplevels(factor(cycle, levels = target_cycles)),
        across(all_of(p_cat),
               ~ droplevels(factor(if_else(is.na(.x), "Missing", as.character(.x)))))
      )
    for (v in c(p_num, "p_fi", "baseline_fi")) {
      d[[v]] <- ifelse(is.na(d[[v]]), num_med[[v]], d[[v]])
    }
    d
  }
  atr <- build_md(atr)

  # keep only predictors that actually vary among at-risk rows; a factor with one
  # level (e.g. a covariate that is all "Missing" after lagging, or the first
  # analytic cycle which is never at risk) breaks glm contrasts.
  fac_all  <- c("cyc_f", p_cat)
  fac_keep <- fac_all[vapply(fac_all, function(v) nlevels(atr[[v]]) >= 2, logical(1))]
  num_all  <- c("cancer_now", "tsc", "p_fi", "baseline_fi", "cum_observed_prior", p_num)
  num_keep <- num_all[vapply(num_all, function(v) {
    x <- atr[[v]]; is.numeric(x) && !all(is.na(x)) && stats::sd(x, na.rm = TRUE) > 0
  }, logical(1))]
  dropped <- setdiff(c(fac_all, num_all), c(fac_keep, num_keep))
  if (length(dropped) > 0)
    message("IPCW: dropped degenerate predictors (no variation among at-risk): ",
            paste(dropped, collapse = ", "))

  age_term  <- sprintf("ns(age_nom, df = %d)", age_df)
  den_terms <- c(age_term, fac_keep, num_keep)
  num_terms <- c(age_term,
                 if ("cyc_f" %in% fac_keep) "cyc_f",
                 intersect(c("cancer_now", "tsc", "baseline_fi"), num_keep))

  den_form <- as.formula(paste("observed ~", paste(den_terms, collapse = " + ")))
  num_form <- as.formula(paste("observed ~", paste(num_terms, collapse = " + ")))

  cat("\nFitting IPCW denominator model on", nrow(atr), "at-risk person-cycles ...\n")
  m_den <- glm(den_form, data = atr, family = binomial())
  m_num <- glm(num_form, data = atr, family = binomial())

  atr$pden <- as.numeric(predict(m_den, atr, type = "response"))
  atr$pnum <- as.numeric(predict(m_num, atr, type = "response"))

  # join predictions back, form per-cycle ratio (1 when not at risk), cumprod
  skel <- skel %>%
    left_join(select(atr, id, cycle, pden, pnum), by = c("id", "cycle")) %>%
    mutate(
      r_t = ifelse(at_risk == 1 & !is.na(pden) & pden > 0, pnum / pden, 1)
    ) %>%
    group_by(id) %>%
    arrange(cycle_order, .by_group = TRUE) %>%
    mutate(sw_cum = cumprod(ifelse(is.na(r_t), 1, r_t))) %>%
    ungroup()

  idcycle <- skel %>%
    select(id, cycle, cycle_order, observed, alive, at_risk,
           age_nom, cancer_now, tsc, pden, pnum, r_t, sw_cum)

  cat("Stabilized id-cycle weight summary (should center near 1):\n")
  print(summary(idcycle$sw_cum[idcycle$observed == 1]))

  invisible(list(idcycle = idcycle, model_den = m_den, model_num = m_num,
                 covariates = list(cat = present_cat, num = present_num)))
}

# ---- 2) Attach post-index conditional IPCW to a matched dataset -------------
attach_ipcw_to_matched <- function(matched_path,
                                   idcycle,
                                   out_path,
                                   results_dir = NULL,
                                   out_prefix = NULL,
                                   trunc = c(0.01, 0.99)) {

  if (!file.exists(matched_path)) stop("Matched dataset not found at ", matched_path)
  ml <- readRDS(matched_path) %>%
    mutate(id = as.character(id), cycle = as.character(cycle))

  w <- idcycle %>% distinct(id, cycle, sw_cum)

  ml <- ml %>%
    left_join(w, by = c("id", "cycle")) %>%
    mutate(
      sw_cum = ifelse(is.na(sw_cum), 1, sw_cum),
      trajectory_id = paste(Cohort, match_set, id, role, sep = "__")
    )

  # index reference = cumulative weight at the last pre-index observed cycle
  idx_ref <- ml %>%
    filter(Age_Centered < 0) %>%
    group_by(trajectory_id) %>%
    arrange(Age_Centered, .by_group = TRUE) %>%
    summarize(sw_index_ref = dplyr::last(sw_cum), .groups = "drop")

  ml <- ml %>%
    left_join(idx_ref, by = "trajectory_id") %>%
    mutate(
      sw_index_ref = ifelse(is.na(sw_index_ref), 1, sw_index_ref),
      sw_ipcw_raw  = ifelse(Age_Centered >= 0, sw_cum / sw_index_ref, 1)
    )

  qs <- quantile(ml$sw_ipcw_raw, probs = trunc, na.rm = TRUE)
  ml <- ml %>% mutate(sw_ipcw = pmin(pmax(sw_ipcw_raw, qs[1]), qs[2]))

  saveRDS(ml, out_path)

  cat("\nAttached IPCW to:", basename(matched_path), "\n")
  cat("Truncation bounds [", round(qs[1], 3), ",", round(qs[2], 3), "]\n")
  cat("Post-index sw_ipcw summary (post-index rows):\n")
  print(summary(ml$sw_ipcw[ml$Age_Centered >= 0]))
  cat("Saved IPCW-augmented matched dataset to:", out_path, "\n")

  if (!is.null(results_dir) && !is.null(out_prefix)) {
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
    write.csv(
      ml %>% summarize(
        n_rows = dplyr::n(),
        mean_w = mean(sw_ipcw, na.rm = TRUE),
        sd_w = sd(sw_ipcw, na.rm = TRUE),
        min_w = min(sw_ipcw, na.rm = TRUE),
        max_w = max(sw_ipcw, na.rm = TRUE)
      ),
      file.path(results_dir, paste0(out_prefix, "_ipcw_weight_summary.csv")),
      row.names = FALSE
    )
  }

  invisible(ml)
}

# ---- 3) Convenience: build + attach in one call -----------------------------
run_ipcw <- function(panel_path, matched_path, out_path,
                     target_cycles, results_dir = NULL, out_prefix = NULL,
                     cat_covs = IPCW_CAT_COVS, num_covs = IPCW_NUM_COVS,
                     age_df = 4, trunc = c(0.01, 0.99)) {
  w <- build_ipcw_weights(panel_path, target_cycles, cat_covs, num_covs, age_df)
  attach_ipcw_to_matched(matched_path, w$idcycle, out_path,
                         results_dir = results_dir, out_prefix = out_prefix,
                         trunc = trunc)
}
