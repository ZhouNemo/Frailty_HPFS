# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  8.1_joint_model_overall.R
# Author:  Nemo Zhou
# Date started:      2026-06-30
# Date last updated: 2026-07-17 (visual-data outputs moved; PNG writes removed)
#
# Purpose:
#   Fits the joint longitudinal-survival model for the OVERALL incident cancer
#   cohort to evaluate informative mortality in the frailty-trajectory analysis,
#   following Documents/Methods/Joint_Modeling_Informative_Mortality.md.
#
#   Key features required by the revised specification:
#     1. DEDUPLICATION to one survival contribution per physical participant.
#        Risk-set matching reuses controls across matched sets, so the raw
#        trajectory_id file would copy a participant's single death across every
#        assignment and count it as several independent events. We keep exactly
#        one assignment per participant: the participant's own case index if it
#        ever becomes a case, otherwise its earliest control assignment (ties
#        broken by smallest match_set). The random-effects / survival grouping
#        variable is therefore the participant id, not trajectory_id.
#     2. A spline basis that is COMMENSURABLE with the primary GLME: the natural
#        cubic-spline knots/boundaries are computed on the FULL +/-16 matched
#        cohort (the GLME recipe) before deduplication and then held fixed, so
#        the joint-model and GLME difference curves share a basis.
#     3. CO-PRIMARY association structures: current value and current value +
#        current slope (value-only is the fallback if value+slope fails).
#     4. A stratified PILOT profile (all cases + a seeded 25% control subsample)
#        for feasibility, and full reporting MCMC settings for inference.
#
#   Conceptual model:
#     - Longitudinal submodel (nlme::lme, REML):
#         fi_score_nocancer ~ Group * ns(Age_Centered; fixed GLME knots)
#                           + index_age_z + base_race + base_marital + base_pckgr
#                           + (1 + Age_Centered | id)
#     - Survival submodel (Cox, one row per participant):
#         time from assigned index to all-cause death, censored at the fixed
#         2020 analytic-cycle date.
#     - Joint model (JMbayes2::jm): mortality hazard linked to the latent current
#         value (and, co-primary, current slope) of no-cancer frailty.
#
# Run profile:
#   Sys.getenv("JM_PROFILE") in {"full" (default), "pilot"}. "pilot" fits the
#   joint model on the stratified subsample with quick MCMC and is labelled
#   non-inferential; the component submodels and descriptive outputs are always
#   produced on the full deduplicated cohort.
#
# Inputs:
#   Data/riskset_matched_overall_long.rds
#   Data/FI_longitudinal_1986_2020_IMPUTED_Cancer.rds
#   Results/cancer/data/7.2_ipcw_glme_spline_model.rds   (optional: knots cross-check)
#
# Outputs (Results/cancer/data, prefix 8.1_joint_model_overall_*):
#   *_longitudinal_dataset.csv        *_survival_dataset.csv
#   *_survival_summary.csv            *_dedup_audit.csv
#   *_cohort_shift.csv                *_lme_summary.txt
#   *_cox_summary.txt                 *_cox_zph.txt
#   *_predicted_trajectories.csv      *_fit.rds
#   *_association_estimates.csv       *_theta.csv
#   *_jm_summary.txt                  *_jm_diagnostic.txt
#   *_fit.rds
#
# Notes:
#   JMbayes2 is required for the linked fit. The component datasets and submodels
#   are always written before jm() is attempted, so a Bayesian-fit failure never
#   costs the prepared components.
# =============================================================================

library(dplyr)
library(ggplot2)
library(nlme)
library(splines)
library(survival)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir    <- file.path(project_dir, "Data")
results_dir <- file.path(project_dir, "Results", "cancer", "data")

matched_path   <- file.path(data_dir, "riskset_matched_overall_long.rds")
panel_path     <- file.path(data_dir, "FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
glme_model_path <- file.path(results_dir, "7.2_ipcw_glme_spline_model.rds")
out_prefix     <- "8.1_joint_model_overall"

window_yrs <- 16      # shared +/-16 analytic window with the primary GLME
spline_df  <- 4       # natural-spline df matching the primary GLME
jm_seed    <- 20260703
pilot_control_frac <- 0.25

# Pre-specified pre/post windows for the slope-contrast theta (matches the GLME
# rate contrast): post-index (0, +8] vs pre-index [-8, 0).
theta_pre  <- -8
theta_post <-  8

run_profile <- tolower(Sys.getenv("JM_PROFILE", "full"))
if (!run_profile %in% c("full", "pilot")) run_profile <- "full"

# Reporting vs pilot MCMC settings (Section 6 of the method file).
mcmc_full  <- list(n_chains = 3L, n_iter = 15000L, n_burnin = 5000L, n_thin = 5L)
mcmc_pilot <- list(n_chains = 1L, n_iter = 1500L,  n_burnin = 500L,  n_thin = 2L)

# Same fixed 2020 censoring convention used by the IPCW utilities: nominal
# mid-year date in months since 1900 for the 2020 analytic cycle.
admin_censor_date <- (2020 - 1900) * 12 + 6

needed_matched_cols <- c(
  "Cohort", "Group", "Age_Centered", "index_age_z",
  "base_race", "base_marital", "base_pckgr", "id", "match_set", "role",
  "index_date", "fi_score_nocancer"
)
needed_panel_cols <- c("id", "dtdth")

if (!file.exists(matched_path)) {
  stop(
    "Overall matched dataset not found at ", matched_path,
    ". Run Code/2_data_analysis/2.5_create_riskset_overall.R first."
  )
}
if (!file.exists(panel_path)) {
  stop("Longitudinal panel not found at ", panel_path, ".")
}
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

matched_long <- readRDS(matched_path)
missing_matched <- setdiff(needed_matched_cols, names(matched_long))
if (length(missing_matched) > 0) {
  stop("Matched dataset is missing required columns: ",
       paste(missing_matched, collapse = ", "))
}

panel <- readRDS(panel_path)
missing_panel <- setdiff(needed_panel_cols, names(panel))
if (length(missing_panel) > 0) {
  stop("Longitudinal panel is missing required columns: ",
       paste(missing_panel, collapse = ", "))
}

# ---------------------------------------------------------------------------
# 1. Full +/-16 matched analytic frame (used for knots + cohort-shift baseline)
# ---------------------------------------------------------------------------
analysis_long_full <- matched_long %>%
  mutate(
    trajectory_id = paste(Cohort, match_set, id, role, sep = "__"),
    Cohort = droplevels(factor(Cohort)),
    Group = factor(Group, levels = c("Control", "Cancer Case")),
    id = as.character(id),
    match_set = as.character(match_set),
    role = as.character(role),
    trajectory_id = as.character(trajectory_id),
    index_date = as.numeric(as.character(index_date)),
    Age_Centered = as.numeric(as.character(Age_Centered)),
    index_age_z = as.numeric(as.character(index_age_z)),
    fi_score_nocancer = as.numeric(as.character(fi_score_nocancer)),
    base_race = factor(base_race),
    base_marital = factor(base_marital),
    base_pckgr = factor(base_pckgr)
  ) %>%
  filter(
    Cohort == "All Cancer Cohort",
    abs(Age_Centered) <= window_yrs,
    !is.na(fi_score_nocancer),
    !is.na(index_date),
    !is.na(index_age_z),
    !is.na(base_race),
    !is.na(base_marital),
    !is.na(base_pckgr)
  ) %>%
  droplevels()

if (nrow(analysis_long_full) == 0) {
  stop("No longitudinal rows remain after applying the joint-model filters.")
}
if (length(unique(analysis_long_full$Group)) < 2) {
  stop("Both Control and Cancer Case rows are required for the joint model.")
}

# ---------------------------------------------------------------------------
# 2. GLME-commensurable spline basis, fixed on the FULL cohort before dedup
# ---------------------------------------------------------------------------
boundary_knots <- range(analysis_long_full$Age_Centered, na.rm = TRUE)
inner_knots <- as.numeric(quantile(
  analysis_long_full$Age_Centered,
  probs = seq(0, 1, length.out = spline_df + 1)[-c(1, spline_df + 1)],
  na.rm = TRUE
))

# Cross-check against the stored primary-GLME knots when available (informational).
glme_knots_note <- "primary-GLME knots artifact not found; reconstructed on full cohort"
if (file.exists(glme_model_path)) {
  glme_obj <- tryCatch(readRDS(glme_model_path), error = function(e) NULL)
  if (!is.null(glme_obj) && !is.null(glme_obj$spline_knots)) {
    glme_knots_note <- paste0(
      "stored GLME (7.2) inner knots [",
      paste(round(glme_obj$spline_knots, 4), collapse = ", "),
      "], boundary [", paste(round(glme_obj$boundary_knots, 4), collapse = ", "),
      "]; reconstructed inner [", paste(round(inner_knots, 4), collapse = ", "),
      "], boundary [", paste(round(boundary_knots, 4), collapse = ", "), "]"
    )
    rm(glme_obj)
  }
}

spline_basis <- ns(
  analysis_long_full$Age_Centered,
  knots = inner_knots,
  Boundary.knots = boundary_knots
)
spline_terms <- paste0("S", seq_len(ncol(spline_basis)))

fixed_formula <- as.formula(paste0(
  "fi_score_nocancer ~ Group * (", paste(spline_terms, collapse = " + "), ") + ",
  "index_age_z + base_race + base_marital + base_pckgr"
))

# ---------------------------------------------------------------------------
# 3. Deduplicate to one contribution per participant
# ---------------------------------------------------------------------------
assignments <- analysis_long_full %>%
  distinct(trajectory_id, id, Group, role, index_date, match_set) %>%
  mutate(
    is_case = as.integer(Group == "Cancer Case"),
    match_set_num = suppressWarnings(as.numeric(match_set))
  )

retained <- assignments %>%
  group_by(id) %>%
  arrange(desc(is_case), index_date, match_set_num, match_set, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup()

if (any(duplicated(retained$id))) {
  stop("Deduplication failed: more than one retained assignment for some id.")
}

analysis_long <- analysis_long_full %>%
  semi_join(retained %>% select(trajectory_id), by = "trajectory_id") %>%
  arrange(id, Age_Centered) %>%
  droplevels()

# Deduplication audit
n_full_assignments <- n_distinct(analysis_long_full$trajectory_id)
n_participants     <- n_distinct(analysis_long_full$id)
n_retained         <- nrow(retained)
dedup_audit <- tibble::tibble(
  metric = c("full_matched_assignments", "distinct_participants",
             "retained_contributions", "discarded_assignments",
             "retained_cancer_case", "retained_control"),
  value = c(n_full_assignments, n_participants, n_retained,
            n_full_assignments - n_retained,
            sum(retained$Group == "Cancer Case"),
            sum(retained$Group == "Control"))
)
write.csv(dedup_audit,
          file.path(results_dir, paste0(out_prefix, "_dedup_audit.csv")),
          row.names = FALSE)

# Cohort-shift: arm composition, matched (assignment-level) vs deduplicated
cohort_shift <- bind_rows(
  assignments %>%
    left_join(distinct(analysis_long_full, trajectory_id, index_age_z),
              by = "trajectory_id") %>%
    group_by(Group) %>%
    summarize(cohort = "matched_assignments", n = n(),
              mean_index_age_z = mean(index_age_z, na.rm = TRUE), .groups = "drop"),
  retained %>%
    left_join(distinct(analysis_long_full, trajectory_id, index_age_z),
              by = "trajectory_id") %>%
    group_by(Group) %>%
    summarize(cohort = "deduplicated", n = n(),
              mean_index_age_z = mean(index_age_z, na.rm = TRUE), .groups = "drop")
) %>%
  group_by(cohort) %>%
  mutate(arm_share = n / sum(n)) %>%
  ungroup() %>%
  select(cohort, Group, n, arm_share, mean_index_age_z)
write.csv(cohort_shift,
          file.path(results_dir, paste0(out_prefix, "_cohort_shift.csv")),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Survival dataset: one row per participant (structurally one death per id)
# ---------------------------------------------------------------------------
death_lookup <- panel %>%
  transmute(
    id = as.character(id),
    dtdth = as.numeric(as.character(dtdth))
  ) %>%
  group_by(id) %>%
  summarize(dtdth = suppressWarnings(max(dtdth, na.rm = TRUE)), .groups = "drop") %>%
  mutate(dtdth = if_else(is.infinite(dtdth), NA_real_, dtdth))

survival_data <- analysis_long %>%
  distinct(
    id, trajectory_id, Cohort, Group, match_set, role, index_date,
    index_age_z, base_race, base_marital, base_pckgr
  ) %>%
  left_join(death_lookup, by = "id") %>%
  mutate(
    death_event = as.integer(!is.na(dtdth) & dtdth > index_date & dtdth <= admin_censor_date),
    surv_end_date = if_else(death_event == 1L, dtdth, admin_censor_date),
    surv_time = (surv_end_date - index_date) / 12
  ) %>%
  filter(!is.na(surv_time), surv_time > 0) %>%
  droplevels()

# Keep the long rows in sync with retained survival participants.
analysis_long <- analysis_long %>%
  semi_join(select(survival_data, id), by = "id") %>%
  droplevels()
survival_data <- survival_data %>%
  semi_join(distinct(analysis_long, id), by = "id") %>%
  arrange(id) %>%
  droplevels()

# Correctness assertions: the death-duplication defect is gone iff there is
# exactly one survival row per participant.
if (nrow(survival_data) == 0) {
  stop("No survival records remain after requiring positive time from index to 2020 censoring.")
}
if (any(duplicated(survival_data$id))) {
  stop("Internal check failed: a participant appears more than once in the survival data.")
}
if (nrow(survival_data) != n_distinct(analysis_long$id)) {
  stop("Internal check failed: survival rows do not equal the number of distinct participants.")
}
if (length(unique(survival_data$Group)) < 2) {
  stop("Both Control and Cancer Case records are required in the survival data.")
}
if (any(survival_data$death_event == 1L & survival_data$dtdth <= survival_data$index_date)) {
  stop("Internal check failed: at least one death event occurs on/before index.")
}
if (any(survival_data$death_event == 1L & survival_data$dtdth > admin_censor_date)) {
  stop("Internal check failed: at least one death event occurs after the 2020 censor date.")
}

cat("\nRun profile:", run_profile, "\n")
cat("Knot note:", glme_knots_note, "\n")
cat("Joint-model longitudinal rows:", nrow(analysis_long), "\n")
cat("Deduplicated participants (survival rows):", nrow(survival_data), "\n")
print(table(Group = analysis_long$Group))
print(table(Group = survival_data$Group, Death = survival_data$death_event))

survival_summary <- survival_data %>%
  group_by(Group) %>%
  summarize(
    participants = n(),
    deaths_by_2020 = sum(death_event),
    percent_deaths_by_2020 = 100 * mean(death_event),
    mean_surv_time = mean(surv_time),
    median_surv_time = median(surv_time),
    min_surv_time = min(surv_time),
    max_surv_time = max(surv_time),
    .groups = "drop"
  )

# Materialize the fixed spline basis on the deduplicated long data.
analysis_long <- bind_cols(
  analysis_long,
  setNames(
    as.data.frame(predict(spline_basis, analysis_long$Age_Centered)),
    spline_terms
  )
)

write.csv(analysis_long,
          file.path(results_dir, paste0(out_prefix, "_longitudinal_dataset.csv")),
          row.names = FALSE)
write.csv(survival_data,
          file.path(results_dir, paste0(out_prefix, "_survival_dataset.csv")),
          row.names = FALSE)
write.csv(survival_summary,
          file.path(results_dir, paste0(out_prefix, "_survival_summary.csv")),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Longitudinal submodel (nlme::lme) with documented random-effect fallback
# ---------------------------------------------------------------------------
lme_control <- nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, niterEM = 25)
fit_lme <- function(random_form) {
  nlme::lme(fixed = fixed_formula, random = random_form, data = analysis_long,
            method = "REML", na.action = na.omit, control = lme_control)
}

cat("\nFitting nlme longitudinal submodel ...\n")
re_structure <- "~ 1 + Age_Centered | id"
m_lme <- tryCatch(
  fit_lme(~ 1 + Age_Centered | id),
  error = function(e) {
    message("Random-slope lme failed (", conditionMessage(e),
            "); falling back to random intercept only.")
    re_structure <<- "~ 1 | id (fallback; slope association loses RE slope info)"
    fit_lme(~ 1 | id)
  }
)
capture.output(summary(m_lme),
               file = file.path(results_dir, paste0(out_prefix, "_lme_summary.txt")))

# ---------------------------------------------------------------------------
# 6. Survival submodel (Cox) + PH diagnostic + KM-by-arm figure
# ---------------------------------------------------------------------------
cat("\nFitting Cox survival submodel ...\n")
m_cox <- survival::coxph(
  survival::Surv(surv_time, death_event) ~
    Group + index_age_z + base_race + base_marital + base_pckgr,
  data = survival_data, x = TRUE, model = TRUE, ties = "efron"
)
capture.output(summary(m_cox),
               file = file.path(results_dir, paste0(out_prefix, "_cox_summary.txt")))

zph <- tryCatch(survival::cox.zph(m_cox), error = function(e) NULL)
capture.output(
  if (is.null(zph)) cat("cox.zph could not be computed.\n") else print(zph),
  file = file.path(results_dir, paste0(out_prefix, "_cox_zph.txt"))
)

# ---------------------------------------------------------------------------
# 7. Adjusted predicted trajectories from the longitudinal submodel
# ---------------------------------------------------------------------------
make_ref_grid <- function(age_values) {
  expand.grid(
    Age_Centered = age_values,
    Group = factor(c("Control", "Cancer Case"),
                   levels = c("Control", "Cancer Case")),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      index_age_z = 0,
      base_race = factor(levels(analysis_long$base_race)[1],
                         levels = levels(analysis_long$base_race)),
      base_marital = factor(levels(analysis_long$base_marital)[1],
                            levels = levels(analysis_long$base_marital)),
      base_pckgr = factor(levels(analysis_long$base_pckgr)[1],
                          levels = levels(analysis_long$base_pckgr))
    )
}

prediction_grid <- make_ref_grid(seq(-window_yrs, window_yrs, by = 0.1))
prediction_grid <- bind_cols(
  prediction_grid,
  setNames(as.data.frame(predict(spline_basis, prediction_grid$Age_Centered)),
           spline_terms)
)
prediction_grid$pred <- as.numeric(predict(m_lme, newdata = prediction_grid, level = 0))
prediction_grid$Cohort <- "All Cancer Cohort"

pred_out <- prediction_grid %>% select(Cohort, Group, Age_Centered, pred)
write.csv(pred_out,
          file.path(results_dir, paste0(out_prefix, "_predicted_trajectories.csv")),
          row.names = FALSE)

group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")
p_traj <- ggplot(pred_out, aes(x = Age_Centered, y = pred, color = Group)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  scale_color_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    title = "Frailty Trajectories Used in Overall Cancer Joint Model",
    subtitle = paste0("nlme natural-spline submodel (df = ", spline_df,
                      ", shared GLME basis); joint model links latent FI to death risk"),
    x = "Years relative to index", y = "Predicted frailty index (no-cancer FI)",
    color = NULL
  ) +
  scale_x_continuous(breaks = seq(-window_yrs, window_yrs, by = 4))
# ---------------------------------------------------------------------------
# 8. theta contrast vector (post (0,+8] vs pre [-8,0) slope of the difference)
# ---------------------------------------------------------------------------
# Cancer-minus-Control design difference at relative time t (covariates cancel).
theta_cvec <- tryCatch({
  fixef_names <- names(nlme::fixef(m_lme))
  design_diff <- function(t) {
    nd <- make_ref_grid(t)
    nd <- bind_cols(nd,
                    setNames(as.data.frame(predict(spline_basis, nd$Age_Centered)),
                             spline_terms))
    X <- model.matrix(delete.response(terms(fixed_formula)), data = nd)
    # rows: Control then Cancer Case (make_ref_grid order for a single t)
    as.numeric(X[nd$Group == "Cancer Case", ] - X[nd$Group == "Control", ])
  }
  cvec <- (design_diff(theta_post) - design_diff(0)) / theta_post -
          (design_diff(0) - design_diff(theta_pre)) / (-theta_pre)
  names(cvec) <- colnames(model.matrix(delete.response(terms(fixed_formula)),
                                       data = bind_cols(
                                         make_ref_grid(0),
                                         setNames(as.data.frame(predict(spline_basis, 0)),
                                                  spline_terms))))
  cvec[fixef_names]
}, error = function(e) {
  message("theta contrast vector could not be built: ", conditionMessage(e))
  NULL
})

theta_glme <- NULL
if (!is.null(theta_cvec)) {
  b <- nlme::fixef(m_lme)
  V <- as.matrix(vcov(m_lme))
  est <- as.numeric(theta_cvec %*% b)
  se  <- sqrt(as.numeric(t(theta_cvec) %*% V %*% theta_cvec))
  theta_glme <- c(estimate = est, se = se,
                  lcl = est - 1.96 * se, ucl = est + 1.96 * se)
}

# ---------------------------------------------------------------------------
# 9. Joint models (co-primary: value, and value + slope)
# ---------------------------------------------------------------------------
if (!requireNamespace("JMbayes2", quietly = TRUE)) {
  writeLines(
    c("JMbayes2 not installed; component datasets and submodels were written.",
      "Install JMbayes2 and rerun to fit the linked joint model."),
    con = file.path(results_dir, paste0(out_prefix, "_jm_diagnostic.txt"))
  )
  saveRDS(
    list(joint_model_value = NULL, joint_model_value_slope = NULL,
         longitudinal_model = m_lme, survival_model = m_cox,
         longitudinal_data = analysis_long, survival_data = survival_data,
         survival_summary = survival_summary, dedup_audit = dedup_audit,
         cohort_shift = cohort_shift, spline_knots = inner_knots,
         boundary_knots = boundary_knots, re_structure = re_structure,
         theta_glme = theta_glme, admin_censor_date = admin_censor_date,
         run_profile = run_profile),
    file.path(results_dir, paste0(out_prefix, "_fit.rds"))
  )
  stop("JMbayes2 is required to fit the joint model. Component outputs written to ", results_dir)
}

mcmc <- if (run_profile == "pilot") mcmc_pilot else mcmc_full

# Data and component submodels passed to jm(): full for reporting, stratified
# subsample for the pilot feasibility profile.
if (run_profile == "pilot") {
  set.seed(jm_seed)
  case_ids    <- survival_data$id[survival_data$Group == "Cancer Case"]
  control_ids <- survival_data$id[survival_data$Group == "Control"]
  keep_ids <- c(case_ids,
                sample(control_ids, size = ceiling(pilot_control_frac * length(control_ids))))
  jm_long <- analysis_long %>% filter(id %in% keep_ids) %>% droplevels()
  jm_surv <- survival_data %>% filter(id %in% keep_ids) %>% arrange(id) %>% droplevels()
  jm_lme  <- tryCatch(
    nlme::lme(fixed = fixed_formula, random = ~ 1 + Age_Centered | id,
              data = jm_long, method = "REML", na.action = na.omit, control = lme_control),
    error = function(e) nlme::lme(fixed = fixed_formula, random = ~ 1 | id,
                                  data = jm_long, method = "REML",
                                  na.action = na.omit, control = lme_control))
  jm_cox  <- survival::coxph(
    survival::Surv(surv_time, death_event) ~
      Group + index_age_z + base_race + base_marital + base_pckgr,
    data = jm_surv, x = TRUE, model = TRUE, ties = "efron")
} else {
  jm_long <- analysis_long; jm_surv <- survival_data
  jm_lme  <- m_lme;         jm_cox  <- m_cox
}

fit_jm <- function(ff) {
  JMbayes2::jm(
    jm_cox, jm_lme, time_var = "Age_Centered",
    data_Surv = jm_surv, id_var = "id",
    functional_forms = ff,
    n_chains = mcmc$n_chains, n_iter = mcmc$n_iter,
    n_burnin = mcmc$n_burnin, n_thin = mcmc$n_thin,
    cores = 1L, save_random_effects = FALSE, seed = jm_seed
  )
}

cat("\nFitting JMbayes2 joint model (value association) ...\n")
jm_err_value <- NULL
jm_value <- tryCatch(
  fit_jm(~ value(fi_score_nocancer)),
  error = function(e) { jm_err_value <<- conditionMessage(e); NULL }
)

cat("Fitting JMbayes2 joint model (value + slope association) ...\n")
jm_err_vs <- NULL
jm_value_slope <- tryCatch(
  fit_jm(~ value(fi_score_nocancer) + slope(fi_score_nocancer)),
  error = function(e) { jm_err_vs <<- conditionMessage(e); NULL }
)

jm_fits <- list(value = jm_value, value_slope = jm_value_slope)
jm_fits <- jm_fits[!vapply(jm_fits, is.null, logical(1))]

# --- Association estimates: HR per 0.1 FI (value) / per 0.1 FI/year (slope) ---
assoc_rows <- list()
for (nm in names(jm_fits)) {
  fit <- jm_fits[[nm]]
  draws <- tryCatch(do.call(rbind, fit$mcmc$alphas), error = function(e) NULL)
  if (is.null(draws) || ncol(draws) == 0) next
  hr01 <- exp(0.1 * draws)
  for (j in seq_len(ncol(draws))) {
    assoc_rows[[length(assoc_rows) + 1]] <- data.frame(
      structure = nm,
      association = colnames(draws)[j],
      alpha_mean = mean(draws[, j]),
      alpha_sd = sd(draws[, j]),
      hr_per_0.1_mean = mean(hr01[, j]),
      hr_per_0.1_lcl = as.numeric(quantile(hr01[, j], 0.025)),
      hr_per_0.1_ucl = as.numeric(quantile(hr01[, j], 0.975)),
      stringsAsFactors = FALSE
    )
  }
}
if (length(assoc_rows) > 0) {
  write.csv(do.call(rbind, assoc_rows),
            file.path(results_dir, paste0(out_prefix, "_association_estimates.csv")),
            row.names = FALSE)
}

# --- theta_JM (from longitudinal fixed-effect posterior) and attenuation ------
theta_out <- tryCatch({
  rows <- list()
  if (!is.null(theta_glme)) {
    rows[[length(rows) + 1]] <- data.frame(
      quantity = "theta_GLME_dedup", structure = NA_character_,
      estimate = unname(theta_glme["estimate"]),
      lcl = unname(theta_glme["lcl"]), ucl = unname(theta_glme["ucl"]),
      stringsAsFactors = FALSE)
  }
  if (!is.null(theta_cvec)) {
    for (nm in names(jm_fits)) {
      betas <- tryCatch(do.call(rbind, jm_fits[[nm]]$mcmc$betas1),
                        error = function(e) NULL)
      if (is.null(betas)) next
      cv <- theta_cvec
      common <- intersect(names(cv), colnames(betas))
      if (length(common) < length(cv)) next
      td <- as.numeric(betas[, names(cv), drop = FALSE] %*% cv)
      att <- if (!is.null(theta_glme)) mean(td) - unname(theta_glme["estimate"]) else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        quantity = "theta_JM", structure = nm,
        estimate = mean(td),
        lcl = as.numeric(quantile(td, 0.025)),
        ucl = as.numeric(quantile(td, 0.975)),
        attenuation_vs_glme = att,
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) > 0) dplyr::bind_rows(rows) else NULL
}, error = function(e) { message("theta computation failed: ", conditionMessage(e)); NULL })
if (!is.null(theta_out)) {
  write.csv(theta_out, file.path(results_dir, paste0(out_prefix, "_theta.csv")),
            row.names = FALSE)
}

# --- Joint-model summaries + convergence diagnostics --------------------------
jm_summary_con <- file.path(results_dir, paste0(out_prefix, "_jm_summary.txt"))
if (length(jm_fits) > 0) {
  capture.output(
    for (nm in names(jm_fits)) {
      cat("\n==================== association:", nm, "====================\n")
      print(summary(jm_fits[[nm]]))
    },
    file = jm_summary_con
  )
}

diag_lines <- c(
  paste0("Run profile: ", run_profile,
         if (run_profile == "pilot") " (NON-INFERENTIAL pilot on stratified subsample)" else ""),
  paste0("MCMC: chains=", mcmc$n_chains, " iter=", mcmc$n_iter,
         " burnin=", mcmc$n_burnin, " thin=", mcmc$n_thin, " seed=", jm_seed),
  paste0("Random-effects structure: ", re_structure),
  paste0("Spline knots: ", glme_knots_note),
  paste0("value association fit: ", if (is.null(jm_value)) paste0("FAILED - ", jm_err_value) else "OK"),
  paste0("value+slope association fit: ", if (is.null(jm_value_slope)) paste0("FAILED - ", jm_err_vs) else "OK"),
  "",
  "Convergence (target: R-hat < 1.05; ESS > 400 for association alpha):"
)
for (nm in names(jm_fits)) {
  st <- tryCatch(jm_fits[[nm]]$statistics, error = function(e) NULL)
  rhat <- tryCatch({
    all_rhat <- unlist(lapply(st$Rhat, function(x) x))
    if (is.null(all_rhat)) NA else max(all_rhat, na.rm = TRUE)
  }, error = function(e) NA)
  diag_lines <- c(diag_lines,
                  paste0("  [", nm, "] max R-hat across monitored params: ",
                         ifelse(is.na(rhat), "unavailable", round(rhat, 4))))
}
if (length(jm_fits) == 0) {
  diag_lines <- c(diag_lines,
                  "No joint model converged; component submodels and datasets are available as feasibility evidence.")
}
writeLines(diag_lines, con = file.path(results_dir, paste0(out_prefix, "_jm_diagnostic.txt")))

# ---------------------------------------------------------------------------
# 10. Persist everything
# ---------------------------------------------------------------------------
saveRDS(
  list(
    joint_model_value = jm_value,
    joint_model_value_slope = jm_value_slope,
    joint_model_error_value = jm_err_value,
    joint_model_error_value_slope = jm_err_vs,
    longitudinal_model = m_lme,
    survival_model = m_cox,
    longitudinal_data = analysis_long,
    survival_data = survival_data,
    survival_summary = survival_summary,
    dedup_audit = dedup_audit,
    cohort_shift = cohort_shift,
    theta_glme = theta_glme,
    theta_table = theta_out,
    spline_knots = inner_knots,
    boundary_knots = boundary_knots,
    re_structure = re_structure,
    run_profile = run_profile,
    admin_censor_date = admin_censor_date
  ),
  file.path(results_dir, paste0(out_prefix, "_fit.rds"))
)

cat("\nSaved overall joint-model outputs with prefix '", out_prefix,
    "' to: ", results_dir, "\n", sep = "")
