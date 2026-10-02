# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-20
# Purpose: Recompute CR2/CR0 covariance and large-sample inference from accepted,
# cached overall M1 fits. Never invokes lmer, matching, Satterthwaite or HTZ.
# Validates source hashes, fixed-design reconstruction and saved theta before
# inference. Quantifies both row/reuse concentration and theta-specific sandwich
# contributions. The latter are descriptive and are not Satterthwaite df.
# Output: Codex/glme_diagnostics/<cached-run>/asymptotic_robust (diagnostics only).
# Usage from project root:
#   Rscript Code/tests/4.0_diagnose_cached_inference.R
# Optional first argument: path to a different overall diagnostic run directory.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(lme4); library(clubSandwich) })
source('Code/2_data_analysis/4.0_GLME_spline_functions.R')
# Exact nested-participant working covariance, avoiding repeated wide sparse
# row slices in clubSandwich 0.6.2. Retain its covariance/adjustment algorithms.
# Reject crossed groups/weights; verify blocks against lme4's sparse Z Lambda.
fast_nested_engine <- function(fit, cluster) {
  stopifnot(as.character(packageVersion('clubSandwich')) == '0.6.2',
            is.null(fit@call$weights))
  cluster <- droplevels(factor(cluster))
  groups <- getME(fit, 'flist')
  stopifnot(all(vapply(groups, function(g) identical(as.character(g), as.character(cluster)), logical(1))))
  rows <- split(seq_along(cluster), cluster)
  frame <- model.frame(fit)
  vc <- VarCorr(fit, sigma = 1)
  designs <- lapply(vc, function(G) {
    terms <- colnames(G)
    Z <- matrix(1, nrow(frame), length(terms), dimnames = list(NULL, terms))
    for (v in setdiff(terms, '(Intercept)')) Z[, v] <- frame[[v]]
    Z
  })
  target <- lapply(rows, function(i) {
    T <- diag(length(i))
    for (j in seq_along(vc)) {
      Z <- designs[[j]][i, , drop = FALSE]
      T <- T + Z %*% as.matrix(vc[[j]]) %*% t(Z)
    }
    T
  })
  sparse <- Matrix::tcrossprod(getME(fit, 'Z'), getME(fit, 'Lambdat'))
  probe <- unique(c(1L, length(rows), which.max(lengths(rows)),
                    round(seq(1, length(rows), length.out = 10))))
  for (j in probe) {
    Z <- sparse[rows[[j]], , drop = FALSE]
    original <- as.matrix(Matrix::tcrossprod(Z)) + diag(length(rows[[j]]))
    stopifnot(isTRUE(all.equal(unname(target[[j]]), unname(original), tolerance = 1e-10)))
  }
  rm(sparse)
  weights <- lapply(target, function(T) chol2inv(chol(T)))
  core <- getFromNamespace('vcov_CR', 'clubSandwich')
  scope <- new.env(parent = environment(core))
  scope$weightMatrix <- function(obj, cluster) weights
  scope$targetVariance <- function(obj, cluster) target
  # Reuse row indices instead of scanning all observations once per cluster.
  scope$matrix_list <- function(x, fac, dim) {
    stopifnot(identical(fac, cluster), dim == 'row', is.matrix(x))
    lapply(rows, function(i) x[i, , drop = FALSE])
  }
  environment(core) <- scope
  function(type, form = 'sandwich') core(fit, cluster, type = type,
                                         inverse_var = TRUE, form = form)
}

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else paste0('Codex/glme_diagnostics/',
  'riskset_matched_overall_long/364e02ea_0cf6c5ee_f1fa5b8b')
root <- normalizePath(root, mustWork = TRUE)
cohort_dir <- file.path(root, 'All_Cancer_Cohort')
out <- file.path(root, 'asymptotic_robust')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
provenance <- readRDS(file.path(root, 'provenance.rds'))
current <- validate_matching_provenance('Data/riskset_matched_overall_long.rds')
stopifnot(identical(current$input_md5, provenance$input$input_md5))
script_hash <- unname(tools::md5sum('Code/tests/4.0_diagnose_cached_inference.R'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))

# Reconstruct only the fixed-effect spline coordinates, never a fitted model.
d <- readRDS('Data/riskset_matched_overall_long.rds')
if (!'trajectory_id' %in% names(d)) d$trajectory_id <- paste(d$Cohort, d$match_set, d$id, d$role, sep = '__')
d <- d %>% filter(Cohort == 'All Cancer Cohort', !post_own_cancer,
  !is.na(fi_score_nocancer), !is.na(Age_Centered), if_all(all_of(.primary_covars), ~ !is.na(.x))) %>%
  mutate(Group = factor(Group, levels = c('Control', 'Cancer Case')),
         across(all_of(c('id', 'base_race', 'base_marital', 'base_living')), factor)) %>% droplevels()
spec <- .trajectory_model_specs$M1_primary_spline
bases <- list(adjusted_spline = .make_spline_basis(d$Age_Centered, 3, 'S'))
reference <- .make_ref_grid(d, covars = spec$covars, .reference = 1)
scaled <- .scale_model_numeric_columns(.add_model_time_terms(d, spec, bases), bases$adjusted_spline$terms)
expected_X <- model.matrix(.fixed_rhs(bases$adjusted_spline$terms, spec$covars), scaled$data)
reuse <- d %>% group_by(id) %>% summarize(n_rows = n(), n_cycles = n_distinct(cycle),
  n_assignments = n_distinct(trajectory_id), .groups = 'drop')
write.csv(reuse, file.path(out, 'participant_clusters.csv'), row.names = FALSE)
cluster_summary <- bind_rows(lapply(c('n_rows', 'n_cycles', 'n_assignments'), function(v) {
  x <- reuse[[v]]; shares <- sort(x / sum(x), decreasing = TRUE)
  data.frame(quantity = v, n_clusters = length(x), min = min(x), median = median(x),
    p90 = unname(quantile(x, .9)), p95 = unname(quantile(x, .95)),
    p99 = unname(quantile(x, .99)), max = max(x), max_share = shares[1],
    top_1pct_share = sum(head(shares, ceiling(.01 * length(x)))),
    inverse_share_hhi = 1 / sum(shares^2))
}))
write.csv(cluster_summary, file.path(out, 'cluster_summary.csv'), row.names = FALSE)
previous <- read.csv(file.path(cohort_dir, 'fit_comparison.csv'))
results <- list(); coefficients <- list(); concentration <- list(); audit <- list()
models <- c('M1_intercept', 'M1_uncorrelated_scaled', 'M1_attained_age')
for (name in models) {
  message('Reading cached fit: ', name)
  fit_path <- file.path(cohort_dir, paste0(name, '_fit.rds'))
  fit_md5 <- unname(tools::md5sum(fit_path))
  fit <- readRDS(fit_path)$fit
  frame <- model.frame(fit); beta <- fixef(fit); fe <- names(beta)
  stopifnot(identical(dim(getME(fit, 'X')), dim(expected_X[, fe, drop = FALSE])),
    identical(colnames(getME(fit, 'X')), fe),
    identical(as.character(frame$id), as.character(d$id)),
    isTRUE(all.equal(as.numeric(model.response(frame)), d$fi_score_nocancer, tolerance = 1e-12)),
    isTRUE(all.equal(unname(getME(fit, 'X')), unname(expected_X[, fe, drop = FALSE]), tolerance = 1e-10, check.attributes = FALSE)))
  grids <- .spline_grid_factory(reference, spec, bases, scaled$parameters, fe)
  C <- .exact_slope_contrasts(grids$difference, 8)$theta
  theta <- as.numeric(C %*% beta)
  stopifnot(abs(theta - previous$theta[previous$configuration == name]) < 1e-10)
  terms <- .spline_interaction_terms(fe, bases$adjusted_spline$terms)
  A <- diag(length(beta))[match(terms, fe), , drop = FALSE]
  stopifnot(nrow(A) == 3L)
  fast_cov <- fast_nested_engine(fit, frame$id)
  cache_key <- list(fit_md5 = fit_md5, input_md5 = current$input_md5,
    script_md5 = script_hash, clubSandwich = as.character(packageVersion('clubSandwich')))
  for (type in c('CR2', 'CR0')) {
    message('  Computing ', type, ' covariance (no small-sample tests)')
    cache_path <- file.path(out, paste0(name, '_', type, '.rds'))
    cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
    if (!is.null(cached) && identical(cached$key, cache_key)) {
      V <- cached$V
    } else {
      cv <- fast_cov(type)
      V <- as.matrix(cv)
      .validate_covariance(V, fe)
      if (name == 'M1_intercept' && type == 'CR2') {
        baseline_path <- file.path(out, 'public_CR2_reference.rds')
        if (!file.exists(baseline_path)) {
          message('  Computing one public-package reference (may take several minutes)')
          reference_cov <- clubSandwich::vcovCR(fit, cluster = frame$id, type = type)
          saveRDS(list(key = cache_key, V = as.matrix(reference_cov)), baseline_path)
          rm(reference_cov); gc()
        }
        baseline <- readRDS(baseline_path)
        stopifnot(identical(baseline$key$fit_md5, fit_md5),
          isTRUE(all.equal(unname(V), unname(baseline$V), tolerance = 1e-8)))
      }
      # Package z-test independently verifies coefficient SE and reference law.
      ct <- as.data.frame(clubSandwich::coef_test(fit, vcov = cv, test = 'z'))
      stopifnot(isTRUE(all.equal(unname(ct$SE), unname(sqrt(diag(V))), tolerance = 1e-10)))
      saveRDS(list(key = cache_key, V = V), cache_path)
      rm(cv, ct); gc()
    }
    .validate_covariance(V, fe)
    theta_se <- sqrt(as.numeric(C %*% V %*% C))
    S <- A %*% V %*% t(A); z <- as.vector(A %*% beta)
    stopifnot(qr(S)$rank == nrow(A), min(eigen(S, symmetric = TRUE)$values) > 0)
    Q <- as.numeric(crossprod(z, solve(S, z)))
    key <- paste(name, type, sep = '__')
    results[[key]] <- data.frame(model = name, covariance = type, theta = theta,
      se = theta_se, lwr = theta - qnorm(.975) * theta_se,
      upr = theta + qnorm(.975) * theta_se, z = theta / theta_se,
      p = 2 * pnorm(-abs(theta / theta_se)),
      omnibus_chisq = Q, omnibus_df = nrow(A),
      omnibus_p = pchisq(Q, nrow(A), lower.tail = FALSE),
      inference = 'asymptotic normal / chi-square; participant-clustered')
    coefficients[[key]] <- data.frame(model = name, covariance = type, term = fe,
      estimate = unname(beta), se = sqrt(diag(V)),
      p = 2 * pnorm(-abs(beta / sqrt(diag(V)))))
    write.csv(bind_rows(results), file.path(out, 'theta_and_omnibus.csv'), row.names = FALSE)
    write.csv(bind_rows(coefficients), file.path(out, 'coefficient_tests.csv'), row.names = FALSE)
    if (type == 'CR2') {
      message('  Computing CR2 participant contributions to theta variance')
      score_path <- file.path(out, paste0(name, '_CR2_scores.rds'))
      saved <- if (file.exists(score_path)) readRDS(score_path) else NULL
      if (!is.null(saved) && identical(saved$key, cache_key)) {
        E <- saved$E
      } else {
        E <- fast_cov(type, form = 'estfun')
        saveRDS(list(key = cache_key, E = E), score_path)
      }
      stopifnot(nrow(E) == length(beta), ncol(E) == nlevels(droplevels(frame$id)),
        isTRUE(all.equal(unname(tcrossprod(E)), unname(V), tolerance = 1e-8)))
      ids <- levels(droplevels(frame$id))
      if (!is.null(colnames(E))) stopifnot(identical(colnames(E), ids))
      score <- as.vector(C %*% E)
      stopifnot(abs(sum(score^2) / theta_se^2 - 1) < 1e-8)
      shares <- score^2 / sum(score^2)
      rank <- order(shares, decreasing = TRUE)
      contributions <- data.frame(id = ids, theta_score = score, variance_share = shares) %>%
        left_join(mutate(reuse, id = as.character(id)), by = 'id') %>% arrange(desc(variance_share))
      write.csv(contributions, file.path(out, paste0(name, '_theta_cluster_contributions.csv')), row.names = FALSE)
      concentration[[name]] <- data.frame(model = name, n_clusters = length(shares),
        max_variance_share = max(shares), top_10_variance_share = sum(head(shares[rank], 10)),
        top_1pct_variance_share = sum(head(shares[rank], ceiling(.01 * length(shares)))),
        inverse_variance_share_hhi = 1 / sum(shares^2),
        note = 'Descriptive concentration only; not effective degrees of freedom or a leverage test')
      write.csv(bind_rows(concentration), file.path(out, 'theta_cluster_concentration.csv'), row.names = FALSE)
      rm(E, saved); gc()
    }
  }
  audit[[name]] <- data.frame(model = name, cached_fit_md5 = fit_md5,
    n_rows = nobs(fit), n_clusters = nrow(reuse), fixed_design_matches = TRUE,
    theta_matches_previous = TRUE, no_model_refit = TRUE)
  write.csv(bind_rows(audit), file.path(out, 'verification.csv'), row.names = FALSE)
  rm(fit, frame, fast_cov); gc()
}
res <- bind_rows(results)
a <- res[res$covariance == 'CR2', ]
b <- res[res$covariance == 'CR0', ]; b <- b[match(a$model, b$model), ]
comparison <- data.frame(model = a$model, CR2_se = a$se, CR0_se = b$se,
  pct_se_difference = 100 * (a$se / b$se - 1),
  model_based_se = previous$theta_se[match(a$model, previous$configuration)])
comparison$CR2_to_model_based_se <- comparison$CR2_se / comparison$model_based_se
write.csv(comparison, file.path(out, 'covariance_comparison.csv'), row.names = FALSE)
cat('Completed cached-fit inference diagnostics:', out, '\n')
