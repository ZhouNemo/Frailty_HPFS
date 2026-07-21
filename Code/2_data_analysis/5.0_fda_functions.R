# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  5.0_fda_functions.R
# Author:  Nemo Zhou
# Date started:      2026-06-29
# Date last updated: 2026-07-17 (visual-data outputs moved; PNG writes removed)
#
# Purpose:
#   Shared functional data analysis (FDA) utilities used by 5.1+. Converts a
#   risk-set matched long dataset (from 2.1-2.5) into the SPARSE functional-data
#   format required by fdapace (the PACE / sparse-FPCA method), and optionally
#   runs FPCA with a Cancer-vs-Control mean-function comparison.
#
#   Why sparse FPCA: HPFS frailty is measured on a few irregular questionnaire
#   cycles per person, and relative time (Age_Centered) differs across people
#   because the index date varies. This is sparse, irregularly sampled functional
#   data -- the setting fdapace::FPCA is built for. It does NOT require a common
#   time grid, so no binning/interpolation is fabricated (contrast the dense
#   fda:: matrix format, which would).
#
#   Functional unit = trajectory ASSIGNMENT, keyed Cohort__match_set__id__role,
#   so a control reused across match sets contributes as distinct curves (matches
#   the unit used in 1.1_check_timebin_fda_support.R).
#
#   Required fdapace input:
#     Ly : list of FI value vectors, one per curve
#     Lt : list of matching relative-time (Age_Centered) vectors, one per curve
#   plus a `meta` table (trajectory_id, Group, Cohort, id) aligned to names(Ly).
#
# Units:
#   Age_Centered is years relative to each person's own attained age at the
#   assigned index (0 = index). Curve value is fi_score_nocancer.
# =============================================================================

library(dplyr)

# ---- 1) Convert matched long data -> fdapace Lt/Ly/meta --------------------
prepare_fda_input <- function(matched_path,
                              out_path = NULL,
                              min_measures = 3,
                              window_yrs = NULL,
                              builder_script = NULL) {

  req <- c("Cohort", "match_set", "id", "role", "Group",
           "Age_Centered", "fi_score_nocancer")

  if (!file.exists(matched_path)) {
    stop("Matched dataset not found at ", matched_path,
         if (!is.null(builder_script)) paste0(". Run ", builder_script, " first.") else "")
  }
  ml <- readRDS(matched_path)
  miss <- setdiff(req, names(ml))
  if (length(miss) > 0) stop("Matched dataset missing columns: ", paste(miss, collapse = ", "))

  ml <- ml %>%
    mutate(
      trajectory_id = paste(Cohort, match_set, id, role, sep = "__"),
      Age_Centered = as.numeric(Age_Centered),
      fi_score_nocancer = as.numeric(fi_score_nocancer)
    ) %>%
    filter(!is.na(Age_Centered), !is.na(fi_score_nocancer))

  if (!is.null(window_yrs)) ml <- ml %>% filter(abs(Age_Centered) <= window_yrs)

  # one value per (curve, time): average any tied relative times
  fda_long <- ml %>%
    group_by(trajectory_id, Age_Centered) %>%
    summarize(fi = mean(fi_score_nocancer),
              Group = first(Group), Cohort = first(Cohort),
              id = first(id), .groups = "drop") %>%
    arrange(trajectory_id, Age_Centered)

  # FPCA needs enough points per curve to estimate within-curve shape
  counts <- fda_long %>% count(trajectory_id, name = "n_meas")
  keep   <- counts %>% filter(n_meas >= min_measures) %>% pull(trajectory_id)
  n_drop <- nrow(counts) - length(keep)
  fda_long <- fda_long %>% filter(trajectory_id %in% keep)

  sp <- split(fda_long, fda_long$trajectory_id)
  Ly <- lapply(sp, `[[`, "fi")
  Lt <- lapply(sp, `[[`, "Age_Centered")

  meta <- fda_long %>%
    distinct(trajectory_id, Group, Cohort, id)
  meta <- meta[match(names(Ly), meta$trajectory_id), , drop = FALSE]

  fda_input <- list(
    Ly = Ly, Lt = Lt, meta = meta,
    min_measures = min_measures, source = matched_path,
    built = Sys.time()
  )

  cat("\n==== FDA input prepared from:", basename(matched_path), "====\n")
  cat("Curves kept (>=", min_measures, "measures):", length(Ly),
      " | dropped (too sparse):", n_drop, "\n")
  cat("Curves by Group:\n"); print(table(meta$Group))
  cat("Measurements per curve (kept):\n")
  print(summary(lengths(Ly)))
  cat("Relative-time range across curves: [",
      round(min(unlist(Lt)), 1), ",", round(max(unlist(Lt)), 1), "]\n")

  if (!is.null(out_path)) {
    saveRDS(fda_input, out_path)
    cat("Saved fdapace input (Ly/Lt/meta) to:", out_path, "\n")
  }
  invisible(fda_input)
}

# ---- 2) Optional: run sparse FPCA + Cancer-vs-Control mean comparison -------
run_fpca <- function(fda_input,
                     results_dir,
                     out_prefix,
                     fve_threshold = 0.95,
                     max_curves_per_group = NULL,
                     seed = 2026) {

  if (!requireNamespace("fdapace", quietly = TRUE)) {
    message("fdapace not installed; skipping FPCA. Run install.packages('fdapace').")
    return(invisible(NULL))
  }
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  Ly <- fda_input$Ly; Lt <- fda_input$Lt; meta <- fda_input$meta

  # optional down-sampling per group (FPCA covariance smoothing is heavy on
  # tens of thousands of sparse curves; controls usually dominate)
  if (!is.null(max_curves_per_group)) {
    set.seed(seed)
    idx <- unlist(lapply(split(seq_len(nrow(meta)), meta$Group), function(ii) {
      if (length(ii) > max_curves_per_group) sample(ii, max_curves_per_group) else ii
    }))
    idx <- sort(idx)
    Ly <- Ly[idx]; Lt <- Lt[idx]; meta <- meta[idx, , drop = FALSE]
    cat("Down-sampled to", length(Ly), "curves (<=", max_curves_per_group, "per group).\n")
  }

  optns <- list(dataType = "Sparse", error = TRUE, methodSelectK = "FVE",
                FVEthreshold = fve_threshold, nRegGrid = 51, verbose = FALSE)

  # ---- overall FPCA: modes of variation + scores for downstream models ----
  cat("\nFitting overall FPCA on", length(Ly), "curves ...\n")
  fpca_all <- fdapace::FPCA(Ly, Lt, optns = optns)
  K <- fpca_all$selectK
  cat("Selected", K, "components; cumulative FVE:",
      round(utils::tail(fpca_all$cumFVE, 1), 1), "%\n")

  scores <- as.data.frame(fpca_all$xiEst)
  names(scores) <- paste0("FPC", seq_len(ncol(scores)))
  scores <- bind_cols(meta, scores)
  write.csv(scores, file.path(results_dir, paste0(out_prefix, "_fpca_scores.csv")),
            row.names = FALSE)

  # mean function + eigenfunctions on the regular work grid
  eig <- data.frame(rel_time = fpca_all$workGrid, mean = fpca_all$mu)
  for (k in seq_len(K)) eig[[paste0("phi", k)]] <- fpca_all$phi[, k]
  write.csv(eig, file.path(results_dir, paste0(out_prefix, "_fpca_eigenfunctions.csv")),
            row.names = FALSE)

  # ---- group-stratified mean functions: Cancer vs Control comparison ----
  grp_mu <- lapply(c("Control", "Cancer Case"), function(g) {
    ii <- which(meta$Group == g)
    if (length(ii) < 10) return(NULL)
    f <- fdapace::FPCA(Ly[ii], Lt[ii], optns = optns)
    data.frame(Group = g, rel_time = f$workGrid, mean_fi = f$mu)
  })
  grp_mu <- bind_rows(grp_mu)
  write.csv(grp_mu, file.path(results_dir, paste0(out_prefix, "_fpca_mean_by_group.csv")),
            row.names = FALSE)

  group_mean_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(grp_mu) > 0) {
    gcol <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")
    group_mean_plot <- ggplot2::ggplot(grp_mu, ggplot2::aes(rel_time, mean_fi, color = Group)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_color_manual(values = gcol) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(title = "FDA (sparse FPCA) mean frailty function",
                    subtitle = "PACE mean curve by group; centered on own attained age at index",
                    x = "Years relative to index", y = "Mean frailty index (no-cancer FI)",
                    color = NULL)
  }

  cat("Saved FPCA scores, eigenfunctions, and group mean functions to:", results_dir, "\n")
  invisible(list(
    fpca_all = fpca_all,
    scores = scores,
    group_means = grp_mu,
    group_mean_plot = group_mean_plot
  ))
}
