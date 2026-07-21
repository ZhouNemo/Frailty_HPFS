# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.2_plot_saved_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-07-20
# Date last updated: 2026-07-20
#
# Purpose:
#   Create trajectory plots from the saved 4.2 prediction CSV without fitting
#   any models. The plot includes M0 raw spline, M1 primary spline, M2 full
#   spline, and M3 matching-set spline predictions for the smoking-related
#   cancer risk-set matched cohort.
#
# Inputs:
#   Results/cancer/data/4.2_spline_predicted_trajectories.csv
#   Results/cancer/data/4.2_spline_model_status.csv
#
# Outputs:
#   Results/cancer/visuals/4.2_saved_trajectories_M0_M1_M2_M3.pdf
#   Results/cancer/visuals/4.2_saved_trajectories_M0_M1_M2_M3.html
#
#   Persistent PNG files are not created. The HTML contains an inline SVG so
#   the visual is self-contained and can be reviewed without refitting models.
# =============================================================================

library(ggplot2)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Results", "cancer", "data")
visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")

prediction_path <- file.path(data_dir, "4.2_spline_predicted_trajectories.csv")
status_path <- file.path(data_dir, "4.2_spline_model_status.csv")

if (!file.exists(prediction_path)) {
  stop("Saved trajectory predictions not found: ", prediction_path, call. = FALSE)
}

if (!file.exists(status_path)) {
  stop("Saved model-status table not found: ", status_path, call. = FALSE)
}

status <- read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
pred <- read.csv(prediction_path, stringsAsFactors = FALSE, check.names = FALSE)

required_prediction_columns <- c(
  "Cohort", "model_id", "model_label", "model_role", "Group",
  "Age_Centered", "pred", "lwr", "upr"
)
missing_prediction_columns <- setdiff(required_prediction_columns, names(pred))
if (length(missing_prediction_columns)) {
  stop(
    "Saved trajectory predictions are missing required columns: ",
    paste(missing_prediction_columns, collapse = ", "),
    call. = FALSE
  )
}

model_order <- c(
  M0_raw = "M0 raw",
  M1_primary_spline = "M1 primary spline",
  M2_full_spline = "M2 full spline",
  M3_primary_matching_spline = "M3 matching-set spline"
)

missing_models <- setdiff(names(model_order), unique(pred$model_id))
if (length(missing_models)) {
  status_cols <- intersect(
    c("Cohort", "model_id", "status", "n_obs", "prediction_min_time",
      "prediction_max_time", "message"),
    names(status)
  )
  status_detail <- status[status$model_id %in% missing_models, status_cols, drop = FALSE]
  status_text <- paste(capture.output(print(status_detail, row.names = FALSE)), collapse = "\n")
  stop(
    "The saved prediction CSV does not yet contain all requested models.\n",
    "Missing: ", paste(missing_models, collapse = ", "), "\n",
    "Models currently present: ", paste(sort(unique(pred$model_id)), collapse = ", "),
    "\n\nCurrent model-status detail:\n", status_text, "\n\n",
    "Rerun Code/2_data_analysis/4.2_GLME_smoking_related.R, then rerun this plotting script.",
    call. = FALSE
  )
}

if (file.info(status_path)$mtime > file.info(prediction_path)$mtime) {
  warning(
    "The model-status file is newer than the prediction file. Confirm that the ",
    "prediction CSV was regenerated from the same 4.2 run before interpreting the plot.",
    call. = FALSE
  )
}

trajectory_predictions <- pred[pred$model_id %in% names(model_order), , drop = FALSE]
trajectory_predictions$Model <- factor(
  unname(model_order[trajectory_predictions$model_id]),
  levels = unname(model_order)
)
trajectory_predictions$Cohort <- factor(
  trajectory_predictions$Cohort,
  levels = unique(trajectory_predictions$Cohort)
)
trajectory_predictions$Group <- factor(
  trajectory_predictions$Group,
  levels = c("Control", "Cancer Case")
)

if (anyNA(trajectory_predictions$Model) || anyNA(trajectory_predictions$Group)) {
  stop("Unexpected model or Group labels were found in the saved prediction CSV.", call. = FALSE)
}

p_trajectory <- ggplot(
  trajectory_predictions,
  aes(
    x = Age_Centered,
    y = pred,
    group = Group,
    color = Group,
    fill = Group
  )
) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  facet_grid(Cohort ~ Model) +
  scale_color_manual(values = c("Control" = "#3182bd", "Cancer Case" = "#de2d26")) +
  scale_fill_manual(values = c("Control" = "#3182bd", "Cancer Case" = "#de2d26")) +
  scale_x_continuous(breaks = seq(-20, 20, by = 4)) +
  labs(
    title = "Saved frailty trajectories: 4.2 smoking-related cancer",
    subtitle = "Predicted no-cancer frailty index; shaded areas are saved 95% confidence intervals",
    x = "Years relative to index",
    y = "Predicted frailty index",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

if (!dir.exists(visuals_dir)) dir.create(visuals_dir, recursive = TRUE)

pdf_path <- file.path(visuals_dir, "4.2_saved_trajectories_M0_M1_M2_M3.pdf")
ggsave(pdf_path, p_trajectory, width = 16, height = 8, units = "in")

# Render the same plot as a self-contained HTML file with inline SVG.
svg_path <- tempfile(fileext = ".svg")
svglite::svglite(svg_path, width = 16, height = 8)
print(p_trajectory)
grDevices::dev.off()
svg_markup <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
unlink(svg_path)

html_path <- file.path(visuals_dir, "4.2_saved_trajectories_M0_M1_M2_M3.html")
html_text <- c(
  "<!doctype html>",
  "<html lang=\"en\">",
  "<head><meta charset=\"utf-8\"><title>4.2 Saved Trajectories M0 M1 M2 M3</title></head>",
  "<body>",
  svg_markup,
  "</body></html>"
)
writeLines(html_text, html_path, useBytes = TRUE)

message("Wrote: ", pdf_path)
message("Wrote: ", html_path)
