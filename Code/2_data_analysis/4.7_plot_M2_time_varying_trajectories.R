# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script:  4.7_plot_M2_time_varying_trajectories.R
# Author:  Nemo Zhou
# Date started:      2026-07-28
# Date last updated: 2026-07-28
#
# Purpose:
#   Visualize all active cancer-cohort trajectories from the restricted M2
#   time-varying-covariate GLMEs. This script reads saved prediction CSVs only;
#   it does not refit models or calculate additional contrasts.
#
# Inputs:
#   Results/cancer/data/4.2.2_m2_time_varying_predicted_trajectories.csv
#   Results/cancer/data/4.2.2_other_m2_time_varying_predicted_trajectories.csv
#
# Outputs:
#   Results/cancer/visuals/4.7_m2_time_varying_trajectories.pdf
#   Results/cancer/visuals/4.7_m2_time_varying_trajectories.html
#
#   Persistent PNG files are not created. The HTML contains the plot as inline
#   SVG so it can be reviewed without refitting the GLMEs.
# =============================================================================

library(dplyr)
library(ggplot2)

project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
data_dir <- file.path(project_dir, "Results", "cancer", "data")
visuals_dir <- file.path(project_dir, "Results", "cancer", "visuals")

prediction_paths <- c(
  smoking = file.path(data_dir, "4.2.2_m2_time_varying_predicted_trajectories.csv"),
  other = file.path(data_dir, "4.2.2_other_m2_time_varying_predicted_trajectories.csv")
)
if (any(!file.exists(prediction_paths))) {
  stop("Required saved prediction files are missing: ",
       paste(names(prediction_paths)[!file.exists(prediction_paths)], collapse = ", "),
       ". Run the two 4.2.2 model scripts first.", call. = FALSE)
}

required_columns <- c(
  "Cohort", "model_id", "model_label", "Group", "Age_Centered",
  "pred", "se", "lwr", "upr", "spline_df", "vcov_type"
)
read_predictions <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns)) {
    stop("Saved predictions in ", basename(path), " are missing: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  x
}

pred <- bind_rows(lapply(prediction_paths, read_predictions))
if (!all(pred$model_id == "M2_time_varying") || !all(pred$spline_df == 4L)) {
  stop("The combined prediction files are not all the requested M2 df-4 outputs.",
       call. = FALSE)
}
if (!all(pred$vcov_type == "model-based vcov(fit)")) {
  stop("The combined prediction files do not all use model-based covariance.",
       call. = FALSE)
}

cohort_order <- c(
  "Low/Moderate Burden Cohort",
  "High Burden Cohort",
  "Smoking-Related Cancer Cohort",
  "Obesity-Related Cancer Cohort",
  "All Cancer Cohort"
)
pred$Cohort <- factor(pred$Cohort, levels = cohort_order)
pred$Group <- factor(pred$Group, levels = c("Control", "Cancer Case"))
if (anyNA(pred$Cohort) || anyNA(pred$Group)) {
  stop("Unexpected cohort or Group labels were found in the saved predictions.",
       call. = FALSE)
}
if (!dir.exists(visuals_dir)) dir.create(visuals_dir, recursive = TRUE)

group_cols <- c("Control" = "#3182bd", "Cancer Case" = "#de2d26")
p_trajectory <- ggplot(
  pred,
  aes(x = Age_Centered, y = pred, group = Group, color = Group, fill = Group)
) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ Cohort, ncol = 2, scales = "free_y") +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_x_continuous(breaks = seq(-20, 20, by = 4)) +
  labs(
    title = "M2 frailty trajectories with time-varying covariates",
    subtitle = "All active cancer cohorts; natural spline df = 4; model-based 95% CIs",
    x = "Years relative to index",
    y = "Predicted frailty index (no-cancer FI)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

pdf_path <- file.path(visuals_dir, "4.7_m2_time_varying_trajectories.pdf")
ggsave(pdf_path, p_trajectory, width = 11, height = 10, units = "in")

svg_path <- tempfile(fileext = ".svg")
svglite::svglite(svg_path, width = 11, height = 10)
print(p_trajectory)
grDevices::dev.off()
svg_markup <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
unlink(svg_path)

html_path <- file.path(visuals_dir, "4.7_m2_time_varying_trajectories.html")
html_text <- c(
  "<!doctype html>",
  "<html lang=\"en\">",
  "<head><meta charset=\"utf-8\"><title>4.7 M2 Time-Varying Trajectories</title></head>",
  "<body>",
  svg_markup,
  "</body></html>"
)
writeLines(html_text, html_path, useBytes = TRUE)

message("Wrote: ", pdf_path)
message("Wrote: ", html_path)
