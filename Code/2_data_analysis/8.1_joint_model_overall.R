# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Script: 8.1_joint_model_overall.R
# Author: Nemo Zhou
# Date started: 2026-06-30
# Date last updated: 2026-09-29
# Purpose: Overall joint longitudinal-survival mortality sensitivity aligned
# with the canonical 4.5 M2 df-3 spline, nine baseline covariates, and
# attained-age random slope. One assignment
# per participant is selected before missing-data exclusions: case when present
# in the matched ledger, otherwise earliest control (deterministic set ties).
# Index is time zero; pre-index FI remains history, mortality starts at index,
# and follow-up ends at death, own cancer for controls, or the fixed 2020 date.
# Preparation validates data/provenance and writes recoverable input artifacts.
# All profiles currently stop after audits: the month-resolution date policy
# remains unresolved. No runtime override; no statistical fitter is reached.
# After a separately reviewed policy revision, Nemo runs pilot/full analyses.
#
# Required inputs: canonical matched overall RDS and FI panel under Data;
# matching Gate G4, hashes, age scaling, and schema-3 overall metadata with M2 specification and scaling
# under Results/cancer/data. See 8.0_joint_model_functions.R for contracts.
#
# Interface: JM_PROFILE=prepare (default), pilot, or full. Invalid values stop.
# Outputs: Results/cancer/data/8.1_joint_model_overall/<profile>/<unique-run>/.
# Each run preserves previous results. Report: 8.1_joint_model_report.Rmd;
# render to Results/cancer/visuals using its run_dir parameter. No PNG writes.
# Source-safe: tests may source this wrapper without invoking the workflow;
# execution is via Rscript or an explicit jm_run(project_dir, profile) call.
# =============================================================================
project_dir <- "/Users/nemo/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Frailty HPFS"
source(file.path(project_dir, "Code", "2_data_analysis", "8.0_joint_model_functions.R"))
if (sys.nframe() == 0L) jm_run(project_dir, jm_profile())
