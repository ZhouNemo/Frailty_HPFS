# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-29
# Date last updated: 2026-09-29 (FI-independent matching and complete assignment ledgers)
# Purpose: Pure eligibility, assignment expansion, preparation, provenance,
# and retired-entry-point contracts. Synthetic temporary files only.
# No sampling, matching builder, or model fitting is executed.
# Usage: Rscript Code/tests/2.0_test_no_fi_eligibility.R
# =============================================================================
source("Code/2_data_analysis/2.0_riskset_matching_functions.R")
source("Code/2_data_analysis/2.0_matching_provenance.R")
source("Code/2_data_analysis/5.0_ipcw_gee_functions.R")
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e),fixed=TRUE))
}
pool <- data.frame(id=c("case","ok","late","dead","cancer","old","future"),
  dob=c(160,160,160,160,160,124,160),first=c(900,1000,1001,900,900,900,900),
  dth=c(NA,NA,NA,1000,NA,NA,1001),canc=c(1000,NA,NA,NA,1000,NA,1001))
f <- riskset_control_flags(pool,"case",1000,70,2)
stopifnot(identical(which(f$eligible),c(2L,7L)))
for (value in c(NA,900,952,1000,1001)) {
  pool$fi_date <- value; pool$n_fi <- 0; pool$participated <- 0; pool$cycle <- "20"
  stopifnot(identical(f,riskset_control_flags(pool,"case",1000,70,2)))
}
# Cases and controls with zero FI survive assignment expansion, including a
# wholly unobserved set; old FI and post-index-only FI remain observed outcomes.
a <- data.frame(Cohort="test",match_set=rep(1:3,each=2),id=letters[1:6],
  role=rep(c("Case","Control"),3),index_date=1000,index_cycle="00")
p <- data.frame(id=letters[1:6],dbmy09=160,dtdth=NA_real_,
  cancer_dateca=rep(c(1000,NA_real_),3),true_age_at_cancer=rep(c(70,NA_real_),3),first_return=800)
b <- data.frame(id=letters[1:6],base_race="a")
o <- data.frame(id=c("b","b","c"),cycle=c("88","92","00"),
  worked_rtmnyr=c(800,952,1001),age_at_cycle=(c(800,952,1001)-160)/12,fi_score_nocancer=.1)
x <- expand_matched_assignments(a,p,b,o,"test",c("88","92","00"),70,1)
stopifnot(nrow(x$assignments)==6,nrow(x$matched_long)==3,
  sum(!x$assignments$contributes_fi)==4,all(x$set_integrity$valid_set),
  !any(x$outcome_support$both_arms_with_fi),
  all(x$matched_long$trajectory_id %in% x$assignments$trajectory_id),
  x$assignments$n_preindex_fi[x$assignments$id=="c"]==0,
  x$assignments$n_recent_preindex_fi[x$assignments$id=="b"]==1,
  x$assignments$n_preindex_fi[x$assignments$id=="b"]==2,
  any(x$matched_long$worked_rtmnyr==800))
expect_error(expand_matched_assignments(a,p,b,rbind(o,o[1,]),"test",c("88","92","00"),70,1),"Duplicated trajectory")
empty <- expand_matched_assignments(a,p,b,o[FALSE,],"test",c("88","92","00"),70,1)
stopifnot(nrow(empty$matched_long)==0,nrow(empty$assignments)==6,all(empty$assignments$n_analytic_fi_visits==0))
s <- complete_assignment_summary(x$assignments,data.frame(trajectory_id=x$assignments$trajectory_id[2],n_model_rows=2),"trajectory_id","n_model_rows")
stopifnot(nrow(s)==6,sum(s$n_model_rows)==2)
check_preparation <- function() {
  path <- tempfile(fileext=".rds"); on.exit(unlink(path))
  d <- data.frame(id=rep(c("a","b"),each=4),cycle=rep(c("86","88","92","96"),2),
    participated=1,worked_rtmnyr=rep(c(1032,1056,1104,1152),2),
    cancer_dateca=1153,cancer_index_dateca=1153,dtdth=NA_real_,dbmy09=300,
    fi_score_nocancer=c(.1,.1,NA,.2,NA,NA,NA,NA),race="White",marital_status="Married",living_arr="With others",pckgr="Never")
  saveRDS(d,path)
  z <- prepare_riskset_inputs(path,c("88","92","96"))
  stopifnot(nrow(z$person_level)==2,all(z$person_level$first_return==1056),
    z$person_level$n_visits[z$person_level$id=="b"]==0,all(z$person_level$is_case==1),
    identical(z$fi_trajectory$cycle,c("88","96")))
  saveRDS(rbind(d,d[8,]),path)
  expect_error(prepare_riskset_inputs(path,c("88","92","96")),"Duplicate participant-cycle")
}
check_preparation()
check_provenance <- function() {
  root <- tempfile(); dir.create(file.path(root,"Data"),recursive=TRUE)
  on.exit(unlink(root,recursive=TRUE))
  path <- file.path(root,"Data/riskset_matched_overall_long.rds")
  ap <- matching_assignment_path(path)
  saveRDS(x$matched_long,path);saveRDS(x$assignments,ap)
  md <- file.path(root,"Results/cancer/data/matching_diagnostics");dir.create(md,recursive=TRUE)
  stem <- "riskset_matched_overall_long"
  write.csv(data.frame(gate_pass=TRUE),file.path(md,paste0(stem,"_gate_g4.csv")),row.names=FALSE)
  mp <- file.path(md,paste0(stem,"_run_metadata.rds"))
  meta <- list(eligibility_version="no_fi_requirement_v1",fi_requirement="none",
    cohort_entry_rule="first_participated_analytic_cycle_return",control_entry_on_or_before_index=TRUE,
    output_md5=unname(tools::md5sum(path)),assignment_md5=unname(tools::md5sum(ap)))
  saveRDS(meta,mp);stopifnot(nrow(read_matching_assignments(path))==6)
  ipcw_gee_validate_matching_provenance(path)
  derived <- x$matched_long
  attr(derived,"assignment_provenance") <- validate_matching_provenance(path)[c("assignment_path","assignment_md5","input_md5")]
  validate_derived_assignment_provenance(derived)
  bad_derived <- derived;attr(bad_derived,"assignment_provenance")$input_md5 <- "stale"
  expect_error(validate_derived_assignment_provenance(bad_derived),"stale primary")
  for(v in c("legacy","recent_fi_v1")) {
    bad <- meta;bad$eligibility_version <- v;saveRDS(bad,mp)
    expect_error(validate_matching_provenance(path),"2.5_create_riskset_overall.R")
    expect_error(ipcw_gee_validate_matching_provenance(path),"2.5_create_riskset_overall.R")
  }
  saveRDS(meta,mp);saveRDS(x$assignments[1,],ap)
  expect_error(validate_matching_provenance(path),"ledger hash")
  saveRDS(x$assignments,ap);saveRDS(x$matched_long[1,],path)
  expect_error(validate_matching_provenance(path),"RDS hash")
  expect_error(validate_matching_provenance("riskset_matched_overall_exact_cycle_long.rds"),"retired")
  expect_error(validate_matching_provenance("riskset_matched_analysis_cancer_free_full_endpoint_long.rds"),"retired")
}
check_provenance()
assert_assignment_rows(x$assignments)
bad <- x$assignments;bad$first_return[1] <- 1001
expect_error(assert_assignment_rows(bad),"cohort-entry/version")
engine <- paste(readLines("Code/2_data_analysis/2.0_riskset_matching_functions.R"),collapse="\n")
stopifnot(!grepl("fi_lookback_months|min_visits|preindex_support|cycle_caliper",engine))
for(p in list.files("Code/2_data_analysis",pattern="^2[.][1-5]_create",full.names=TRUE))
  stopifnot(!any(grepl("FI_LOOKBACK|MIN_VISITS|exact_cycle|CYCLE_CALIPER",readLines(p))))
retired <- c("9.0_cancer_free_matching_functions.R","9.1_create_cancer_free_high_low_burden.R",
             "10.0_GLME_cancer_free_functions.R","10.1_GLME_cancer_free_high_low_burden.R")
stopifnot(!any(file.exists(file.path("Code/2_data_analysis",retired))))
source("Code/2_data_analysis/4.0.1_GLME_sensitivity_analyses.R")
expect_error(run_selected_glme_sensitivities("S7"),"retired")
expect_error(run_overall_glme_sensitivities("unused","unused",selected="S7"),"retired")
cat("PASS: no-FI eligibility, complete assignments, preparation, provenance, retirement; no matching/fitting.\n")
# Preparation only: no response regression is fitted. Inject a complete ledger
# into the assignment-panel reader; file/hash validation was exercised above.
check_ipcw_panel <- function() {
  original <- read_matching_assignments
  ledger <- x$assignments
  for(v in setdiff(.ipcw_gee_m2_covars,names(ledger))) ledger[[v]] <- 1
  assign("read_matching_assignments",function(...) ledger,envir=.GlobalEnv)
  on.exit(assign("read_matching_assignments",original,envir=.GlobalEnv))
  response <- tidyr::crossing(id=ledger$id,ipcw_gee_cycle_dates())
  response$fi_score_nocancer <- NA_real_
  response$age_at_cycle <- NA_real_;response$age_nom <- (response$cycle_date-160)/12
  response$alive_at_cycle <- TRUE;response$cancer_index_dateca <- NA_real_
  response$cancer_now <- 0L;response$years_since_cancer <- 0
  response$lag_fi <- NA_real_;response$prior_response_count <- 0L;response$observed <- 0L
  panel <- ipcw_gee_build_assignment_panel("unused","test",list(ledger=response))
  stopifnot(n_distinct(panel$trajectory_id)==6,nrow(panel)==6*9,
            all(!panel$fi_observed),all(is.na(panel$fi_score_nocancer)))
}
check_ipcw_panel()
cat("PASS: IPCW panel retains all zero-FI assignment denominators; no fitting.\n")
