# =============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the
#          Health Professionals Follow-up Study
# Author: Nemo Zhou
# Date started: 2026-09-20
# Date last updated: 2026-09-29 (M2 alignment, date-policy block, complete posterior diagnostics)
# Purpose: Model-free regression tests for 8.1 selection, censoring, dynamic
# spline/hazard-time evaluation, metadata integrity, posterior schemas,
# diagnostics, output isolation, and preparation checkpoints. Synthetic data
# only; never invokes lme, coxph, jm, or risk-set matching. Diagnostic artifacts
# live under Codex/jm_m2_revision_2026-09-29, not analytic dataset/result locations.
# Run from project root with Rscript Code/tests/8.1_test_joint_model_contract.R.
# =============================================================================
project_dir <- normalizePath(getwd())
source(file.path(project_dir,"Code","2_data_analysis","8.0_joint_model_functions.R"))
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL},error=function(e) conditionMessage(e))
  stopifnot(!is.null(e),grepl(pattern,e,fixed=TRUE))
}
# Test-only convenience: production selection requires an explicit complete ledger.
fixture_assignments <- function(d) {
  d$trajectory_id <- paste(d$Cohort,d$match_set,d$id,d$role,sep="__")
  d$first_return <- 1000
  d$eligibility_version <- "no_fi_requirement_v1"
  unique(d[c("trajectory_id","id","Cohort","Group","role","match_set","index_date",
             "index_age","cancer_dateca",jm_covars,"first_return","eligibility_version")])
}
production_prepare <- jm_prepare
production_ledger <- jm_assignment_ledger
jm_prepare <- function(matched,panel,assignments=fixture_assignments(matched))
  production_prepare(matched,panel,assignments)
jm_assignment_ledger <- function(d,assignments=fixture_assignments(d))
  production_ledger(d,assignments)
passed <- character()
check <- function(name,code) {force(code); passed <<- c(passed,name); cat("PASS:",name,"\n")}
make_assignment <- function(id,role,index,cancer,match_set,lev="a") {
  t <- c(-8,-4,0,4,8,12)
  data.frame(id=id,Cohort="All Cancer Cohort",role=role,
    Group=if(role=="Case") "Cancer Case" else "Control",match_set=as.character(match_set),
    cycle=c("88","92","96","00","04","08"),index_date=index,index_age=70,
    worked_rtmnyr=index+12*t,Age_Centered=t,fi_score_nocancer=0.1+t/1000,
    age_at_cycle=70+t,index_age_z=0,base_race=lev,base_marital=lev,base_living=lev,base_pckgr=lev,
    base_calor=2000+match(id,LETTERS)*100,base_sat=20+match(id,LETTERS),
    base_diet_chol=200+match(id,LETTERS)*10,base_alco=5+match(id,LETTERS),
    cancer_dateca=cancer,post_own_cancer=role=="Control" & !is.na(cancer) & index+12*t>=cancer)
}
d <- dplyr::bind_rows(make_assignment("A","Control",1164,1200,1),
                     make_assignment("A","Case",1200,1200,2),
                     make_assignment("B","Control",1200,NA,10,"b"),
                     make_assignment("B","Control",1200,NA,3,"b"),
                     make_assignment("C","Case",1200,1200,4,"b"),
                     make_assignment("D","Control",1200,1272,5))
panel <- data.frame(id=c("A","B","C","D"),dtdth=c(1320,NA,1446,1272))
t <- seq(-15,20,length.out=80)
b <- splines::ns(t,knots=c(-4,4),Boundary.knots=c(-20,30))
m <- list(schema_version=3L,cohort="All Cancer Cohort",outcome="fi_score_nocancer",
  primary_model_id="M1_primary_spline",spline_reference_model_id="M1_primary_spline",
  spline_df=3L,knots=c(-4,4),boundary_knots=c(-20,30),
  primary_covariates=jm_primary_covars,random_clock="attained_age",random_time_center=60,random_time_scale=4,
  model_specification=list(M2_full_spline=list(covars=jm_covars,basis_key="adjusted_spline",
    time_structure="natural spline",spline_df=3L,random="(1 + .random_time | id)",matching_set_random=FALSE)),
  model_scaling=list(M2_full_spline=data.frame(column=c(paste0("S",1:3),jm_diet_covars),
    center=c(colMeans(b),2000,20,200,5),scale=c(apply(b,2,sd),500,10,100,10))))

# No fitted objects: independent synthetic chains and explicit package schemas.
mock_joint <- function(bn,gn,structure="value_slope",diagonal=FALSE,n=1500L) {
  ch <- function(nms,positive=FALSE) lapply(1:3,function(i) {
    x<-matrix(rnorm(n*length(nms),sd=.03),n,length(nms),dimnames=list(NULL,nms))
    if(positive) x<-exp(x)
    x
  })
  an<-if(structure=="value") "value(fi_score_nocancer)" else c("value(fi_score_nocancer)","slope(fi_score_nocancer)")
  dc<-ch(c("D[1, 1]","D[2, 1]","D[2, 2]"))
  for(i in seq_along(dc)) {dc[[i]][,c(1,3)]<-exp(dc[[i]][,c(1,3)]);if(diagonal)dc[[i]][,2]<-0}
  z<-list(mcmc=list(betas1=ch(bn),gammas=ch(gn),alphas=ch(an),D=dc,
    sigmas=ch("sigmas_1",TRUE),bs_gammas=ch(paste0("bs_gammas_",1:4)),tau_bs_gammas=ch("tau_bs_gammas",TRUE)),
    model_data=list(U_H=list(matrix(0,1,length(an),dimnames=list(NULL,an))),W0_H=matrix(0,2,4),
      ind_zero_D=if(diagonal)matrix(c(1L,2L),1,2) else matrix(integer(),0,2)),
    initial_values=list(D=diag(2)),priors=list(test_double=TRUE),
    control=c(jm_hazard_control,list(knots=list(seq(0,20,length.out=4)))),
    fit_stats=list(marginal=list(DIC=100,WAIC=101,LPML=-52)))
  for(v in c("Wlong_bar","Wlong_sds","Wlong_std","W_bar","W_sds","W_std"))z[[v]]<-matrix(1,1,1)
  z
}

check("profile validation", {
  stopifnot(jm_profile("prepare")=="prepare")
  expect_error(jm_profile("FULL"),"JM_PROFILE")
})
check("case priority and numeric tie break independent of row order", {
  a <- jm_assignment_ledger(d)$retained
  stopifnot(a$role[a$id=="A"]=="Case",a$match_set[a$id=="B"]=="3")
  set.seed(4); a2 <- jm_assignment_ledger(d[sample(nrow(d)),])$retained
  stopifnot(identical(a$trajectory_id,a2$trajectory_id))
})
check("death dates, endpoint hierarchy, and no duplicated survival", {
  p <- jm_prepare(d,panel)
  s <- p$survival
  stopifnot(!anyDuplicated(s$id),s$death_event[s$id=="C"]==1,
            s$endpoint_reason[s$id=="D"]=="own_cancer",s$death_event[s$id=="D"]==0,
            s$cancer_death_tie[s$id=="D"],s$surv_time[s$id=="D"]==6,
            any(p$long$Age_Centered<0),all(p$long$Age_Centered<=p$long$surv_time))
  stopifnot(!any(p$long$id=="D" & p$long$Age_Centered>=6))
  z <- panel; z$dtdth[z$id=="A"] <- 1200
  expect_error(jm_prepare(d,z),"Death on/before index")
  expect_error(jm_prepare(d,rbind(panel,data.frame(id="A",dtdth=1350))),"Conflicting death dates")
  late <- jm_assignment_ledger(d)$retained
  late$index_date[late$id=="B"] <- jm_admin_date
  expect_error(jm_survival(late,jm_deaths(panel)),"Nonpositive index-to-endpoint")
})
check("missingness cannot replace selected assignment", {
  z <- d; z$fi_score_nocancer[z$id=="A" & z$role=="Case"] <- NA
  a <- jm_assignment_ledger(z)$retained
  stopifnot(a$role[a$id=="A"]=="Case")
  pp <- jm_prepare(z,panel)
  stopifnot(!"A" %in% as.character(pp$survival$id),"A" %in% pp$excluded_people$id)
  z <- d; z$base_living[z$id=="B" & z$match_set=="3"] <- NA
  pp <- jm_prepare(z,panel)
  stopifnot(!"B" %in% as.character(pp$survival$id))
})
check("case priority persists when case has no outcome rows", {
  complete <- fixture_assignments(d)
  observed <- d[!(d$id=="A" & d$role=="Case"),]
  pp <- jm_prepare(observed,panel,complete)
  stopifnot(pp$retained$role[pp$retained$id=="A"]=="Case",
            "A" %in% pp$excluded_people$id, !"A" %in% as.character(pp$long$id))
})
check("duplicates and clock inconsistency fail", {
  expect_error(jm_prepare(rbind(d,d[1,]),panel),"Duplicate assignment-cycle")
  z <- d; z$Age_Centered[7] <- z$Age_Centered[7]+1
  expect_error(jm_prepare(z,panel),"clocks do not share")
  z <- d; z$post_own_cancer[7] <- TRUE
  expect_error(jm_prepare(z,panel),"flag disagrees")
})
check("metadata corruption fails", {
  jm_validate_metadata(m)
  z <- m; z$spline_df <- 4L
  expect_error(jm_validate_metadata(z),"current overall M2")
  z <- m; z$model_scaling <- NULL
  expect_error(jm_validate_metadata(z),"incomplete")
  z <- m; z$knots <- c(4,-4)
  expect_error(jm_validate_metadata(z),"Invalid shared M1-reference spline")
})
check("dynamic spline and installed JM negative-time contract", {
  z <- jm_check_clock_basis(m)
  stopifnot(z$value_max_error<1e-8,z$slope_max_error<1e-7)
  f <- jm_formula(m); bfun <- environment(f)$jm_time_basis
  points <- c(-10.5,-3.2,0.7,9.1)
  raw <- splines::ns(points,knots=m$knots,Boundary.knots=m$boundary_knots)
  target <- sweep(sweep(raw,2,m$model_scaling$M2_full_spline$center[1:3],"-"),2,m$model_scaling$M2_full_spline$scale[1:3],"/")
  stopifnot(max(abs(bfun(points)-target))<1e-12)
  tf <- tempfile(fileext=".rds"); saveRDS(f,tf); f2 <- readRDS(tf)
  stopifnot(max(abs(environment(f2)$jm_time_basis(points)-target))<1e-12)
})
check("support gating and endpoint theta equivalence", {
  supported <- expand.grid(id=1:300,Age_Centered=seq(-19,19,2))
  supported$Group <- ifelse(supported$id<=50,"Cancer Case","Control")
  s <- jm_support(supported)
  stopifnot(identical(s$window,c(-20,20)),s$theta_supported)
  supported <- supported[!(supported$Group=="Cancer Case" & supported$Age_Centered==7),]
  stopifnot(!jm_support(supported)$theta_supported)
  p <- jm_prepare(d,panel)
  design <- jm_design(jm_scale_data(p,m),jm_formula(m),list(window=c(-20,20),theta_supported=TRUE))
  coef <- seq_along(design$theta)/100
  stopifnot(abs(sum((design$theta-design$numeric_theta)*coef))<1e-6)
})
check("chain schema, extraction, and convergence gating", {
  p <- jm_prepare(d,panel)
  design <- jm_design(jm_scale_data(p,m),jm_formula(m),list(window=c(-20,20),theta_supported=TRUE))
  bn <- names(design$theta); gn <- c("GroupCancer Case","index_age_z")
  chains <- function(nms,n=1000,shift=0) lapply(1:3,function(i) {
    x <- matrix(rnorm(n*length(nms),mean=shift*(i-1)),ncol=length(nms)); colnames(x)<-nms; x
  })
  set.seed(2026)
  fake <- mock_joint(bn,gn,n=1000L)
  draws <- jm_posterior(fake,bn,gn,"value_slope",3)
  stopifnot(all(jm_diagnostics(draws,"full")$pass),!any(jm_diagnostics(draws,"pilot")$pass))
  z <- fake; z$mcmc$betas1[[2]] <- z$mcmc$betas1[[2]][,-1]
  expect_error(jm_posterior(z,bn,gn,"value_slope",3),"Incomplete/nonfinite")
  z <- fake; z$mcmc$alphas <- NULL
  expect_error(jm_posterior(z,bn,gn,"value_slope",3),"chain count")
  z <- draws; z$alphas[,2,1] <- z$alphas[,2,1]+5
  stopifnot(!all(jm_diagnostics(z,"full")$pass))
  flat <- jm_flatten(draws$betas)
  stopifnot(identical(unname(flat[1:1000,]),unname(fake$mcmc$betas1[[1]])))
  sums <- jm_joint_summaries(draws,design,list(theta=data.frame(estimate=0)))
  stopifnot(nrow(sums$curves)==4*length(design$times),is.finite(sums$theta$estimate))
  design$theta_supported <- FALSE
  stopifnot(is.na(jm_joint_summaries(draws,design,list(theta=data.frame(estimate=0)))$theta$estimate))
})
check("all profiles and direct fitting adapters blocked after complete audits", {
  diagnostic_root <- file.path(project_dir,"Codex","jm_m2_revision_2026-09-29")
  dir.create(diagnostic_root,recursive=TRUE,showWarnings=FALSE)
  fixture_root <- tempfile("blocked_",diagnostic_root); dir.create(fixture_root)
  original_validate <- jm_validate_inputs; original_fit <- jm_fit_components; original_link <- jm_fit_linked
  # Two distinct and overlapping date flags, plus missing M2 covariate.
  td <- d; td$base_alco[td$id=="A"]<-NA_real_
  tp <- panel; tp$dtdth[tp$id=="A"]<-1200
  jm_validate_inputs <- function(...) list(matched=td,panel=tp,assignments=fixture_assignments(td),metadata=m,
    provenance=list(hashes=c(fixture="synthetic")))
  calls <- 0L
  jm_fit_components <- jm_fit_linked <- function(...) {calls<<-calls+1L;stop("FITTER MUST NOT BE REACHED")}
  for(profile in c("prepare","pilot","full")) {
    expect_error(jm_run(fixture_root,profile),"Fitting blocked")
    run <- list.dirs(file.path(fixture_root,"Results","cancer","data","8.1_joint_model_overall",profile),recursive=FALSE)
    st <- readRDS(file.path(run,"run_status.rds"))
    stopifnot(st$status=="blocked_date_policy",st$alignment_model_id==jm_alignment_model,
      file.exists(file.path(run,"date_audit.csv")),file.exists(file.path(run,"M2_covariate_availability.csv")),
      !file.exists(file.path(run,"prepared_data.rds")),!file.exists(file.path(run,"components.rds")))
    audit <- read.csv(file.path(run,"date_audit.csv"));stopifnot(audit$same_month_death[audit$id=="A"],!audit$m2_complete_covariates[audit$id=="A"])
  }
  stopifnot(calls==0L)
  jm_validate_inputs <- original_validate; jm_fit_components <- original_fit; jm_fit_linked <- original_link
  expect_error(jm_fit_components(NULL,NULL,NULL),"Fitting blocked")
  expect_error(jm_fit_linked(NULL,NULL,NULL,NULL),"Fitting blocked")
  expect_error(jm_validate_inputs(fixture_root),"Missing prerequisite")
})
check("Gate G4, source hashes, and stale M2 metadata fail closed", {
  base <- tempfile("provenance_",file.path(project_dir,"Codex","jm_m2_revision_2026-09-29")); dir.create(base)
  data_dir <- file.path(base,"Data"); dir.create(data_dir)
  result_dir <- file.path(base,"Results","cancer","data"); dir.create(result_dir,recursive=TRUE)
  md_dir <- file.path(result_dir,"matching_diagnostics"); dir.create(md_dir)
  stem <- "riskset_matched_overall_long"
  matched_path <- file.path(data_dir,paste0(stem,".rds"))
  panel_path <- file.path(data_dir,"FI_longitudinal_1986_2020_IMPUTED_Cancer.rds")
  saveRDS(d,matched_path); saveRDS(panel,panel_path)
  ref <- d[!d$post_own_cancer,]
  mm <- m
  mm$knots <- as.numeric(quantile(ref$Age_Centered,c(1/3,2/3)))
  mm$boundary_knots <- range(ref$Age_Centered)
  bb <- splines::ns(ref$Age_Centered,knots=mm$knots,Boundary.knots=mm$boundary_knots)
  mm$model_scaling$M2_full_spline$center <- c(colMeans(bb),colMeans(ref[jm_diet_covars]))
  mm$model_scaling$M2_full_spline$scale <- c(apply(bb,2,sd),vapply(ref[jm_diet_covars],sd,numeric(1)))
  mm$n_obs <- nrow(ref); mm$n_id <- length(unique(ref$id))
  saveRDS(mm,file.path(result_dir,"glme_spline_metadata.rds"))
  ap <- file.path(data_dir,"riskset_matched_overall_assignments.rds")
  saveRDS(fixture_assignments(d),ap)
  mm$input_md5 <- unname(tools::md5sum(matched_path))
  mm$assignment_md5 <- unname(tools::md5sum(ap))
  saveRDS(mm,file.path(result_dir,"glme_spline_metadata.rds"))
  provenance <- list(eligibility_version = "no_fi_requirement_v1", fi_requirement = "none",
    cohort_entry_rule = "first_participated_analytic_cycle_return", control_entry_on_or_before_index = TRUE,
    assignment_md5=unname(tools::md5sum(ap)), output_md5=unname(tools::md5sum(matched_path)),input_md5=unname(tools::md5sum(panel_path)),
    target_cycles=c("88","92","96","00","04","08","12","16","20"))
  saveRDS(provenance,file.path(md_dir,paste0(stem,"_run_metadata.rds")))
  saveRDS(data.frame(index_age_mean=70,index_age_sd=8),file.path(md_dir,paste0(stem,"_scaling_metadata.rds")))
  gate_path <- file.path(md_dir,paste0(stem,"_gate_g4.csv"))
  write.csv(data.frame(gate_pass=TRUE),gate_path,row.names=FALSE)
  stopifnot(nrow(jm_validate_inputs(base)$matched)==nrow(d))
  write.csv(data.frame(gate_pass=FALSE),gate_path,row.names=FALSE)
  expect_error(jm_validate_inputs(base),"Gate G4 failed")
  write.csv(data.frame(gate_pass=TRUE),gate_path,row.names=FALSE)
  broken <- mm; broken$n_obs <- broken$n_obs+1L
  saveRDS(broken,file.path(result_dir,"glme_spline_metadata.rds"))
  expect_error(jm_validate_inputs(base),"does not match the current primary")
  saveRDS(mm,file.path(result_dir,"glme_spline_metadata.rds"))
  changed <- panel; changed$dtdth[2] <- 1400
  saveRDS(changed,panel_path)
  expect_error(jm_validate_inputs(base),"Canonical panel differs")
  saveRDS(panel,panel_path)
  changed <- d; changed$fi_score_nocancer[1] <- 0.8
  saveRDS(changed,matched_path)
  expect_error(jm_validate_inputs(base),"Matched input hash differs")
})
check("association checkpoint survives next-fit failure without any real fitting", {
  test_root <- tempfile("checkpoint_",file.path(project_dir,"Codex","jm_m2_revision_2026-09-29")); dir.create(test_root)
  original_validate <- jm_validate_inputs; original_components <- jm_fit_components
  original_link <- jm_fit_linked; original_gate <- jm_require_resolved_dates
  jm_require_resolved_dates <- function() invisible(TRUE) # Test-only; both fitters below are mocks.
  jm_validate_inputs <- function(...) list(matched=d,panel=panel,assignments=fixture_assignments(d),metadata=m,provenance=list(hashes=c(fixture="synthetic")))
  jm_fit_components <- function(p,formula,run_dir) {
    bn <- colnames(model.matrix(formula,p$long))
    fixed <- setNames(rep(0.01,length(bn)),bn)
    vv <- diag(0.01,length(bn)); dimnames(vv) <- list(bn,bn)
    mock_lme <- structure(list(coefficients=list(fixed=fixed),varFix=vv),class="lme")
    mock_cox <- structure(list(coefficients=c("GroupCancer Case"=0.1,"index_age_z"=0.2)),class="coxph")
    # These are minimal test doubles, not fitted component models.
    cp <- list(longitudinal=mock_lme,survival=mock_cox)
    jm_save(cp,file.path(run_dir,"components.rds"))
    cp
  }
  calls <- 0L
  jm_fit_linked <- function(components,p,config,functional_form) {
    calls <<- calls+1L
    run <- list.dirs(file.path(test_root,"Results","cancer","data","8.1_joint_model_overall","full"),recursive=FALSE)
    stopifnot(file.exists(file.path(run,"components.rds")))
    if(calls==2L) {
      stopifnot(file.exists(file.path(run,"value","joint_model.rds")),
                file.exists(file.path(run,"value","mcmc_diagnostics.csv")))
      stop("SYNTHETIC_SECOND_ASSOCIATION_FAILURE")
    }
    beta_names <- names(nlme::fixef(components$longitudinal))
    gamma_names <- names(coef(components$survival))
    mock_joint(beta_names,gamma_names,"value")
  }
  expect_error(jm_run(test_root,"full"),"One or more joint models failed")
  run <- list.dirs(file.path(test_root,"Results","cancer","data","8.1_joint_model_overall","full"),recursive=FALSE)
  stopifnot(calls==2L,file.exists(file.path(run,"value","joint_model.rds")),
            readRDS(file.path(run,"value_slope","status.rds"))$error=="SYNTHETIC_SECOND_ASSOCIATION_FAILURE")
  jm_validate_inputs <- original_validate; jm_fit_components <- original_components; jm_fit_linked <- original_link; jm_require_resolved_dates <- original_gate
})
check("M2 scaling and reference profile applied once", {
  p <- jm_prepare(d,panel); raw <- p
  z <- jm_scale_data(p,m)
  stopifnot(identical(z$long$base_calor_raw,raw$long$base_calor),is.factor(z$long$base_pckgr))
  expect_error(jm_scale_data(z,m),"already scaled")
  design <- jm_design(z,jm_formula(m),list(window=c(-8,8),theta_supported=TRUE))
  for(v in jm_diet_covars) {
    j <- match(v,z$scaling$column)
    stopifnot(abs(design$raw_reference[[v]]-mean(raw$long[[v]]))<1e-12,
      abs(design$reference[[v]]-(mean(raw$long[[v]])-z$scaling$center[j])/z$scaling$scale[j])<1e-12)
    long_base <- z$long[!duplicated(z$long$id),c("id",v)]
    stopifnot(identical(long_base[[v]],z$survival[[v]]))
  }
  miss <- d;miss$base_calor[miss$id=="A" & miss$role=="Case"]<-NA_real_
  pp<-jm_prepare(miss,panel)
  stopifnot(!"A" %in% as.character(pp$long$id),pp$retained$role[pp$retained$id=="A"]=="Case")
  stopifnot(any(pp$filter_log$stage=="M2_complete_case"))
  bad <- m;bad$schema_version<-2L;expect_error(jm_validate_metadata(bad),"current overall M2")
  bad <- m;bad$model_specification$M2_full_spline$covars<-jm_primary_covars
  expect_error(jm_validate_metadata(bad),"current overall M2")
  bad <- m;bad$model_scaling$M2_full_spline$scale[4]<-0
  expect_error(jm_validate_metadata(bad),"Invalid M2")
})
check("complete date issue classification and overlap without date changes", {
  a <- jm_assignment_ledger(d)$retained
  a$index_date[a$id=="C"] <- jm_admin_date
  a$cancer_dateca[a$id=="C"] <- jm_admin_date
  p <- panel;p$dtdth[p$id=="A"]<-1200;p$dtdth[p$id=="D"]<-1199
  p <- rbind(p,data.frame(id="B",dtdth=1300),data.frame(id="B",dtdth=1350))
  original <- a
  au <- jm_date_audit(a,p)
  stopifnot(identical(a,original),au$records$same_month_death[au$records$id=="A"],
    au$records$conflicting_death_dates[au$records$id=="B"],
    au$records$death_before_index[au$records$id=="D"],
    au$records$index_on_after_admin[au$records$id=="C"],
    au$records$same_month_death[au$records$id=="C"],
    sum(au$overlap$participants)==nrow(a))
  p$dtdth[1]<- -1
  stopifnot(jm_date_audit(a,p)$records$invalid_death_source[au$records$id=="A"])
})
check("all free posterior blocks gated; diagonal zero constraints verified", {
  set.seed(14)
  fake <- mock_joint(c("(Intercept)","GroupCancer Case"),"GroupCancer Case",diagonal=TRUE)
  args<-list(beta_names=c("(Intercept)","GroupCancer Case"),gamma_names="GroupCancer Case",structure="value_slope",expected_chains=3)
  dr<-do.call(jm_posterior,c(list(fit=fake),args))
  stopifnot(all(jm_diagnostics(dr,"full")$pass),nrow(attr(dr,"constrained_parameters"))==1,
    !"D[2, 1]" %in% dimnames(dr$D)[[3]])
  for(block in c("sigmas","D","bs_gammas","tau_bs_gammas","betas","gammas")) {
    bad<-dr;bad[[block]][,2,1]<-bad[[block]][,2,1]+5
    stopifnot(!all(jm_diagnostics(bad,"full")$pass))
  }
  bad<-fake;bad$mcmc$D[[1]][,2]<-0.01
  expect_error(do.call(jm_posterior,c(list(fit=bad),args)),"Constrained covariance")
  bad<-fake;bad$mcmc$bs_gammas<-NULL
  expect_error(do.call(jm_posterior,c(list(fit=bad),args)),"chain count")
  stopifnot(nrow(jm_fit_criteria(fake,"SYNTHETIC"))==3,
    jm_fit_details(fake)$alignment_model_id==jm_alignment_model)
  bad<-fake;bad$control$Bsplines_degree<-3L
  expect_error(jm_fit_details(bad),"Realized baseline hazard")
})
check("failed baseline-hazard mixing suppresses all full posterior inference", {
  root<-tempfile("failed_diagnostics_",file.path(project_dir,"Codex","jm_m2_revision_2026-09-29"));dir.create(root)
  original_validate<-jm_validate_inputs;original_components<-jm_fit_components
  original_link<-jm_fit_linked;original_gate<-jm_require_resolved_dates
  jm_require_resolved_dates<-function() invisible(TRUE) # Test only, both fitters mocked.
  jm_validate_inputs<-function(...)list(matched=d,panel=panel,assignments=fixture_assignments(d),metadata=m,provenance=list(hashes=c(fixture="synthetic")))
  jm_fit_components<-function(p,formula,run_dir) {
    bn<-colnames(model.matrix(formula,p$long));b<-setNames(rep(.01,length(bn)),bn)
    V<-diag(.01,length(bn));dimnames(V)<-list(bn,bn)
    list(longitudinal=structure(list(coefficients=list(fixed=b),varFix=V),class="lme"),
      survival=structure(list(coefficients=c("GroupCancer Case"=.1)),class="coxph"))
  }
  jm_fit_linked<-function(components,p,config,functional_form) {
    structure<-if(grepl("slope",paste(deparse(functional_form),collapse="")))"value_slope" else "value"
    z<-mock_joint(names(nlme::fixef(components$longitudinal)),"GroupCancer Case",structure)
    z$mcmc$bs_gammas[[2]][,1]<-z$mcmc$bs_gammas[[2]][,1]+5
    z
  }
  expect_error(jm_run(root,"full"),"One or more joint models failed")
  run<-list.dirs(file.path(root,"Results","cancer","data","8.1_joint_model_overall","full"),recursive=FALSE)
  for(nm in c("value","value_slope"))stopifnot(file.exists(file.path(run,nm,"joint_model.rds")),
    !file.exists(file.path(run,nm,"posterior_summaries.rds")),!file.exists(file.path(run,nm,"joint_curves.csv")),
    !file.exists(file.path(run,nm,"marginal_fit_statistics.csv")))
  jm_validate_inputs<-original_validate;jm_fit_components<-original_components
  jm_fit_linked<-original_link;jm_require_resolved_dates<-original_gate
})
check("source-safe runner and original factor-level coding", {
  original_run<-jm_run;jm_run<-function(...)stop("RUNNER EXECUTED WHILE SOURCING")
  source(file.path(project_dir,"Code","2_data_analysis","8.1_joint_model_overall.R"))
  # Wrapper sources helpers but never calls the workflow when sourced.
  z<-d;z$base_pckgr<-factor(z$base_pckgr,levels=c("b","a"))
  pp<-production_prepare(z,panel,fixture_assignments(z))
  stopifnot(identical(levels(pp$long$base_pckgr),c("b","a")))
  jm_run<-original_run
})
cat("\n",length(passed),"model-free contract groups passed. No model fitted.\n")

