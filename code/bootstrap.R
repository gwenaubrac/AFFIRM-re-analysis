## ---------------------------
##
## Program: AFFIRM Re-Analysis Bootstrapping
## Author: Gwen Aubrac
##
## Date Created: 2025-07-09
##
## ---------------------------
##
## Notes: 5,000 bootstrap iterations
## Wrapped models in tryCatch to resolve issues when running bootstrap:
## - removed smoking due to positivity violations in samples
## - removed any pre-term due to positivity violations in samples
## - used stabilized IOSW weights to fit weighted IPTW for stability
##
## ---------------------------

#### 1. LOAD PACKAGES ####

options(scipen = 999)

library(boot)
library(dplyr)
library(magrittr)
library(data.table)
library(writexl)
library(glmtoolbox)
library(splines)
library(ggplot2)

#### 2. READ DATA ####

# set working directory to analysis, among:
# main_results
# target_multi_site
# aspirin
# no_aspirin
# restrict_GA_10w

setwd("Z:/RPLATT/GWEN-OHRI/main_results")
data <- readRDS("../data_clean.R")

data %<>% # for sensitivity where target multi-site, switch this indicator
  mutate(S = if_else(single == 1, 0, 1))

# for subgroup analyses by aspirin, add this filter
# and also remove aspirin cov from models
# data %<>%
#   filter(aspirin == 0)

setDT(data)

# keep track of failed iterations:
failures <<- 0
bootstrap_samples <<- list()
OM_values <<- numeric() 
DR1_values <<- numeric() 
DR2_values <<- numeric() 
GEE_sw_values <<- numeric()

#### 3. BOOSTRAP FUNCTION ####

bs <- function(data, indices) {
  
  d <- data[indices,]
  d[, boot_id := .I]
  results <- c()
  
  #### a) CRUDE ATE PRIOR TO RESTRICTIONS ####
  
  # all patients:
  all_crude_prior_1 <- sum(d$outcome[d$trt == 1])/sum(d$trt == 1)
  all_crude_prior_0 <- sum(d$outcome[d$trt == 0])/sum(d$trt == 0)
  all_crude_prior <- all_crude_prior_1 - all_crude_prior_0
  
  # S==1 patients:
  S1crude_prior_1 <- sum(d$outcome[d$trt == 1 & d$S == 1])/sum(d$trt == 1 & d$S == 1)
  S1crude_prior_0 <- sum(d$outcome[d$trt == 0 & d$S == 1])/sum(d$trt == 0 & d$S == 1)
  S1crude_prior <- S1crude_prior_1 - S1crude_prior_0
  
  # S == 0 patients:
  S0crude_prior_1 <- sum(d$outcome[d$trt == 1 & d$S == 0])/sum(d$trt == 1 & d$S == 0)
  S0crude_prior_0 <- sum(d$outcome[d$trt == 0 & d$S == 0])/sum(d$trt == 0 & d$S == 0)
  S0crude_prior <- S0crude_prior_1 - S0crude_prior_0
  
  #### b) CRUDE ATE AFTER APPLY RESTRICTIONS ####
  
  to_exclude <- d %>%
    filter(
      hx_any_loss != 0 |
        gravida_cat != 2 |
        live_births_cat != 1 |
        any_vte == 1 |
        #GA_at_start_cat == "10 to <16 weeks" | # UNCOMMENT FOR SENSITIVITY ANALYSIS
        GA_at_start_cat == "16 to <20 weeks" |
        GA_at_start_cat == ">=20 weeks"
    )

  d %<>%
    filter(!boot_id %in% to_exclude$boot_id)
  
  # crude in all
  all_crude_1 <- sum(d$outcome[d$trt == 1])/sum(d$trt == 1)
  all_crude_0 <- sum(d$outcome[d$trt == 0])/sum(d$trt == 0)
  all_crude <- all_crude_1 - all_crude_0
  
  # crude in S == 1
  S1crude_1 <- sum(d$outcome[d$trt == 1 & d$S == 1])/sum(d$trt == 1 & d$S == 1)
  S1crude_0 <- sum(d$outcome[d$trt == 0 & d$S == 1])/sum(d$trt == 0 & d$S == 1)
  S1crude <- S1crude_1 - S1crude_0
  
  # crude in S == 0
  S0crude_1 <- sum(d$outcome[d$trt == 1 & d$S == 0])/sum(d$trt == 1 & d$S == 0)
  S0crude_0 <- sum(d$outcome[d$trt == 0 & d$S == 0])/sum(d$trt == 0 & d$S == 0)
  S0crude <- S0crude_1 - S0crude_0
  
  
  #### c) OUTCOME MODEL ATE ####
  
  # use info from S==1 to estimate outcome in S==0
  S1data_A1 <- subset(d, S == 1 & trt == 1)
  S1data_A0 <- subset(d, S == 1 & trt == 0)
  
  OM1mod <- tryCatch({
    glm(
      outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin, 
      data = S1data_A1,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with OM1mod \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })

  OM0mod <- tryCatch({
    glm(formula =       
          outcome ~
          age__years_ +
          any_pe +
          any_preterm +
          any_sga +
          bmi_above_30 +
          smoking_bin,
          aspirin,
        data = S1data_A0,
        family = binomial (link = "logit"))
  }, error = function(e) {
    cat("Error with OM0mod \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(OM1mod) | is.null(OM0mod)) {
    
    S1_OM_1 <- NA
    S1_OM_0 <- NA
    S1_OM <- NA
    
  } else {
    
    d$p1 <- predict(OM1mod, newdata = d, type = "response")
    d$p0 <- predict(OM0mod, newdata = d, type = "response")
    S0_data <- subset(d, S == 0)
    S1_OM_1 <- mean(S0_data$p1)
    S1_OM_0 <- mean(S0_data$p0)
    S1_OM <- mean(S0_data$p1) - mean(S0_data$p0)  
    
    }
  
  
  # in all patients
  # (use data from all patients to estimate outcome)
  data_A1 <- d %>% filter(trt == 1)
  data_A0 <- d %>% filter(trt == 0)
  
  all_OM1mod <- tryCatch({
    glm(
      outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin, 
        aspirin, 
      data = data_A1,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with all OM1mod \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  all_OM0mod <- tryCatch({
    glm(formula =       
          outcome ~
          age__years_ +
          any_pe +
          any_preterm +
          any_sga +
          bmi_above_30 +
          smoking_bin,
          aspirin,
        data = data_A0,
        family = binomial (link = "logit"))
  }, error = function(e) {
    cat("Error with all OM0mod \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(all_OM1mod) | is.null(all_OM0mod)) {
    
    all_OM_1 <- NA
    all_OM_0 <- NA
    all_OM <- NA
    
  } else {
    
    d$all_p1 <- predict(all_OM1mod, newdata = d, type = "response")
    d$all_p0 <- predict(all_OM0mod, newdata = d, type = "response")
    all_OM_1 <- mean(d$all_p1)
    all_OM_0 <- mean(d$all_p0)
    all_OM <- mean(d$all_p1) - mean(d$all_p0)  
    
  }
  
  ####  d) ESTIMATE IOSW WEIGHTS ####
  
  w_reg <- tryCatch({
    glm(
      S ~ 
        ns(age__years_, df = 3) + 
        any_pe + 
        any_preterm +
        any_sga + 
        bmi_above_30 +
        smoking_bin, 
        aspirin +
        aspirin * ns(age__years_, df = 2),
      family = binomial(),
      data = d
    )
  }, error = function(e) {
    cat("Error with w_reg -- could not compute IOSW so return NA for all \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(w_reg)) return(NA)
  

  d$ps <- predict(w_reg, newdata = d, type = "response")
  d[, iosw := fifelse(S == 1, (1 - ps) / ps, 0)] # only S1 patients get a weight
  d[, iosw_all := fifelse(S == 1, (1 - ps) / ps, 1)] # all patients get a weight (w=1 for S0)
  
  marg_s_model <- tryCatch({
    glm(S ~ 1, family = binomial(), data = d)
  }, error = function(e) {
    cat("Error with marg s model \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  marg_ps <- predict(marg_s_model, newdata = d, type = "response")
  odds_s <- marg_ps/(1-marg_ps)
  d[, siosw := fifelse(S == 1, odds_s*iosw, 1)]
  d[, siosw_all := fifelse(S == 1, odds_s*iosw_all, 1)]
  
  
  #### e) ESTIMATE IPTW ####

  wmulti_d <- d %>% filter(S == 1)
  
  trt_model <- tryCatch({
    glm(
      trt ~ 
        age__years_ +
        any_pe + 
        any_preterm +
        any_sga + 
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = wmulti_d,
      weights = siosw, # use SIOSW rather than IOSW for improved stability
      family = binomial(link = "logit")
    )
  }, error = function(e) {
    cat("Error with trt model -- could not compute IPTW so return NA for all \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(trt_model)) return(NA)
  
  wmulti_d$pa <- predict(trt_model, type = "response")
  wmulti_d[, iptw := fifelse(trt == 1, 1/pa, 1/(1-pa))]
  
  marg_trt_model <- tryCatch({
    glm(trt ~ 1, data = wmulti_d, family = binomial(link = "logit"))
  }, error = function(e) {
    cat("Error with marg trt model \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  wmulti_d$marg_pa <- predict(marg_trt_model, type = "response")
  wmulti_d[, siptw := fifelse(trt == 1, iptw * marg_pa, iptw * (1-marg_pa))]

  d[, `:=`(iptw = NA_real_, siptw = NA_real_)]
  
  d[S == 1, `:=`(
    iptw = wmulti_d$iptw,
    siptw = wmulti_d$siptw
  )]
  
  
  single_d <- d %>% 
    filter(S == 0)
  
  trt_model_single <- tryCatch({
    glm(
      trt ~ 
        age__years_ +
        any_pe + 
        any_preterm +
        any_sga + 
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = single_d,
      family = binomial(link = "logit")
    )
  }, error = function(e) {
    cat("Error with trt model single -- could not compute IPTW so return NA for all \n")
    failures <- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  if (is.null(trt_model_single)) return(NA)
  
  single_d$pa <- predict(trt_model_single, type = "response")
  single_d[, iptw := fifelse(trt == 1, 1/pa, 1/(1-pa))]
  
  marg_trt_model <- tryCatch({
    glm(trt ~ 1, data = single_d, family = binomial())
  }, error = function(e) {
    cat("Error with marg trt model \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  single_d$marg_pa <- predict(marg_trt_model, type = "response")
  single_d[, siptw := fifelse(trt == 1, iptw * marg_pa, iptw * (1-marg_pa))]
  
  d[S == 0, `:=`(
    iptw = single_d$iptw,
    siptw = single_d$siptw
  )]
  
  d[, `:=`(
    w = iosw * iptw,
    sw = siosw * siptw,
    w_all = iosw_all * iptw,
    sw_all = siosw_all * siptw
  )]
  
  
  
  #### f) IOSW-WEIGHTED ATE ####
  
  A <- d$trt
  S <- d$S
  w <- d$w
  sw <- d$sw
  Y <- d$outcome
  
  # in weighted multi-site patients:
  S1_IOW1_1 <- (sum(w)^-1) * sum(A*S*w*Y)
  S1_IOW1_0 <- (sum(w)^-1) * sum((1-A)*S*w*Y)
  S1_IOW1 <- S1_IOW1_1 - S1_IOW1_0
  
  S1_IOW2_1 <- (sum(sw)^-1) * sum(A*S*sw*Y)
  S1_IOW2_0 <- (sum(sw)^-1) * sum((1-A)*S*sw*Y)
  S1_IOW2 <- S1_IOW2_1 - S1_IOW2_0
  
  # in all patients:
  all_IOW1_1 <- (sum(w_all)^-1) * sum(A*w_all*Y)
  all_IOW1_0 <- (sum(w_all)^-1) * sum((1-A)*w_all*Y)
  all_IOW1 <- all_IOW1_1 - all_IOW1_0
  
  all_IOW2_1 <- (sum(sw_all)^-1) * sum(A*sw_all*Y)
  all_IOW2_0 <- (sum(sw_all)^-1) * sum((1-A)*sw_all*Y)
  all_IOW2 <- all_IOW2_1 - all_IOW2_0
  
  
  
  #### g) DOUBLY ROBUST ESTIMATOR ####
  
  # in weighted S==1 patients
  if (is.null(OM0mod) | is.null(OM1mod)) {
    
    S1_DR1_1 <- NA
    S1_DR1_0 <- NA
    S1_DR1 <- NA
    S1_DR2_1 <- NA
    S1_DR2_0 <- NA
    S1_DR2 <- NA
    
  } else {
    
    p1 <- d$p1
    p0 <- d$p0
    
    S1_DR1_1 <- (sum(w)^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
    S1_DR1_0 <- (sum(w)^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)
    S1_DR1 <- S1_DR1_1 - S1_DR1_0
    
    S1_DR2_1 <- (sum(sw)^-1) * sum(S*A*sw*(Y-p1) + (1-S) * p1)
    S1_DR2_0 <- (sum(sw)^-1) * sum(S*(1-A)*sw*(Y-p0) + (1-S) * p0)
    S1_DR2 <- S1_DR2_1 - S1_DR2_0
    
  }
  
  S1data_A1 <- subset(d, S == 1 & trt == 1)
  S1data_A0 <- subset(d, S == 1 & trt == 0)
  
  DR1mod_w <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = S1data_A1,
      weights = S1data_A1$w,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with DR1mod w \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    failures <<- failures + 1
    return(NULL)
  })

  DR0mod_w <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = S1data_A0,
      weights = S1data_A0$w,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with DR0mod w \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(DR1mod_w) | is.null(DR0mod_w)) {
    
    S1_DR3_1 <- NA
    S1_DR3_0 <- NA
    S1_DR3 <- NA
    
  } else {
    
    p1_w <- predict(DR1mod_w, newdata = S0data, type = "response")
    p0_w <- predict(DR0mod_w, newdata = S0data, type = "response")
    S1_DR3_1 <- mean(p1_w)
    S1_DR3_0 <- mean(p0_w)
    S1_DR3 <- S1_DR3_1 - S1_DR3_0
    
    }
  
  DR1mod_sw <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = S1data_A1,
      weights = S1data_A1$sw,
      family = binomial (link = "logit")
    )
  },
  error = function(e) {
    cat("Error with DR1mod sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  DR0mod_sw <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = S1data_A0,
      weights = S1data_A0$sw,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with DR0mod sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null (DR1mod_sw) | is.null(DR0mod_sw)) {
    
    S1_DR4_1 <- NA
    S1_DR4_0 <- NA
    S1_DR4 <- NA
   
  } else {
    p1_sw <- predict(DR1mod_sw, newdata = S0data, type = "response")
    p0_sw <- predict(DR0mod_sw, newdata = S0data, type = "response")
    
    S1_DR4_1 <- mean(p1_sw)
    S1_DR4_0 <- mean(p0_sw)
    S1_DR4 <- S1_DR4_1 - S1_DR4_0
  }
  
  # in all patients
  if (is.null(all_OM0mod) | is.null(all_OM1mod)) {
    
    all_DR1_1 <- NA
    all_DR1_0 <- NA
    all_DR1 <- NA
    all_DR2_1 <- NA
    all_DR2_0 <- NA
    all_DR2 <- NA
    
  } else {
    
    all_p1 <- d$all_p1
    all_p0 <- d$all_p0
    
    all_DR1_1 <- (sum(w_all)^-1) * sum(A*w_all*(Y-all_p1) + all_p1)
    all_DR1_0 <- (sum(w_all)^-1) * sum((1-A)*w_all*(Y-all_p0) + all_p0)
    all_DR1 <- all_DR1_1 - all_DR1_0
    
    all_DR2_1 <- (sum(sw_all)^-1) * sum(A*sw_all*(Y-all_p1) + all_p1)
    all_DR2_0 <- (sum(sw_all)^-1) * sum((1-A)*sw_all*(Y-all_p0) + all_p0)
    all_DR2 <- all_DR2_1 - all_DR2_0
    
  }
  
  all_DR1mod_w <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = data_A1,
      weights = data_A1$w_all,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with all DR1mod w \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    failures <<- failures + 1
  })
  
  all_DR0mod_w <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = data_A0,
      weights = data_A0$w_all,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with all DR0mod w \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(all_DR1mod_w) | is.null(all_DR0mod_w)) {
    
    all_DR3_1 <- NA
    all_DR3_0 <- NA
    all_DR3 <- NA
    
  } else {
    
    all_p1_w <- predict(all_DR1mod_w, newdata = d, type = "response")
    all_p0_w <- predict(all_DR0mod_w, newdata = d, type = "response")
    all_DR3_1 <- mean(all_p1_w)
    all_DR3_0 <- mean(all_p0_w)
    all_DR3 <- all_DR3_1 - all_DR3_0
    
  }
  
  all_DR1mod_sw <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = data_A1,
      weights = data_A1$sw_all,
      family = binomial (link = "logit")
    )
  },
  error = function(e) {
    cat("Error with all DR1mod sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  all_DR0mod_sw <- tryCatch({
    glm(
      formula =       
        outcome ~
        age__years_ +
        any_pe +
        any_preterm +
        any_sga +
        bmi_above_30 +
        smoking_bin,
        aspirin,
      data = data_A0,
      weights = data_A0$sw_all,
      family = binomial (link = "logit")
    )
  }, error = function(e) {
    cat("Error with all DR0mod sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null (all_DR1mod_sw) | is.null(all_DR0mod_sw)) {
    
    all_DR4_1 <- NA
    all_DR4_0 <- NA
    all_DR4 <- NA
    
  } else {
    
    all_p1_sw <- predict(all_DR1mod_sw, newdata = d, type = "response")
    all_p0_sw <- predict(all_DR0mod_sw, newdata = d, type = "response")
    
    all_DR4_1 <- mean(all_p1_sw)
    all_DR4_0 <- mean(all_p0_sw)
    all_DR4 <- all_DR4_1 - all_DR4_0
    
  }
  
  
  #### h) GEE ####
  
  S1_data <- d %>% 
    filter(S == 1)
  
  S0_data <- d %>% 
    filter(S == 0)
  
  ## crude GEE ##
  
  # S1 patients
  S1gee_crude_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      id = unique_id,
      data = S1_data,
      family = binomial("logit"),
      corstr = "exchangeable"
    )
  }, error = function(e) {
    cat("Error with S1 GEE crude \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    failures <<- failures + 1
    return(NULL)
  })
  
  if (is.null(S1gee_crude_model)) {
    
    S1_GEE_crude_1 <- NA
    S1_GEE_crude_0 <- NA
    S1_GEE_crude <- NA
    
  } else {
    
    treated <- S1data
    treated$trt <- 1
    
    untreated <- S1data
    untreated$trt <- 0
    
    treated$pred <- predict(S1gee_crude_model, newdata = S1treated, type = "response")
    untreated$pred <- predict(S1gee_crude_model, newdata = S1untreated, type = "response")
    
    S1_GEE_crude_1 <- mean(S1treated$pred)
    S1_GEE_crude_0 <- mean(S1untreated$pred)
    S1_GEE_crude <- S1_GEE_crude_1 - S1_GEE_crude_0
    
    rm(treated, untreated)
    
  }
  
  # S0 patients
  S0gee_crude_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      id = unique_id,
      data = S0_data,
      family = binomial("logit"),
      corstr = "exchangeable"
    )
  }, error = function(e) {
    cat("Error with S0 GEE crude \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    failures <<- failures + 1
    return(NULL)
  })
  if (is.null(S0gee_crude_model)) {
    
    S0_GEE_crude_1 <- NA
    S0_GEE_crude_0 <- NA
    S0_GEE_crude <- NA
    
  } else {
    
    treated <- S0_data
    treated$trt <- 1
    
    untreated <- S0_data
    untreated$trt <- 0
    
    treated$pred <- predict(S0gee_crude_model, newdata = treated, type = "response")
    untreated$pred <- predict(S0gee_crude_model, newdata = untreated, type = "response")
    
    S0_GEE_crude_1 <- mean(treated$pred)
    S0_GEE_crude_0 <- mean(untreated$pred)
    S0_GEE_crude <- S0_GEE_crude_1 - S0_GEE_crude_0
    
    rm(treated, untreated)
    
  }
  
  # all patients
  gee_crude_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      id = unique_id,
      data = d,
      family = binomial("logit"),
      corstr = "exchangeable"
    )
  }, error = function(e) {
    cat("Error with GEE crude \n")
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    failures <<- failures + 1
    return(NULL)
  })
  if (is.null(gee_crude_model)) {
    
    all_GEE_crude_1 <- NA
    all_GEE_crude_0 <- NA
    all_GEE_crude <- NA
    
  } else {
    
    treated <- d
    treated$trt <- 1
    
    untreated <- d
    untreated$trt <- 0
    
    treated$pred <- predict(gee_crude_model, newdata = treated, type = "response")
    untreated$pred <- predict(gee_crude_model, newdata = untreated, type = "response")
    
    all_GEE_crude_1 <- mean(treated$pred)
    all_GEE_crude_0 <- mean(untreated$pred)
    all_GEE_crude <- all_GEE_crude_1 - all_GEE_crude_0
    
    rm(treated, untreated)
    
  }
  
  
  
  ## unstabilized weighted GEE
  
  # S=1 patients
  gee_w_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = S1_data,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = S1_data$w
    )
  }, error = function(e) {
    cat("Error with GEE w \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(gee_w_model)) {
    
    S1_GEE_IOW1_1 <- NA
    S1_GEE_IOW1_0 <- NA
    S1_GEE_IOW1 <- NA
    
  } else {
    
    treated <- S1_data
    treated$trt <- 1
    
    untreated <- S1_data
    untreated$trt <- 0
    
    treated$pred <- predict(gee_w_model, newdata = treated, type = "response")
    untreated$pred <- predict(gee_w_model, newdata = untreated, type = "response")
    
    S1_GEE_IOW1_1 <- mean(treated$pred)
    S1_GEE_IOW1_0 <- mean(untreated$pred)
    S1_GEE_IOW1 <- S1_GEE_IOW1_1 - S1_GEE_IOW1_0
    
    rm(treated, untreated)
  }
    
  
  # S=0 patients
  gee_w_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = S0_data,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = S1_data$w
    )
  }, error = function(e) {
    cat("Error with GEE w \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(gee_w_model)) {
    
    S0_GEE_IOW1_1 <- NA
    S0_GEE_IOW1_0 <- NA
    S0_GEE_IOW1 <- NA
    
  } else {
    
    treated <- S0_data
    treated$trt <- 1
    
    untreated <- S0_data
    untreated$trt <- 0
    
    treated$pred <- predict(gee_w_model, newdata = treated, type = "response")
    untreated$pred <- predict(gee_w_model, newdata = untreated, type = "response")
    
    S0_GEE_IOW1_1 <- mean(treated$pred)
    S0_GEE_IOW1_0 <- mean(untreated$pred)
    S0_GEE_IOW1 <- S0_GEE_IOW1_1 - S0_GEE_IOW1_0
    
    rm(treated, untreated)
  }
  
  
  # all patients
  all_gee_w_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = d,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = d$w_all
    )
  }, error = function(e) {
    cat("Error with all GEE w \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(all_gee_w_model)) {
    
    all_GEE_IOW1_1 <- NA
    all_GEE_IOW1_0 <- NA
    all_GEE_IOW1 <- NA
    
  } else {
    
    treated <- d
    treated$trt <- 1
    
    untreated <- d
    untreated$trt <- 0
    
    treated$pred <- predict(all_gee_w_model, newdata = treated, type = "response")
    untreated$pred <- predict(all_gee_w_model, newdata = untreated, type = "response")
    
    all_GEE_IOW1_1 <- mean(treated$pred)
    all_GEE_IOW1_0 <- mean(untreated$pred)
    all_GEE_IOW1 <- all_GEE_IOW1_1 - all_GEE_IOW1_0
    
    rm(untreated, treated)
  }
  
  ## stabilized weighted GEE
  
  # S1 patients
  gee_sw_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = S1_data,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = S1_data$sw
    )
  }, error = function(e) {
    cat("Error with GEE sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(gee_sw_model)) {
    
    S1_GEE_IOW2_1 <- NA
    S1_GEE_IOW2_0 <- NA
    S1_GEE_IOW2 <- NA
    
  } else {
    
    treated <- S1_data
    treated$trt <- 1
    
    untreated <- S1_data
    untreated$trt <- 0
    
    treated$pred <- predict(gee_sw_model, newdata = treated, type = "response")
    untreated$pred <- predict(gee_sw_model, newdata = untreated, type = "response")
    
    S1_GEE_IOW2_1 <- mean(treated$pred)
    S1_GEE_IOW2_0 <- mean(untreated$pred)
    S1_GEE_IOW2 <- S1_GEE_IOW2_1 - S1_GEE_IOW2_0
    rm(treated, untreated)
    
  }
  
  # S0 patients
  gee_sw_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = S0_data,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = S0_data$sw
    )
  }, error = function(e) {
    cat("Error with GEE sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(gee_sw_model)) {
    
    S0_GEE_IOW2_1 <- NA
    S0_GEE_IOW2_0 <- NA
    S0_GEE_IOW2 <- NA
    
  } else {
    
    treated <- S0_data
    treated$trt <- 1
    
    untreated <- S0_data
    untreated$trt <- 0
    
    treated$pred <- predict(gee_sw_model, newdata = treated, type = "response")
    untreated$pred <- predict(gee_sw_model, newdata = untreated, type = "response")
    
    S0_GEE_IOW2_1 <- mean(treated$pred)
    S0_GEE_IOW2_0 <- mean(untreated$pred)
    S0_GEE_IOW2 <- S0_GEE_IOW2_1 - S0_GEE_IOW2_0
    rm(treated, untreated)
    
  }
  
  
  
  # in all patients
  all_gee_sw_model <- tryCatch({
    glmgee(
      outcome ~ trt,
      data = d,
      id = unique_id,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      weights = d$all_sw
    )
  }, error = function(e) {
    cat("Error with all GEE sw \n")
    failures <<- failures + 1
    saveRDS(d, file = "bad_bootstrap_sample.rds")
    return(NULL)
  })
  
  if (is.null(all_gee_sw_model)) {
    
    all_GEE_IOW2_1 <- NA
    all_GEE_IOW2_0 <- NA
    all_GEE_IOW2 <- NA
    
  } else {
    
    treated <- d
    treated$trt <- 1
    
    untreated <- d
    untreated$trt <- 0
    
    treated$pred <- predict(all_gee_sw_model, newdata = treated, type = "response")
    untreated$pred <- predict(all_gee_sw_model, newdata = untreated, type = "response")
    
    all_GEE_IOW2_1 <- mean(treated$pred)
    all_GEE_IOW2_0 <- mean(untreated$pred)
    all_GEE_IOW2 <- all_GEE_IOW2_1 - all_GEE_IOW2_0
    
    rm(treated, untreated)
    
  }
  
  
  
  #### i) ADDITIONAL ANALYSES ####
  
  ## by subgroups of S=0
  
  S1_data <- d %>% 
    filter(S == 1)
  
  A <- S1data$trt
  S <- S1data$S
  Y <- S1data$outcome
  iosw <- S1data$iosw
  siosw <- S1data$siosw
  iptw <- S1data$iptw
  siptw <- S1data$siptw
  w <- S1data$w
  sw <- S1data$sw
  
  S1_iosw_1 <- (sum(iosw*A)^-1) * sum(A*iosw*Y)
  S1_iosw_0 <- (sum(iosw*(1-A))^-1) * sum((1-A)*iosw*Y)
  S1_iosw = S1_iosw_1 - S1_iosw_0
  
  S1_siosw_1 <- (sum(siosw*A)^-1) * sum(A*siosw*Y)
  S1_siosw_0 <- (sum(siosw*(1-A))^-1) * sum((1-A)*siosw*Y)
  S1_siosw = S1_siosw_1 - S1_siosw_0
  
  S1_iptw_1 <- (sum(iptw*A)^-1) * sum(A*iptw*Y)
  S1_iptw_0 <- (sum(iptw*(1-A))^-1) * sum((1-A)*iptw*Y)
  S1_iptw = S1_iptw_1 - S1_iptw_0
  
  S1_siptw_1 <- (sum(siptw*A)^-1) * sum(A*siptw*Y)
  S1_siptw_0 <- (sum(siptw*(1-A))^-1) * sum((1-A)*siptw*Y)
  S1_siptw = S1_siptw_1 - S1_siptw_0
  
  S1_w_1 <- (sum(w*A)^-1) * sum(A*w*Y)
  S1_w_0 <- (sum(w*(1-A))^-1) * sum((1-A)*w*Y)
  S1_w = S1_w_1 - S1_w_0
  
  S1_sw_1 <- (sum(sw*A)^-1) * sum(A*sw*Y)
  S1_sw_0 <- (sum(sw*(1-A))^-1) * sum((1-A)*sw*Y)
  S1_sw = S1_sw_1 - S1_sw_0
  
  
  
  ## by subgroups of S=0
  
  S0data <- d %>% 
    filter(S == 0)
  
  A <- S0data$trt
  S <- S0data$S
  Y <- S0data$outcome
  iosw <- S0data$iosw
  siosw <- S0data$siosw
  iptw <- S0data$iptw
  siptw <- S0data$siptw
  w <- S0data$w
  sw <- S0data$sw
  
  S0_iosw_1 <- (sum(iosw*A)^-1) * sum(A*iosw*Y)
  S0_iosw_0 <- (sum(iosw*(1-A))^-1) * sum((1-A)*iosw*Y)
  S0_iosw = S0_iosw_1 - S0_iosw_0
  
  S0_siosw_1 <- (sum(siosw*A)^-1) * sum(A*siosw*Y)
  S0_siosw_0 <- (sum(siosw*(1-A))^-1) * sum((1-A)*siosw*Y)
  S0_siosw = S0_siosw_1 - S0_siosw_0
  
  S0_iptw_1 <- (sum(iptw*A)^-1) * sum(A*iptw*Y)
  S0_iptw_0 <- (sum(iptw*(1-A))^-1) * sum((1-A)*iptw*Y)
  S0_iptw = S0_iptw_1 - S0_iptw_0
  
  S0_siptw_1 <- (sum(siptw*A)^-1) * sum(A*siptw*Y)
  S0_siptw_0 <- (sum(siptw*(1-A))^-1) * sum((1-A)*siptw*Y)
  S0_siptw = S0_siptw_1 - S0_siptw_0
  
  S0_w_1 <- (sum(w*A)^-1) * sum(A*w*Y)
  S0_w_0 <- (sum(w*(1-A))^-1) * sum((1-A)*w*Y)
  S0_w = S0_w_1 - S0_w_0
  
  S0_sw_1 <- (sum(sw*A)^-1) * sum(A*sw*Y)
  S0_sw_0 <- (sum(sw*(1-A))^-1) * sum((1-A)*sw*Y)
  S0_sw = S0_sw_1 - S0_sw_0
  

  # whole population
  
  A <- d$trt
  S <- d$S
  Y <- d$outcome
  iosw <- d$iosw
  siosw <- d$siosw
  iptw <- d$iptw
  siptw <- d$siptw
  w <- d$w_all
  sw <- d$sw_all
  
  all_iosw_1 <- (sum(iosw)^-1) * sum(A*S*iosw*Y)
  all_iosw_0 <- (sum(iosw)^-1) * sum((1-A)*S*iosw*Y)
  all_iosw = all_iosw_1 - all_iosw_0
  
  all_siosw_1 <- (sum(siosw)^-1) * sum(A*S*siosw*Y)
  all_siosw_0 <- (sum(siosw)^-1) * sum((1-A)*S*siosw*Y)
  all_siosw = all_siosw_1 - all_siosw_0
  
  all_iptw_1 <- (sum(iptw)^-1) * sum(A*S*iptw*Y)
  all_iptw_0 <- (sum(iptw)^-1) * sum((1-A)*S*iptw*Y)
  all_iptw = all_iptw_1 - all_iptw_0
  
  all_siptw_1 <- (sum(siptw)^-1) * sum(A*S*siptw*Y)
  all_siptw_0 <- (sum(siptw)^-1) * sum((1-A)*S*siptw*Y)
  all_siptw = all_siptw_1 - all_siptw_0
  
  all_w_1 <- (sum(w)^-1) * sum(A*S*w*Y)
  all_w_0 <- (sum(w)^-1) * sum((1-A)*S*w*Y)
  all_w = all_w_1 - all_w_0
  
  all_sw_1 <- (sum(sw)^-1) * sum(A*S*sw*Y)
  all_sw_0 <- (sum(sw)^-1) * sum((1-A)*S*sw*Y)
  all_sw = all_sw_1 - all_sw_0
  
  
  
  #### j) RESULTS ####
  
  results <- c(
    S1crude_prior_1,
    S1crude_prior_0,
    S1crude_prior,
    S0crude_prior_1,
    S0crude_prior_0,
    S0crude_prior,
    all_crude_prior_1,
    all_crude_prior_0,
    all_crude_prior,
    S1crude_1,
    S1crude_0,
    S1crude,
    S0crude_1,
    S0crude_0,
    S0crude,
    all_crude_1,
    all_crude_0,
    all_crude,
    S1_OM_1,
    S1_OM_0,
    S1_OM,
    all_OM_1,
    all_OM_0,
    all_OM,
    S1_IOW1_1,
    S1_IOW1_0,
    S1_IOW1,
    S1_IOW2_1,
    S1_IOW2_0,
    S1_IOW2,
    all_IOW1_1,
    all_IOW1_0,
    all_IOW1,
    all_IOW2_1,
    all_IOW2_0,
    all_IOW2,
    S1_DR1_1,
    S1_DR1_0,
    S1_DR1,
    S1_DR2_1,
    S1_DR2_0,
    S1_DR2,
    S1_DR3_1,
    S1_DR3_0,
    S1_DR3,
    S1_DR4_1,
    S1_DR4_0,
    S1_DR4,
    all_DR1_1,
    all_DR1_0,
    all_DR1,
    all_DR2_1,
    all_DR2_0,
    all_DR2,
    all_DR3_1,
    all_DR3_0,
    all_DR3,
    all_DR4_1,
    all_DR4_0,
    all_DR4,
    S1_GEE_crude_1,
    S1_GEE_crude_0,
    S1_GEE_crude,
    S0_GEE_crude_1,
    S0_GEE_crude_0,
    S0_GEE_crude,
    all_GEE_crude_1,
    all_GEE_crude_0,
    all_GEE_crude,
    S1_GEE_IOW1_1,
    S1_GEE_IOW1_0,
    S1_GEE_IOW1,
    S1_GEE_IOW2_1,
    S1_GEE_IOW2_0,
    S1_GEE_IOW2,
    S0_GEE_IOW1_1,
    S0_GEE_IOW1_0,
    S0_GEE_IOW1,
    S0_GEE_IOW2_1,
    S0_GEE_IOW2_0,
    S0_GEE_IOW2,
    all_GEE_IOW1_1,
    all_GEE_IOW1_0,
    all_GEE_IOW1,
    all_GEE_IOW2_1,
    all_GEE_IOW2_0,
    all_GEE_IOW2,
    S1_iosw_1,
    S1_iosw_0,
    S1_iosw,
    S1_siosw_1,
    S1_siosw_0,
    S1_siosw,
    S1_iptw_1,
    S1_iptw_0,
    S1_iptw,
    S1_siptw_1,
    S1_siptw_0,
    S1_siptw,
    S1_w_1,
    S1_w_0,
    S1_w,
    S1_sw_1,
    S1_sw_0,
    S1_sw,
    S0_iosw_1,
    S0_iosw_0,
    S0_iosw,
    S0_siosw_1,
    S0_siosw_0,
    S0_siosw,
    S0_iptw_1,
    S0_iptw_0,
    S0_iptw,
    S0_siptw_1,
    S0_siptw_0,
    S0_siptw,
    S0_w_1,
    S0_w_0,
    S0_w,
    S0_sw_1,
    S0_sw_0,
    S0_sw,
    all_iosw_1,
    all_iosw_0,
    all_iosw,
    all_siosw_1,
    all_siosw_0,
    all_siosw,
    all_iptw_1,
    all_iptw_0,
    all_iptw,
    all_siptw_1,
    all_siptw_0,
    all_siptw,
    all_w_1,
    all_w_0,
    all_w,
    all_sw_1,
    all_sw_0,
    all_sw
  )
  
  OM_values <<- c(OM_values, OM)
  DR1_values <<- c(DR1_values, DR1)
  DR2_values <<- c(DR2_values, DR2)
  GEE_sw_values <<- c(GEE_sw_values, GEE_IOW2)
  
  bootstrap_samples[[length(OM_values)]] <<- d
  
  return(results)
  
}



#### 3. RUN THE BOOTSTRAP FUNCTION ####

set.seed(1)
R = 5000

system.time({
  bs_results_boot <- boot(
    data = data,
    statistic = bs,
    R = R
    #parallel = "snow"
  )
})

saveRDS(bs_results_boot, "bs_results_boot.rds")

# in main sample, oM model not fitted in 35 iterations
# let us check what the near-0 models lead towards
# and set "extreme" values for those iterations for the percentile CIs
# OM: indices 7, 8, 9
# DR: indices 16, 17, 18, 19, 20, 21
# GEE sw: 34, 35, 36
max_index <- which.max(OM_values)
min_index <- which.min(OM_values)

max_data <- bootstrap_samples[[max_index]]
min_data <- bootstrap_samples[[min_index]]

to_exclude <- data %>%
  filter(
    hx_any_loss != 0 |
      gravida_cat != 2 |
      live_births_cat != 1 |
      any_vte == 1 |
      #GA_at_start_cat == "10 to <16 weeks" | # UNCOMMENT FOR SENSITIVITY ANALYSIS
      GA_at_start_cat == "16 to <20 weeks" |
      GA_at_start_cat == ">=20 weeks"
  )

restr_data <- data %>% 
  filter(!new_id %in% to_exclude$new_id)

summary(restr_data[S==1]$smoking_bin)
summary(restr_data[S==1]$any_preterm)

summary(max_data[S==1]$smoking_bin)
summary(max_data[S==1]$any_preterm)

summary(min_data[S==1]$smoking_bin)
summary(min_data[S==1]$any_preterm)

test_val <- data.frame(test_value = OM_values)

p <- ggplot(test_val, aes(x = test_value)) +
  geom_histogram(
    bins = 30,            
    fill = "#4C72B0",         
    color = "white",          
    alpha = 0.9
  ) +
  labs(
    title = "Distribution of Bootstrapped OM Estimate",
    x = "Estimate",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
p

ggsave("bootsrapped_OM_dist.png", plot = p, width = 8, height = 6, dpi = 300)

# assign 99 to half and -99 to other
# to failed iterations (so they are not skipped when building CIs)
# OM: indices 19, 20, 21
# DR1 & DR2: indices 37-42 (won't work if OM did not work)
# DR3 & DR4 (w and sw): 43-48
# DR all: 49-60
# GEE sw: 73, 74, 75

# note: failures to represent failed MODELS, but can be multiple in same iteration

test <- bs_results_boot$t
test2 <- as.data.frame(test) %>% filter(is.na(V19))

cols_to_replace <- c(19, 20, 21, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 73, 74, 75)
test3 <- bs_results_boot$t[, cols_to_replace]
test4 <- as.data.frame(test) %>%
  filter(if_any(everything(), is.na))

for (col in cols_to_replace) {
  
  na_idx <- which(is.na(bs_results_boot$t[, col]))
  if (length(na_idx) > 1) {
    
    bs_results_boot$t[sample(na_idx, length(na_idx) %/% 2), col] <- 99
    bs_results_boot$t[is.na(bs_results_boot$t[, col]), col] <- -99
    
  } else if (length(na_idx) == 1) {
    
    bs_results_boot$t[na_idx, col] <- sample(c(99, -99), 1)
    
  }
}

#### 4. GET CONFIDENCE INTERVALS ####

bs_col_names <- c(
  "all_crude_prior_1",
  "all_crude_prior_0",
  "all_crude_prior",
  "S1crude_prior_1",
  "S1crude_prior_0",
  "S1crude_prior",
  "S0crude_prior_1",
  "S0crude_prior_0",
  "S0crude_prior",
  "all_crude_1",
  "all_crude_0",
  "all_crude",
  "S1_crude_1",
  "S1_crude_0",
  "S1_crude",
  "S0_crude_1",
  "S0_crude_0",
  "S0_crude",
  "OM_1",
  "OM_0",
  "OM",
  "all_OM_1",
  "all_OM_0",
  "all_OM",
  "IOW1_1",
  "IOW1_0",
  "IOW1",
  "IOW2_1",
  "IOW2_0",
  "IOW2",
  "all_IOW1_1",
  "all_IOW1_0",
  "all_IOW1",
  "all_IOW2_1",
  "all_IOW2_0",
  "all_IOW2",
  "DR1_1",
  "DR1_0",
  "DR1",
  "DR2_1",
  "DR2_0",
  "DR2",
  "DR3_1",
  "DR3_0",
  "DR3",
  "DR4_1",
  "DR4_0",
  "DR4",
  "all_DR1_1",
  "all_DR1_0",
  "all_DR1",
  "all_DR2_1",
  "all_DR2_0",
  "all_DR2",
  "all_DR3_1",
  "all_DR3_0",
  "all_DR3",
  "all_DR4_1",
  "all_DR4_0",
  "all_DR4",
  "GEE_crude_1",
  "GEE_crude_0",
  "GEE_crude",
  "S1GEE_crude_1",
  "S1GEE_crude_0",
  "S1GEE_crude",
  "S0GEE_crude_1",
  "S0GEE_crude_0",
  "S0GEE_crude",
  "GEE_IOW1_1",
  "GEE_IOW1_0",
  "GEE_IOW1",
  "GEE_IOW2_1",
  "GEE_IOW2_0",
  "GEE_IOW2",
  "all_GEE_IOW1_1",
  "all_GEE_IOW1_0",
  "all_GEE_IOW1",
  "all_GEE_IOW2_1",
  "all_GEE_IOW2_0",
  "all_GEE_IOW2",
  "S1_iosw_1",
  "S1_iosw_0",
  "S1_iosw",
  "S1_siosw_1",
  "S1_siosw_0",
  "S1_siosw",
  "S1_iptw_1",
  "S1_iptw_0",
  "S1_iptw",
  "S1_siptw_1",
  "S1_siptw_0",
  "S1_siptw",
  "S1_w_1",
  "S1_w_0",
  "S1_w",
  "S1_sw_1",
  "S1_sw_0",
  "S1_sw",
  "S0_iosw_1",
  "S0_iosw_0",
  "S0_iosw",
  "S0_siosw_1",
  "S0_siosw_0",
  "S0_siosw",
  "S0_iptw_1",
  "S0_iptw_0",
  "S0_iptw",
  "S0_siptw_1",
  "S0_siptw_0",
  "S0_siptw",
  "S0_w_1",
  "S0_w_0",
  "S0_w",
  "S0_sw_1",
  "S0_sw_0",
  "S0_sw",
  "iosw_1",
  "iosw_0",
  "iosw",
  "siosw_1",
  "siosw_0",
  "siosw",
  "iptw_1",
  "iptw_0",
  "iptw",
  "siptw_1",
  "siptw_0",
  "siptw",
  "w_1",
  "w_0",
  "w",
  "sw_1",
  "sw_0",
  "sw",
  "test_iosw_1",
  "test_iosw_0",
  "test_iosw",
  "test_siosw_1",
  "test_siosw_0",
  "test_siosw",
  "test_iptw_1",
  "test_iptw_0",
  "test_iptw",
  "test_siptw_1",
  "test_siptw_0",
  "test_siptw",
  "test_w_1",
  "test_w_0",
  "test_w",
  "test_sw_1",
  "test_sw_0",
  "test_sw"
)

bs_ci_boot <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(bs_ci_boot) <- c('estimate', 'lower_ci', 'upper_ci')

for (i in 1:length(bs_results_boot$t0)) {

  stat_ci <- data.frame(estimate = NA,
                        lower_ci = NA,
                        upper_ci = NA)
  
  conf.int <- boot.ci(bs_results_boot, type = 'perc', index = i)
  CI <- conf.int$percent
  
  stat_ci$estimate <- bs_results_boot$t0[i]
  stat_ci$lower_ci <- t(CI[, 4])[[1]]
  stat_ci$upper_ci <- t(CI[, 5])[[1]]
  bs_ci_boot <- rbind(bs_ci_boot, stat_ci)
  
}

rownames(bs_ci_boot) <- bs_col_names

bs_ci_boot %<>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

bs_ci_boot %<>%
  mutate(across(where(is.numeric), ~ .x*100))

write.table (bs_ci_boot , "bs_ci.csv", col.names = T, row.names=F, append= F, sep=',')


