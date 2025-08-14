## ---------------------------
##
## Program: AFFIRM Re-Analysis
## Author: Gwen Aubrac
##
## Date Created: 2025-06-23
##
## ---------------------------
##
## Notes: 
## Weighted GEE tutorial: 
## https://rstudio-pubs-static.s3.amazonaws.com/415697_d6fed7f7b9754c0999cc2b4075013c81.html
##
## Augmented odds weights adapted from:
## Dahabreh IJ, Robertson SE, Stuart EA, Hernan MA. 
## Transporting inferences from a randomized trial to a new target population. 
## Published online May 1, 2018. doi:10.48550/arXiv.1805.00550
##
## For full code to get estimates in manuscript, please see bootstrap code.
## This is a more commented version that also includes the data preparation process,
## but not all of the different estimators.
##
## ---------------------------

#### 1. LOAD PACKAGES ####

options(scipen = 999)

library(PSweight)
library(dplyr)
library(magrittr)
library(table1)
library(writexl)
library(tidysmd)
library(cobalt)
library(ggplot2)
library(glmtoolbox)

#### 2. READ AND EXPLORE DATA ####

# set working directory to analysis, among:
# main_results
# target_multi_site
# aspirin
# no_aspirin
# restrict_GA_10w

setwd("Z:/RPLATT/GWEN-OHRI/main_results")
data <- read.csv("../AFFIRM_GWEN.csv")

data %<>%
  mutate(new_id = row_number())

table(data$unique_id)

# we are not analyzing patients from ETHIGII, so we remove them
data %<>%
  filter(unique_id != "ETHIGII")

# explore the data a little
data %>%
  group_by(unique_id) %>%
  summarize(
    min_value = min(ga_at_start__weeks_, na.rm = TRUE),
    mean_value = mean(ga_at_start__weeks_, na.rm = TRUE),
    max_value = max(ga_at_start__weeks_, na.rm = TRUE)
  )

#### 3. EXCLUDE PATIENTS WITH COMPETING EVENT OR MISSING OUTCOME DATA ####

# the outcome is a composite of:
# early-onset PE (<34 weeks) or severe PE
# birth of SGA neonate (<5th percentile)
# placental abruption
# late pregnancy loss (>= 20 weeks)

# patients by trial prior to exclusions
data %>% 
  group_by(unique_id) %>% 
  summarize(
    N = n()
  )

# first we have to create some variables to help us define the competing event

# determine gestational age at PE onset using trial start and onset dates
# note: will be NA for patients who did not have PE
data$active_trial_start <- as.Date(data$active_trial_start, format = "%d%b%Y")
data$pre_eclampsia_onset <- as.Date(data$pre_eclampsia_onset, format = "%d%b%Y")

data$PE_days_elapsed <- as.numeric(data$pre_eclampsia_onset - data$active_trial_start)
data$PE_weeks_elapsed <- round(data$PE_days_elapsed/7, 2)
data$GA_at_PE <- data$ga_at_start__weeks_ + data$PE_weeks_elapsed

# determine GA at pregnancy loss using trial start and "delivery"/pregnancy outcome date
# for now, considering any type of pregnancy loss
data$pregnancy_outcome_date <- as.Date(data$pregnancy_outcome_date, format = "%d%b%Y")
data$del_days_elapsed <- as.numeric(data$pregnancy_outcome_date - data$active_trial_start)
data$del_weeks_elapsed <- data$del_days_elapsed / 7
data$GA_at_delivery <- data$ga_at_start__weeks_ + data$del_weeks_elapsed

table(data$pregnancy_outcome)

data %<>%
  mutate(
    pregnancy_loss = if_else( # 1 = delivery
      pregnancy_outcome == 2 | # 2 = unexplained loss
        pregnancy_outcome == 3 | # 3 = explained loss
        pregnancy_outcome == 4 | # 4 = social termination
        pregnancy_outcome == 5, # 5 = medical termination
      1, 
      0
    )
  )

table(data$pregnancy_loss)

# let us exclude women with a competing event 
# namely those with any loss prior to 20 weeks, termination unrelated to outcome,
# or missing data on outcome timing

data %>% 
  group_by(unique_id) %>% 
  summarize(
    missing_GA_start = sum(is.na(ga_at_start__weeks_)),
    missing_outcome_date = sum(is.na(pregnancy_outcome_date))
  )

data %>% 
  group_by(unique_id) %>% 
  summarize(
    explained = sum(pregnancy_outcome == 2 & GA_at_delivery <20 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))),
    unexplained = sum(pregnancy_outcome == 3 & GA_at_delivery < 20 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))),
    social = sum(pregnancy_outcome == 4 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))),
    medical = sum(pregnancy_outcome == 5 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))),
    missing = sum(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date)),
    any = sum(explained + unexplained + social + medical + missing)
  )

# let us exclude women who had explained or unexplained early loss,
# medical or social termination, or missing data on outcome date or GA at start
test <- data %>% 
  filter(
    (pregnancy_outcome == 2 & GA_at_delivery <20 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))) |
      (pregnancy_outcome == 3 & GA_at_delivery < 20 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))) |
      (pregnancy_outcome == 4 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))) |
      (pregnancy_outcome == 5 & !(is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))) |
      (is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))
  )

# apply those exclusions
data %<>%
  filter(!(new_id %in% test$new_id))

# patients by trial after exclusions
data %>% 
  group_by(unique_id) %>% 
  summarize(
    N = n()
  )

#### 4. DEFINE OUTCOME ####

# the outcome is a composite of:
# early-onset PE (<34 weeks) or severe PE
# birth of SGA neonate (<5th percentile)
# placental abruption
# late pregnancy loss (>= 20 weeks)

# see missing outcome data by trial
data %>% 
  group_by(unique_id) %>% 
  summarize(
    missing_PE = sum(is.na(pre_eclampsia) | pre_eclampsia == 99),
    missing_GA_PE = sum(is.na(GA_at_PE)),
    missing_severe_PE = sum(is.na(severe_pre_eclampsia) | severe_pre_eclampsia == 99),
    missing_sga5 = sum(is.na(sga__5) | sga__5 == 99),
    missing_PA = sum(is.na(abruption_no_delivery) | abruption_no_delivery == 99),
    missing_PA_del = sum(is.na(abruption_with_delivery) | abruption_with_delivery == 99),
    missing_pregnancy = sum(is.na(pregnancy_outcome) | pregnancy_outcome == 99)
  )

# assume that missing means did not occur
data$sga__5[data$sga__5 == 99 | is.na(data$sga__5)] <- 1
data$abruption_no_delivery[data$abruption_no_delivery == 99 | is.na(data$abruption_no_delivery)] <- 1
data$abruption_with_delivery[data$abruption_with_delivery == 99 | is.na(data$abruption_with_delivery)] <- 1

# define early PE as PE occurring prior to 34 weeks GA
data %<>%
  mutate(early_PE = if_else(pre_eclampsia == 2 & GA_at_PE < 34, 1, 0))

table(data$early_PE)

data %>% 
  group_by(unique_id) %>% 
  summarize(
    events = sum(early_PE == 1)
  )

table(data$severe_pre_eclampsia)
table(data$sga__5)
table(data$abruption_no_delivery)
table(data$abruption_with_delivery)

data %<>%
  mutate(late_loss = if_else((pregnancy_outcome == 2 | pregnancy_outcome == 3) & GA_at_delivery >=20, 1, 0))

table(data$late_loss)

# finally define composite outcome as any one of the above
data %<>%
  mutate(
    outcome = if_else(
      early_PE == 1 |
        severe_pre_eclampsia == 2 |
        sga__5 == 2 |
        abruption_no_delivery == 2 |
        abruption_with_delivery == 2 |
        late_loss == 1, 
      1, 
      0
    )
  )

table(data$outcome)



#### 5. DEFINE TREATMENT ####

# let us define treatment (LMWH or not)
# and also create a variable for single-site trials (NOH-AP and NOH-PE trials)

data %<>%
  mutate(
    single = if_else(unique_id == "NOH-AP" | unique_id == "NOH-PE", 1, 0),
    trt = if_else(allocation == 1 | allocation == 2, 1, 0), # 1= LMWH, 2 = LMWH + Aspirin
  )

table(data$single)
table(data$trt)

# create variable for concomitant aspirin use (covariate and EMM/subgroup analysis)
data %<>%
  mutate(
    aspirin = if_else(allocation == 2 | allocation == 3, 1, 0)
  )

table(data$aspirin)



#### 6. DESCRIBE STUDY POPULATION ####

data %<>%
  mutate(
    bmi = round((baseline_weight_kg_ / height__cm_ / height__cm_) * 10000, 2)
  )


tab1_data <- data

tab1_data %<>%
  mutate(across(c(
    hx_late_loss_12_weeks, 
    hx_late_loss_16_weeks, 
    hx_sga_10th,
    hx_pre_eclampsia,
    hx_abruption,
    hx_abruption_delivery,
    fvl, 
    pgm, 
    anti_thrombin,
    protein_c,
    protein_s,
    apla,
    other_thrombo_philia,
    smoking, 
    chronic_hyper_tension,
    diabetes,
    vte,
    family_history_of_vte,
    family_history_of_arterial_disea,
    hx_losses,
    hx_late_loss_20,
    hx_sga_5,
    hx_sga_3,
    hx_preterm_34,
    hx_preterm_37,
    hx_severe_pre_eclampsia,
    hx_early_pre_eclampsia,
    hx_abruption,
    aspirin
    ), as.factor))

# describe study population before formatting our covariates
tab1_covs <- c(
  "ga_at_start__weeks_",
  "age__years_",
  "race",
  "bmi",
  "hx_late_loss_12_weeks", 
  "hx_late_loss_16_weeks",
  "hx_sga_10th",
  "hx_pre_eclampsia",
  "hx_abruption_delivery",
  "fvl",
  "pgm",
  "anti_thrombin",
  "protein_c",
  "protein_s",
  "apla",
  #"other_thrombo_philia",
  "smoking",
  "chronic_hyper_tension",
  "diabetes",
  "systolic_baseline",
  "diastolic_baseline", 
  "vte", 
  "family_history_of_vte",
  "family_history_of_arterial_disea",
  "gravida",
  "hx_losses",
  "hx_late_loss_20",
  "live_births",
  "hx_sga_5",
  "hx_sga_3",
  "hx_preterm_34",
  "hx_preterm_37",
  "hx_severe_pre_eclampsia",
  "hx_early_pre_eclampsia",
  "hx_abruption",
  "aspirin"
)

caption <- 'Baseline Characteristics of Cohort by Single- vs Multi-Site Trials'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|single"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1

#write.table (tab1 , "./results/table1_by_site_prior.csv", col.names = T, row.names=F, append= F, sep=',')

# get smds
factor_covs <- c(
  "race",
  "hx_late_loss_12_weeks", 
  "hx_late_loss_16_weeks",
  "hx_sga_10th",
  "hx_pre_eclampsia",
  "hx_abruption_delivery",
  "fvl",
  "pgm",
  "anti_thrombin",
  "protein_c",
  "protein_s",
  "apla",
  "smoking",
  "chronic_hyper_tension",
  "diabetes",
  "vte", 
  "family_history_of_vte",
  "family_history_of_arterial_disea",
  "gravida",
  "hx_losses",
  "hx_late_loss_20",
  "live_births",
  "hx_sga_5",
  "hx_sga_3",
  "hx_preterm_34",
  "hx_preterm_37",
  "hx_severe_pre_eclampsia",
  "hx_early_pre_eclampsia",
  "hx_abruption",
  "aspirin"
)

smd_data <- data %>% 
  mutate(across(all_of(factor_covs), as.factor))

print(tidy_smd(smd_data, tab1_covs, .group = single, na.rm = TRUE), n = 40)

caption <- 'Baseline Characteristics of Cohort by Trial'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|unique_id"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1

#write.table (tab1 , "table1_by_trial_prior.csv", col.names = T, row.names=F, append= F, sep=',')


caption <- 'Baseline Characteristics of Cohort by Treatment'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|trt"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1

#write.table (tab1 , "table1_by_trt_prior.csv", col.names = T, row.names=F, append= F, sep=',')



#### 7. FORMAT COVARIATES FOR ANALYSES ####

# missing covariate data by trial
data %>% 
  group_by(unique_id) %>% 
  summarize(
    N = n()
  )

data %>% 
  group_by(unique_id) %>% 
  summarize(
    misS_BMI = sum(is.na(bmi)),
    miss_race = sum(is.na(race) | race == ""),
    miss_paternity = sum(is.na(paternity) | paternity == 99),
    miss_loss_12 = sum(is.na(hx_late_loss_12_weeks) | hx_late_loss_12_weeks == 99),
    miss_loss_16 = sum(is.na(hx_late_loss_16_weeks) | hx_late_loss_16_weeks == 99),
    misS_sga10 = sum(is.na(hx_sga_10th) | hx_sga_10th == 99),
    miss_sga3 = sum(is.na(hx_sga_3) | hx_sga_3 == 99),
    miss_sga5 = sum(is.na(hx_sga_5) | hx_sga_5 == 99),
    miss_smoking = sum(is.na(smoking) | smoking == 99),
    miss_hypertension = sum(is.na(chronic_hyper_tension) | chronic_hyper_tension == 99),
    misS_vte = sum(is.na(vte) | vte == 99),
    miss_hx_vte = sum(is.na(family_history_of_vte) | family_history_of_vte == 99),
    miss_hx_art = sum(is.na(family_history_of_arterial_disea) | family_history_of_arterial_disea == 99),
    miss_hx_loss20 = sum(is.na(hx_late_loss_20) | hx_late_loss_20 == 99),
    miss_pt_34 = sum(is.na(hx_preterm_34) | hx_preterm_34 == 99),
    miss_pt_37 = sum(is.na(hx_preterm_37) | hx_preterm_37 == 99),
    miss_PE = sum(is.na(hx_pre_eclampsia) | hx_pre_eclampsia == 99),
    miss_severe_PE = sum(is.na(hx_severe_pre_eclampsia) | hx_severe_pre_eclampsia == 99),
    miss_PA = sum(is.na(hx_abruption) | hx_abruption == 99),
    miss_PA_del = sum(is.na(hx_abruption_delivery) | hx_abruption_delivery == 99),
  ) %>% 
  print(width = Inf)

# remove 1 patient with missing age
data %<>% filter(!is.na(age__years_))

data %>% 
  summarize(
    misS_BMI = sum(is.na(bmi)),
    miss_race = sum(is.na(race) | race == ""),
    miss_paternity = sum(is.na(paternity) | paternity == 99),
    miss_loss_12 = sum(is.na(hx_late_loss_12_weeks) | hx_late_loss_12_weeks == 99),
    miss_loss_16 = sum(is.na(hx_late_loss_16_weeks) | hx_late_loss_16_weeks == 99),
    misS_sga10 = sum(is.na(hx_sga_10th) | hx_sga_10th == 99),
    misS_sga5 = sum(is.na(hx_sga_5) | hx_sga_5 == 99),
    miss_sga3 = sum(is.na(hx_sga_3) | hx_sga_3 == 99),
    miss_smoking = sum(is.na(smoking) | smoking == 99),
    miss_hypertension = sum(is.na(chronic_hyper_tension) | chronic_hyper_tension == 99),
    misS_vte = sum(is.na(vte) | vte == 99),
    miss_hx_vte = sum(is.na(family_history_of_vte) | family_history_of_vte == 99),
    miss_hx_art = sum(is.na(family_history_of_arterial_disea) | family_history_of_arterial_disea == 99),
    miss_hx_loss20 = sum(is.na(hx_late_loss_20) | hx_late_loss_20 == 99),
    miss_pt_34 = sum(is.na(hx_preterm_34) | hx_preterm_34 == 99),
    miss_pt_37 = sum(is.na(hx_preterm_37) | hx_preterm_37 == 99),
    miss_PE = sum(is.na(hx_pre_eclampsia) | hx_pre_eclampsia == 99),
    miss_severe_PE = sum(is.na(hx_severe_pre_eclampsia) | hx_severe_pre_eclampsia == 99),
    miss_PA = sum(is.na(hx_abruption) | hx_abruption == 99),
    miss_PA_del = sum(is.na(hx_abruption_delivery) | hx_abruption_delivery == 99),
  )

data %>% 
  group_by(unique_id) %>% 
  summarize(
    miss_fvl = sum(is.na(fvl) | fvl == 4 | fvl == 99),
    miss_pgm = sum(is.na(pgm) | pgm == 4 | pgm == 99),
    miss_antithrombin = sum(is.na(anti_thrombin) | anti_thrombin == 3 | anti_thrombin == 99),
    miss_protein_c = sum(is.na(protein_c) | protein_c == 3 | protein_c == 99),
    miss_protein_s = sum(is.na(protein_s) | protein_s == 3 | protein_s == 99),
    miss_apla = sum(is.na(apla) | apla == 3 | apla == 99),
    miss_other = sum(is.na(other_thrombo_philia) | other_thrombo_philia == 3 | other_thrombo_philia == 99)
  )

data %>% 
  summarize(
    miss_fvl = sum(is.na(fvl) | fvl == 4 | fvl == 99),
    miss_pgm = sum(is.na(pgm) | pgm == 4 | pgm == 99),
    miss_antithrombin = sum(is.na(anti_thrombin) | anti_thrombin == 3 | anti_thrombin == 99),
    miss_protein_c = sum(is.na(protein_c) | protein_c == 3 | protein_c == 99),
    miss_protein_s = sum(is.na(protein_s) | protein_s == 3 | protein_s == 99),
    miss_apla = sum(is.na(apla) | apla == 3 | apla == 99),
    miss_other = sum(is.na(other_thrombo_philia) | other_thrombo_philia == 3 | other_thrombo_philia == 99)
  )

# we will not use covariates that were systematically not collected by some trials
# namely race, paternity, hypertension, arterial disease, severe PE, other thrombophilia, family history of VTE
# also remove diabetes (since no one has it)

# for the others, assume that if missing then did not happen
# please refer to AFFIRM data dictionary for explanation of ref category
data$hx_late_loss_12_weeks[is.na(data$hx_late_loss_12_weeks) == T] <- 0
data$hx_late_loss_16_weeks[is.na(data$hx_late_loss_16_weeks) == T] <- 0
data$hx_sga_10th[is.na(data$hx_sga_10th) == T] <- 1
data$hx_sga_10th[data$hx_sga_10th == 99] <- 1
data$smoking[is.na(data$smoking) == T] <- 1
data$vte[is.na(data$vte) == T] <- 1
data$hx_late_loss_20[is.na(data$hx_late_loss_20) == T] <- 0
data$hx_sga_5[is.na(data$hx_sga_5) == T] <- 1
data$hx_sga_5[data$hx_sga_5 == 99] <- 1
data$hx_sga_3[is.na(data$hx_sga_3) == T] <- 1
data$hx_sga_3[data$hx_sga_3 == 99] <- 1
data$hx_preterm_34[is.na(data$hx_preterm_34) == T] <- 0
data$hx_preterm_37[is.na(data$hx_preterm_37) == T] <- 0
data$hx_pre_eclampsia[is.na(data$hx_pre_eclampsia) == T] <- 1
data$hx_early_pre_eclampsia[is.na(data$hx_early_pre_eclampsia) == T] <- 1
data$hx_abruption[is.na(data$hx_abruption) == T] <- 1

# for testing, mark missing to be same as "not tested"
# and not tested as negative
# here we are making big assumption that if no data, did not occur
data$fvl[is.na(data$fvl) == T] <- 4
data$fvl[data$fvl == 4] <- 1

data$pgm[is.na(data$pgm) == T] <- 4
data$pgm[data$pgm == 4] <- 1

data$anti_thrombin[is.na(data$anti_thrombin) == T] <- 3
data$anti_thrombin[data$anti_thrombin == 3] <- 1

data$protein_c[is.na(data$protein_c) == T] <- 3
data$protein_c[data$protein_c == 3] <- 1

data$protein_s[is.na(data$protein_s) == T] <- 3
data$protein_s[data$protein_s == 3] <- 1

data$apla[is.na(data$apla) == T] <- 3
data$apla[data$apla == 3] <- 1

data$other_thrombo_philia[is.na(data$other_thrombo_philia) == T] <- 3
data$other_thrombo_philia[data$other_thrombo_philia == 3] <- 1

# for thrombophilia, count how many abnormal results for each patient
data %<>%
  mutate(
    fvl_flag = if_else(fvl == 2 | fvl == 3, 1, 0),
    pgm_flag = if_else(pgm == 2 | pgm == 3, 1, 0),
    anti_thrombin_flag = if_else(anti_thrombin == 2, 1, 0),
    protein_c_flag = if_else(protein_c == 2, 1, 0),
    protein_s_flag = if_else(protein_s == 2, 1, 0),
    apla_flag = if_else(apla == 2, 1, 0)
    )

data %<>%
  group_by(new_id) %>% 
  mutate(
    thrombo_count = sum(
      fvl_flag, 
      pgm_flag,
      anti_thrombin_flag,
      protein_c_flag,
      protein_s_flag,
      apla_flag
    )
  ) %>% 
  ungroup()

table(data$thrombo_count)

# covariate cleaning / categorization

data <- data %>%
  mutate(
    any_preterm = if_else(hx_preterm_34 >= 1 | hx_preterm_37 >=1, 1, 0),
    age_cat = case_when(
      age__years_ < 25 ~ "<25",
      age__years_ >= 25 & age__years_ < 35 ~ "25 to <35",
      age__years_ >= 35 ~ ">=35"
    ),
    above_30 = case_when(
      age__years_ < 30 ~ 0,
      age__years_ >= 30 ~ 1,
    ),
    any_pe = if_else(hx_pre_eclampsia == 2 | hx_early_pre_eclampsia == 2, 1, 0),
    any_sga = if_else(hx_sga_10th == 2 | hx_sga_5 == 2 | hx_sga_3 == 2, 1, 0),
    bmi_cat = case_when(
      is.na(bmi) ~ NA,
      bmi < 18.5 ~ "<18.5",
      bmi >= 18.5 & bmi < 24.9 ~ "18.5-24.9",
      bmi >= 25 & bmi < 29.9 ~ "25-29.9",
      bmi >= 30 & bmi < 34.9 ~ "30-34.9",
      bmi >= 34.9 & bmi <39.9 ~ "34.9-39.9",
      bmi >=40 ~ ">=40"
    ),
    bmi_above_30 = case_when(
      is.na(bmi) | bmi < 30 ~ 0,
      bmi >= 30 ~ 1
    ),
    GA_at_start_cat = case_when(
      ga_at_start__weeks_ < 10 ~ "<10 weeks",
      ga_at_start__weeks_ >= 10 & ga_at_start__weeks_ < 16 ~ "10 to <16 weeks",
      ga_at_start__weeks_ >= 16 & ga_at_start__weeks_< 20 ~ "16 to <20 weeks",
      ga_at_start__weeks_ >= 20 ~ ">=20 weeks"
    ),
    hx_loss_12_cat = case_when(
      hx_late_loss_12_weeks == 0 ~ "0",
      hx_late_loss_12_weeks == 1 ~ "1",
      hx_late_loss_12_weeks >= 2 ~ "2+"
    ),
    hx_loss_16_cat = case_when(
      hx_late_loss_16_weeks == 0 ~ "0",
      hx_late_loss_16_weeks == 1 ~ "1",
      hx_late_loss_16_weeks >= 2 ~ "2+"
    ),
    gravida_cat = case_when(
      gravida == 0 ~ "0",
      gravida == 1 ~ "1",
      gravida == 2 ~ "2",
      gravida == 3 ~ "3",
      gravida >= 4 ~ "4+"
    ),
    hx_losses_cat = case_when(
      hx_losses == 0 ~ "0",
      hx_losses == 1 ~ "1",
      hx_losses >= 2 ~ "2+"
    ),
    live_births_cat = case_when(
      live_births == 0 ~ "0",
      live_births == 1 ~ "1",
      live_births >= 2 ~ "2+"
    ),
    any_preterm_34 = case_when(
      hx_preterm_34 == 0 ~ 0,
      hx_preterm_34 >= 1 ~ 1
    ),
    any_preterm_37 = case_when(
      hx_preterm_37 == 0 ~ 0,
      hx_preterm_37 >= 1 ~ 1
    ),
    hx_loss_20_cat = case_when(
      hx_late_loss_20 == 0 ~ "0",
      hx_late_loss_20 == 1 ~ "1",
      hx_late_loss_20 >= 2 ~ "2+"
    ),
    hx_any_loss = if_else(hx_late_loss_20 >= 1 | hx_late_loss_12_weeks >= 1 | hx_late_loss_16_weeks >=1 | hx_losses >= 1, 1, 0), 
    thrombophilia = case_when(
      fvl == 1 & pgm == 1 & anti_thrombin == 1 & protein_c == 1 & protein_s == 1 & apla == 1 & other_thrombo_philia == 1 ~ "none",
      fvl == 2 | pgm == 2 ~ "weak",
      protein_c == 2 | protein_s == 2 | other_thrombo_philia == 2 ~ "moderate",
      anti_thrombin == 2 | apla == 2 | fvl == 3 | pgm == 3 | thrombo_count > 1 ~ "strong"
    ),
    thrombophilia_bin = case_when(
      thrombophilia == "none" ~ "no",
      thrombophilia == "weak" | thrombophilia == "moderate" | thrombophilia == "strong" ~ "yes"
    ),
    smoking_bin = case_when(
      smoking == 1 | smoking == 3 ~ "no",
      smoking == 2 ~ "yes"
    ),
    hx_any_abruption = if_else(hx_abruption == 2 | hx_abruption_delivery == 2, 1, 0),
    any_vte = if_else(vte != 1, 1, 0)
  )

data %<>%
  mutate(across(c(
    any_preterm,
    age_cat,
    above_30,
    any_pe,
    any_sga,
    bmi_cat,
    bmi_above_30,
    GA_at_start_cat,
    hx_loss_12_cat,
    hx_loss_16_cat,
    gravida_cat,
    hx_losses_cat,
    live_births_cat,
    any_preterm_34,
    any_preterm_37,
    hx_loss_20_cat,
    hx_any_loss,
    thrombophilia,
    thrombophilia_bin,
    race,
    hx_sga_10th,
    hx_pre_eclampsia,
    hx_abruption_delivery,
    thrombophilia,
    smoking,
    smoking_bin,
    vte,
    gravida,
    hx_sga_5,
    hx_sga_3,
    hx_early_pre_eclampsia,
    hx_abruption,
    aspirin,
    hx_any_abruption,
    any_vte
  ), as.factor))

# get SMDs
# tidy_smd(data, tab1_covs, .group = single, na.rm = T)
# rm(caption, tab, tab_formula, tab_data)

# save cleaned data (for bootstrapping later)
# data_clean <- data
# saveRDS(data_clean, file = "../data_clean.R")
data <- readRDS("../data_clean.R")

#### 8. RESTRICTIONS PRIOR TO ANALYSES ####

# examine whether some patients are systematically not represented in single-site trials
# as will violate positivity in the weighting model
covs_list <- c(
  "any_preterm",
  "age_cat",
  "above_30",
  "any_pe",
  "any_sga",
  "bmi_above_30",
  "GA_at_start_cat",
  "gravida_cat",
  "live_births_cat",
  "hx_any_loss",
  "thrombophilia_bin",
  "hx_any_abruption",
  "smoking_bin",
  "any_vte",
  "aspirin"
)

# covariate distribution by trial
caption <- 'Baseline Characteristics (Cleaned Covs) of Cohort by Trial'
tab2_formula <- as.formula(paste("~", paste(covs_list, collapse = "+"), "|unique_id"))
tab2 <- table1(
  tab2_formula,
  data = data,
  overall = c(right = 'Total'),
  caption = caption
)

tab2

#write.table (tab2 , "table1_by_trt_clean_covs.csv", col.names = T, row.names=F, append= F, sep=',')

# exclusion 1: any pregnancy loss
# exclusion 2: any gravida >2
# exclusion 3: GA at active start >16 weeks
# note: will assume GA <10 and 10-16 are equivalent

analytic_data <- data

table(analytic_data$hx_any_loss != 0)
table(analytic_data$gravida_cat != 2)
table(analytic_data$live_births_cat != 1)
table(analytic_data$any_vte == 1)
table(analytic_data$GA_at_start_cat == "10 to <16 weeks")
table(analytic_data$GA_at_start_cat == "16 to <20 weeks")
table(analytic_data$GA_at_start_cat == ">=20 weeks")


to_exclude <- analytic_data %>% 
  filter(
    hx_any_loss != 0 |
    gravida_cat != 2 |
    live_births_cat != 1 |
    any_vte == 1 |
    #GA_at_start_cat == "10 to <16 weeks" |
    GA_at_start_cat == "16 to <20 weeks" |
    GA_at_start_cat == ">=20 weeks"
  )

analytic_data %<>%
  filter(!new_id %in% to_exclude$new_id)

covs_list <- c(
  "age__years_",
  "age_cat",
  "above_30",
  "any_pe",
  "any_preterm",
  "any_sga",
  "hx_any_abruption",
  "bmi_above_30",
  "thrombophilia_bin",
  "smoking_bin",
  "aspirin",
  "hx_any_abruption"
)

# covariates by single vs multi-site trials after restrictions
caption <- "Baseline Characteristics by Site After Restrictions"
tab3_formula <- as.formula(paste("~", paste(covs_list, collapse = "+"), "|single"))
tab3 <- table1(
  tab3_formula,
  data = analytic_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab3

print(tidy_smd(analytic_data, covs_list, .group = single, na.rm = TRUE), n = 40)

#write.table (tab3 , "table1_by_site_restr.csv", col.names = T, row.names=F, append= F, sep=',')

# covariates by trial after applying exclusions
caption <- "Baseline Characteristcs by Trial After Restrictions"
tab3_formula <- as.formula(paste("~", paste(covs_list, collapse = "+"), "|unique_id"))
tab3 <- table1(
  tab3_formula,
  data = analytic_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab3
#write.table (tab3 , "table1_by_trial_restr.csv", col.names = T, row.names=F, append= F, sep=',')



#### 9. COMPUTE IOSW WEIGHTS AND CHECK BALANCE ####

# adapted from Dahabreh et al.
# thinking of S as "selection into multi-site trial", so:
# single site patients (single = 1 -> S = 0) -- this is our "target" population
# multi site patients (single = 0 -> S = 1) is our "study" population
# so our goal is to extend inference from multi-site to single-site patients
analytic_data %<>%
  mutate(S = if_else(single == 1, 0, 1))

# for sensitivity where reverse S, do this instead:
# analytic_data %<>%
#   mutate(S = if_else(single == 1, 1, 0))

analytic_covs_list <- c(
  "age__years_",
  "any_pe",
  "any_preterm",
  "any_sga",
  "bmi_above_30",
  "smoking_bin",
  "aspirin"
)

# we tested different selection models to achieve balance
# and settled on using splines to model age
library(splines)
selection_model <- as.formula("S ~ 
                              ns(age__years_, df = 3) + 
                              any_pe + 
                              any_preterm + 
                              any_sga + 
                              bmi_above_30 + 
                              smoking_bin + 
                              aspirin +
                              aspirin * ns(age__years_, df = 2)
                              ")

# generate weights
w_reg <- glm(selection_model, family = binomial(), data = analytic_data)
summary(w_reg)

ps <- predict(w_reg, newdata = analytic_data, type = "response")
analytic_data$ps <- ps
analytic_data %<>%
  mutate(
    iosw= if_else(
      S == 1, (1-ps)/ps,
      0 # weight of 0 for S0 patients
    )
  )

analytic_data %<>%
  mutate(
    iosw_all= if_else(
      S == 1, (1-ps)/ps,
      1 # weight of 1 for S0 patients
    )
  )

# stabilize weights: multiply by marginal ODDS of selection for S=1 patients
marg_s_model <- glm(S ~ 1, family = binomial(), data = analytic_data)
marg_ps <- predict(marg_s_model, newdata = analytic_data, type = "response")
odds_s <- marg_ps/(1-marg_ps)
odds_s

analytic_data %<>%
  mutate(
    siosw = if_else(S == 1, odds_s*iosw, 0), # weight of 0 for S0 patients
    siosw_all = if_else(S == 1, odds_s*iosw, 1) # weight of 1 for S0 patients
  )

# check N of pseudo populations, which should be:
# N = 359 for IOSW in S=0 (weight of 1 for all)
# N = 359 for SIOSW in S=0 (weight of 1 for all)
# N = 359 for IOSW in S=1 (size of target S=0)
# N = 136 for SIOSW in S=1 (original size of S=1)
analytic_data %>%
  group_by(S) %>% 
  summarize(
    N_iosw = sum(iosw),
    N_siosw = sum(siosw)
  )

# examine probability of selection distribution by trial
analytic_data %>% 
  group_by(unique_id) %>% 
  summarize(
    mean_ps = mean(ps)
  ) %>% 
  ungroup()

# examine weight distribution by trial
analytic_data %>% 
  group_by(unique_id) %>% 
  summarize(
    mean_w = mean(iosw)
  ) %>% 
  ungroup()

# check for extreme weights
test <- analytic_data %>% 
  filter(iosw > 30)

summary(test[analytic_covs_list])

p <- ggplot(analytic_data, aes(x = ps, fill = as.factor(S))) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Probability of Selection into Multi-Site Trial",
    y = "Density",
    fill = "Multi-site trial",
    title = "Probability of Selection into Multi-Site Trial by Trial Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

p

#ggsave("p_selection_density_plot.png", plot = p, width = 8, height = 6, dpi = 300)

# compare weighted multi-site patients vs unweighted single-site patients
# want SMD < 0.1
covs <- analytic_data %>%
  select(all_of(analytic_covs_list))

tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(iosw_all))

bal_site_iosw <- bal.tab(S ~ covs, data = analytic_data,
                      weights = 'iosw_all',
                      binary = "std", continuous = "std", 
                      stats = c('mean.diffs'),
                      s.d.denom = 'pooled',
                      un = TRUE, 
                      thresholds = c(m = 0.1)
) 

bal_tab_site_iosw <- bal_site_iosw$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_site_iosw)
row.names(bal_tab_site_iosw) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

var_labels <- c(
  "age__years_" = "Maternal age",
  "any_pe" = "Any pre-eclampsia",
  "any_preterm" = "Any pre-term births",
  "any_sga"= "Any SGA neonates",
  "bmi_above_30" = "BMI>30",
  "smoking_bin_yes" = "Smoking",
  "aspirin" = "Aspirin"
)

#write.csv(bal_tab_site_iosw, "bal_tab_site_iosw.csv", row.names = TRUE)

loveplot_site_iosw <- love.plot(
  bal_site_iosw,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Site with IOSW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_site_iosw <- loveplot_site_iosw + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

loveplot_site_iosw

#ggsave("loveplot_site_iosw.png", plot = loveplot_site_iosw, width = 8, height = 6, dpi = 300, bg = "white")
dev.off()

# check SMD by treatment in weighted multi-site patients
wmulti_data <- analytic_data %>% filter(S == 1)
wmulti_covs <- wmulti_data %>%
  select(all_of(analytic_covs_list))

tidy_smd(wmulti_data, analytic_covs_list, .group = trt, .wts = c(iosw_all))

bal_trt_iosw <- bal.tab(trt ~ wmulti_covs, data = wmulti_data,
                      weights = 'iosw_all',
                      binary = "std", continuous = "std", 
                      stats = c('mean.diffs'),
                      s.d.denom = 'pooled',
                      un = TRUE, 
                      thresholds = c(m = 0.1)
) 

bal_tab_trt_iosw <- bal_trt_iosw$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_trt_iosw)
row.names(bal_tab_trt_iosw) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

#write.csv(bal_tab_trt_iosw, "bal_tab_trt_iosw.csv", row.names = TRUE)

loveplot_trt_iosw <- love.plot(
  bal_trt_iosw,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Treatment with IOSW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_trt_iosw <- loveplot_trt_iosw + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

loveplot_trt_iosw
#ggsave("loveplot_trt_iosw.png", plot = loveplot_trt_iosw, width = 8, height = 6, dpi = 300, bg = "white")

# although iosw balances cov distribution between single and multi-site patients
# the weighs also introduce confouding between trt groups in weighted multi-site patients
# so we will compute iptw to balance by treatment in the weighted pseudopopulation

## IPTW IN WEIGHTED MULTI-SITE PATIENTS ##
# fit trt model to estimate probability of trt in iosw-weighted pseudopopulation

wmulti_data <- analytic_data %>% filter(S == 1)

trt_model <- glm(
  trt ~ age__years_ +
    any_pe +
    any_preterm +
    any_sga +
    bmi_above_30 +
    smoking_bin +
    aspirin,
  data = wmulti_data,
  weights = iosw_all,
  family = binomial(link = "logit")
)

wmulti_data$pa <- predict(trt_model, type = "response")

wmulti_data$iptw <- with(wmulti_data, ifelse(trt == 1,
                                             1 / pa,
                                             1 / (1 - pa)))

# stabilized weights
# (same regardless of whether use iosw or siows for pseudo population)
marg_trt_model <- glm(trt ~ 1,
                    data = wmulti_data,
                    family = binomial(link = "logit"))

wmulti_data$marg_pa <- predict(marg_trt_model, type = "response")

wmulti_data$siptw <- with(wmulti_data, ifelse(trt == 1,
                                             iptw * marg_pa,
                                             iptw * (1-marg_pa)
                                             ))

# check N of pseudo populations in S=1, which should be:
# N = 272 for IPTW in S=1 (since 50:50 randomization, twice the size)
# N = 136 for SIPTW in S=1 (original size of S=1)
wmulti_data %>%
  summarize(
    N_iptw = sum(iptw),
    N_siptw = sum(siptw)
  )

analytic_data$iptw <- NA
analytic_data$siptw <- NA
analytic_data$iptw[analytic_data$S == 1] <-wmulti_data$iptw
analytic_data$siptw[analytic_data$S == 1] <-wmulti_data$siptw

## IPTW FOR SINGLE-SITE PATIENTS ## 
single_data <- analytic_data %>% 
  filter(S == 0)

trt_model_single <- glm(trt ~ age__years_ + 
                      any_pe + 
                      any_preterm + 
                      any_sga + 
                      bmi_above_30 + 
                      smoking_bin + 
                      aspirin,
                    data = single_data,
                    family = binomial(link = "logit"))

single_data$pa <- predict(trt_model_single, type = "response")

single_data$iptw <- with(single_data, ifelse(trt == 1,
                                             1 / pa,
                                             1 / (1 - pa)))

# stabilize weights
marg_trt_model <- glm(trt ~ 1, data = single_data, family = binomial())
single_data$marg_pa <- predict(marg_trt_model, type = "response")

single_data$siptw <- with(single_data, ifelse(trt == 1,
                                              iptw * marg_pa,
                                              iptw * (1-marg_pa)))


# check N of pseudo populations, which should be:
# N = 718 for IPTW in S=1 (twice size of S=0)
# N = 359 for SIPTW in S=1 (original size of S=0)
single_data %>%
  summarize(
    N_iptw = sum(iptw),
    N_siptw = sum(siptw)
  )

analytic_data$iptw[analytic_data$S == 0] <-single_data$iptw
analytic_data$siptw[analytic_data$S == 0] <-single_data$siptw

analytic_data %<>%
  mutate(
    w = iptw * iosw, # IOSW of 0 for S0 
    sw = siptw * siosw, # IOSW of 0 for S0 
    w_all = iptw * iosw_all, # IOSW of 1 for S0
    sw_all= siptw * siosw_all # IOSW of 1 for S0
  )

# check trt balance with combined both weights
tidy_smd(analytic_data, analytic_covs_list, .group = trt, .wts = c(w_all))

print(tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(w_all, iosw_all, iptw)), n=40)

# get means/N/% in weighted sample
library(survey)

design <- svydesign(
  ids = ~1,
  data = analytic_data,
  weights = ~w_all
)

results <- lapply(analytic_covs_list, function(var) {
  f <- as.formula(paste0("~", var))
  svyby(f, ~S, design, svymean, na.rm = TRUE)
})

svyby(~aspirin, ~S, design, svytotal, na.rm = TRUE)

data_sub <- analytic_data[analytic_data$S == 1, ]

design_sub <- svydesign(
  ids = ~1,
  data = data_sub,
  weights = ~w_all
)

svymean(~age__years_, design_sub)
svyvar(~age__years_, design_sub) |> sqrt()
svyquantile(~age__years_, design_sub, quantiles = c(0, 0.5, 1), na.rm = TRUE)

bal_trt_w <- bal.tab(trt ~ covs, data = analytic_data,
                   weights = 'w_all',
                   binary = "std", continuous = "std", 
                   stats = c('mean.diffs'),
                   s.d.denom = 'pooled',
                   un = TRUE, 
                   thresholds = c(m = 0.1)
) 

bal_tab_trt_w <- bal_trt_w$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_trt_w)
row.names(bal_tab_trt_w) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

#write.csv(bal_tab_trt_w, "bal_tab_trt_w.csv", row.names = TRUE)

loveplot_trt_w <- love.plot(
  bal_trt_w,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Treatment with IOSW + IPTW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_trt_w <- loveplot_trt_w + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

loveplot_trt_w

#ggsave("loveplot_trt_w.png", plot = loveplot_trt_w, width = 8, height = 6, dpi = 300, bg = "white")

# check site balance with both weights
tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(w_all))

bal_site_w <- bal.tab(S ~ covs, data = analytic_data,
                         weights = 'w_all',
                         binary = "std", continuous = "std", 
                         stats = c('mean.diffs'),
                         s.d.denom = 'pooled',
                         un = TRUE, 
                         thresholds = c(m = 0.1)
) 

bal_tab_site_w <- bal_site_w$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_site_w)
row.names(bal_tab_site_w) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

#write.csv(bal_tab_site_w, "bal_tab_site_w.csv", row.names = TRUE)

loveplot_site_w <- love.plot(
  bal_site_w,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Site with IOSW + IPTW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_site_w <- loveplot_site_w + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

loveplot_site_w
#ggsave("loveplot_site_w.png", plot = loveplot_site_w, width = 8, height = 6, dpi = 300, bg = "white")

# recap of all balance plots
library(gridExtra)

loveplot_grid <- grid.arrange(
  loveplot_site_iosw, 
  loveplot_trt_iosw, 
  loveplot_site_w, 
  loveplot_trt_w, 
  ncol = 2
  )

loveplot_grid
#ggsave("loveplot_grid.png", plot = loveplot_grid, width = 8, height = 6, dpi = 300, bg = "white")

p_weight_distribution <- ggplot(analytic_data, aes(x = w_all, fill = as.factor(S))) +
  geom_density(alpha = 0.3) +
  labs(
    x = "IOSW + IPTW",
    y = "Density",
    fill = "Multi-site trial",
    title = "Weight Distribution (IOSW + IPTW) by Trial Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
p_weight_distribution

#ggsave("p_weight_distribution.png", plot = p_weight_distribution, width = 8, height = 6, dpi = 300, bg = "white")

summary(analytic_data$w)

analytic_data %>% 
  group_by(unique_id) %>% 
  summarize(
    w = mean(w)
  )



#### 10. COMPUTE DIFFERENT ATE ESTIMATORS ####

## A) CRUDE ATE (as a rate difference)

# let's check crude prior to applying restrictions
crude_prior_1 <- sum(data$outcome[data$trt == 1])/sum(data$trt == 1)
crude_prior_0 <- sum(data$outcome[data$trt == 0])/sum(data$trt == 0)
crude_prior <- crude_prior_1 - crude_prior_0

round(crude_prior_1, 4)*100
round(crude_prior_0, 4)*100
round(crude_prior, 4)*100

# and crude in analytic sample with restrictions
crude_1 <- sum(analytic_data$outcome[analytic_data$trt == 1])/sum(analytic_data$trt == 1)
crude_0 <- sum(analytic_data$outcome[analytic_data$trt == 0])/sum(analytic_data$trt == 0)
crude <- crude_1 - crude_0

round(crude_1, 4)*100
round(crude_0, 4)*100
round(crude, 5)*100

## B) OUTCOME MODELING

# model probability of outcome in treated multi-site patients (S1A1)
outcome_model <- as.formula(paste("outcome ~", paste(analytic_covs_list, collapse = "+")))
outcome_model
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
summary(S1data_A1[analytic_covs_list])

OM1mod <- glm(formula = outcome_model, data = S1data_A1, family = "binomial" (link = "logit"))
summary(OM1mod)

# predict probability of outcome in S0 patients based on outcome model in S1A1
S0_data <- analytic_data %>% filter(S==0)
p1 <- predict(OM1mod, newdata = S0_data, type = "response")
S0_data$p1 <- p1
analytic_data$p1[S==0] <- p1

# model probability of outcome in untreated multi-site patients (S1A0)
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
summary(S1data_A0[analytic_covs_list])

OM0mod <- glm(formula = outcome_model, data = S1data_A0, family = "binomial" (link = "logit"))

# predict probability of outcome in whole population based on outcome model in S1A0
p0 <- predict(OM0mod, newdata = S0_data, type = "response")
analytic_data$p0[S==0] <- p0

OM_1 <- mean(S0_data$p1)
OM_0 <- mean(S0_data$p0)
OM <- mean(S0_data$p1) - mean(S0_data$p0)
OM

round(OM_1, 4)*100
round(OM_0, 4)*100
round(OM, 4)*100

## C) INVERSE ODDS WEIGHTING

A <- analytic_data$trt
S <- analytic_data$S
w <- analytic_data$w
sw <- analytic_data$sw # stabilized weights
Y <- analytic_data$outcome

# unstabilized weights
IOW1_1 <- (sum(A*S*w)^-1) * sum(A*S*w*Y)
IOW1_0 <- (sum((1-A)*S*w)^-1) * sum((1-A)*S*w*Y)

IOW1 = IOW1_1 - IOW1_0
IOW1

round(IOW1_1, 4)*100
round(IOW1_0, 4)*100
round(IOW1, 4)*100

# stabilized weights
IOW2_1 <- (sum(A*S*sw)^-1) * sum(A*S*sw*Y)
IOW2_0 <- (sum((1-A)*S*sw)^-1) * sum((1-A)*S*sw*Y)

IOW2 = IOW2_1 - IOW2_0
IOW2

round(IOW2_1, 4)*100
round(IOW2_0, 4)*100
round(IOW2, 4)*100


## D) DOUBLY ROBUST ESTIMATOR
A <- analytic_data$trt
S <- analytic_data$S
w <- analytic_data$w
w_all <- analytic_data$w_all # S0 patients have IOSW of 1 and their respective IPTW weight
sw <- analytic_data$sw # stabilized weights
sw_all <- analytic_data$sw_all
Y <- analytic_data$outcome

# model probability of outcome in treated multi-site patients (S1A1)
outcome_model <- as.formula(paste("outcome ~", paste(analytic_covs_list, collapse = "+")))
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
OM1mod <- glm(formula = outcome_model, data = S1data_A1, family = "binomial" (link = "logit"))
p1 <- predict(OM1mod, newdata = analytic_data, type = "response")
analytic_data$p1 <- p1

S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
OM0mod <- glm(formula = outcome_model, data = S1data_A0, family = "binomial" (link = "logit"))
p0 <- predict(OM0mod, newdata = analytic_data, type = "response")
analytic_data$p0 <- p0

p1 <- analytic_data$p1
p0 <- analytic_data$p0

# using unstabilized weights
# DR1_1 <- (sum(S*A*w)^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
# DR1_0 <- (sum(S*(1-A)*w)^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)

numerator_A1 <- sum(S*A*w_all*Y) - sum(S*A*w_all*p1) + sum((1-S)*A*w_all*Y)
denominator_A1 <- sum(S*A*w_all)
DR1_1 <- numerator_A1/denominator_A1

numerator_A0 <- sum(S*(1-A)*w_all*Y) - sum(S*(1-A)*w_all*p0) + sum((1-S)*A*w_all*Y)
denominator_A0 <- sum(S*(1-A)*w_all)
DR1_0 <- numerator_A0/denominator_A0

DR1 <- DR1_1 - DR1_0
DR1

round(DR1_1, 4)*100
round(DR1_0, 4)*100
round(DR1, 4)*100

# stabilized weights
numerator_sA1 <- sum(S*A*sw_all*Y) - sum(S*A*sw_all*p1) + sum((1-S)*A*sw_all*Y)
denominator_sA1 <- sum(S*A*sw_all)
DR2_1 <- numerator_sA1/denominator_sA1

numerator_sA0 <- sum(S*(1-A)*sw_all*Y) - sum(S*(1-A)*sw_all*p0) + sum((1-S)*A*sw_all*Y)
denominator_sA0 <- sum(S*(1-A)*sw_all)
DR2_0 <- numerator_A0/denominator_sA0

DR2 <- DR2_1 - DR2_0
DR2

round(DR2_1, 4)
round(DR2_0, 4)
round(DR2, 4)

# DR using weighted regression (unstabilized weights)
S0data <- subset(analytic_data, S == 0)

S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
DR1mod_w <- glm(formula = outcome_model, data = S1data_A1, weights = S1data_A1$w, family = "binomial" (link = "logit"))
p1_w <- predict(DR1mod_w, newdata = S0data, type = "response")

S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
DR0mod_w <- glm(formula = outcome_model, data = S1data_A0, weights = S1data_A0$w, family = "binomial" (link = "logit"))
p0_w <- predict(DR0mod_w, newdata = S0data, type = "response")

DR3_1 <- mean(p1_w)
DR3_0 <- mean(p0_w)
DR3 <- DR3_1 - DR3_0
DR3

round(DR3_1, 4)*100
round(DR3_0, 4)*100
round(DR3, 4)*100

# DR stabilized weights
DR1mod_sw <- glm(formula = outcome_model, data = S1data_A1, weights = S1data_A1$sw, family = "binomial" (link = "logit"))
p1_sw <- predict(DR1mod_sw, newdata = S0data, type = "response")

DR0mod_sw <- glm(formula = outcome_model, data = S1data_A0, weights = S1data_A0$sw, family = "binomial" (link = "logit"))
p0_sw <- predict(DR0mod_sw, newdata = S0data, type = "response")

DR4_1 <- mean(p1_sw)
DR4_0 <- mean(p0_sw)
DR4 <- DR4_1 - DR4_0
DR4

round(DR4_1, 4)*100
round(DR4_0, 4)*100
round(DR4, 4)*100

#### 11. FIT WEIGHTED GEE MODEL FOR OUTCOME ####  

## A) CRUDE GEE

# crude from whole population before restrictions

gee_crude_prior_model <- glmgee(
  outcome ~ trt,
  id = unique_id, 
  data = data,
  family = binomial("logit"),
  corstr = "exchangeable" # can try independence, exchangeable, unstructured
)

treated <- data
treated$trt <- 1

untreated <- data
untreated$trt <- 0

treated$pred <- predict(gee_crude_prior_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_crude_prior_model, newdata = untreated, type = "response")

GEE_crude_prior_1 <- mean(treated$pred)
GEE_crude_prior_0 <- mean(untreated$pred)

GEE_crude_prior <- GEE_crude_prior_1 - GEE_crude_prior_0
GEE_crude_prior
round(GEE_crude_prior_1, 4)*100
round(GEE_crude_prior_0, 4)*100
round(GEE_crude_prior, 4)*100
rm(treated, untreated)

# crude from whole population after restrictions

# GEE coefficients have same interpretation as in logistic model
# i.e., odds ratio
# see: https://library.virginia.edu/data/articles/getting-started-with-generalized-estimating-equations 

gee_crude_model <- glmgee(
  outcome ~ trt,
  id = unique_id, 
  data = analytic_data,
  family = binomial("logit"),
  corstr = "exchangeable" # can try independence, exchangeable, unstructured
)

#summary(gee_crude_model)
exp(gee_crude_model$coefficients[2]) # this is the OR

# we got the OR but we want the outcome probability in trt1, trt0, and the ATE
# can follow steps from here:
# https://www.r-bloggers.com/2021/06/estimating-a-risk-difference-and-confidence-intervals-using-logistic-regression/
# alternatively, just use identity link function to get RD directly

treated <- analytic_data
treated$trt <- 1

untreated <- analytic_data
untreated$trt <- 0

treated$pred <- predict(gee_crude_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_crude_model, newdata = untreated, type = "response")

GEE_crude_1 <- mean(treated$pred)
GEE_crude_0 <- mean(untreated$pred)

GEE_crude <- GEE_crude_1 - GEE_crude_0
GEE_crude
round(GEE_crude_1, 4)*100
round(GEE_crude_0, 4)*100
round(GEE_crude, 4)*100
rm(treated, untreated)

## B) IOSW GEE

# only in weighted multi-site (S = 1) patients
# because GEE model cannot handle weight of 0 for single-site patients

analytic_data_S1 <- analytic_data %>% 
  filter(S == 1)

# first use logit to get separate probability by treatment group
gee_w_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable", 
  weights = analytic_data_S1$w
)

treated <- analytic_data_S1
treated$trt <- 1

untreated <- analytic_data_S1
untreated$trt <- 0

treated$pred <- predict(gee_w_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_w_model, newdata = untreated, type = "response")

GEE_IOW1_1 <- mean(treated$pred)
GEE_IOW1_0 <- mean(untreated$pred)
round(GEE_IOW1_1, 4)*100
round(GEE_IOW1_0, 4)*100
round(GEE_IOW1_1 - GEE_IOW1_0, 4)*100

rm(treated, untreated)

# now use identity to get the risk difference directly
# and apply Mancl & deRouen small sample bias variance correction after
gee_w_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "identity"),
  corstr = "exchangeable", 
  weights = analytic_data_S1$w
)

vcov_bias_corrected <- vcov(gee_w_model, type = "bias-corrected")
se_bias_corrected <- sqrt(diag(vcov_bias_corrected))

cbind(
  estimate = round(gee_w_model$coefficients[2]*100, 2), 
  std_error = se_bias_corrected[2], 
  lower = round(gee_w_model$coefficients[2]*100 - se_bias_corrected[2]*1.95, 2), 
  upper = round(gee_w_model$coefficients[2]*100 + se_bias_corrected[2]*1.95, 2)
  )

GEE_IOW1 <- coef(gee_w_model)[2]

# C) STABILZIED IOW GEE

# first use logit to get separate probability by treatment group
gee_sw_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable",
  weights = analytic_data_S1$sw
)

summary(gee_sw_model)
exp(gee_sw_model$coefficients[2]) # this is the OR

treated <- analytic_data_S1
treated$trt <- 1

untreated <- analytic_data_S1
untreated$trt <- 0

treated$pred <- predict(gee_sw_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_sw_model, newdata = untreated, type = "response")

GEE_IOW2_1 <- mean(treated$pred)
GEE_IOW2_0 <- mean(untreated$pred)

rm(treated, untreated)

gee_sw_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "identity"),
  corstr = "exchangeable",
  weights = analytic_data_S1$sw
)

vcov_bias_corrected_sw <- vcov(gee_sw_model, type = "bias-corrected")
se_bias_corrected_sw <- sqrt(diag(vcov_bias_corrected_sw))

cbind(
  Estimate = coef(gee_sw_model),
  `Std.error` = se_bias_corrected_sw
)

cbind(
  estimate = gee_sw_model$coefficients[2]*100, 
  std_error = se_bias_corrected_sw[2], 
  lower = gee_sw_model$coefficients[2]*100 - se_bias_corrected_sw[2]*1.95, 
  upper = gee_sw_model$coefficients[2]*100 + se_bias_corrected_sw[2]*1.95
)

GEE_IOW2 <- coef(gee_sw_model)[2]

#### 12. FORMAT RESULTS ####

results <- data.frame(
  "measure" = c('Proportion Events (A=1)', 'Proportion Events (A=0)', 'ATE'),
  "Crude prior to restrictions" = c(crude_prior_1, crude_prior_0, crude_prior),
  "Crude" = c(crude_1, crude_0, crude),
  "Outcome Model" = c(OM_1, OM_0, OM),
  "IOSW" = c(IOW1_1, IOW1_0, IOW1),
  "Stabilized IOSW" = c(IOW2_1, IOW2_0, IOW2),
  "Doubly Robust" = c(DR1_1, DR1_0, DR1),
  "Doubly Robust (stabilized)" = c(DR2_1, DR2_0, DR2),
  "Doubly Robust (weighted regression)" = c(DR3_1, DR3_0, DR3),
  "Doubly Robust (stabilized weighted regression)" = c(DR4_1, DR4_0, DR4),
  "Crude (GEE)" = c(GEE_crude_1, GEE_crude_0, GEE_crude),
  "IOSW (GEE)"= c(GEE_IOW1_1, GEE_IOW1_0, GEE_IOW1),
  "Stabilized IOSW (GEE)"= c(GEE_IOW2_1, GEE_IOW2_0, GEE_IOW2)
)

results %<>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

results %<>%
  mutate(across(where(is.numeric), ~ .x*100))

# transpose results data frame for easier reading
results_transposed <- as.data.frame(t(results[,-1]))
colnames(results_transposed) <- results$measure
results_transposed$Method <- rownames(results_transposed)
results_transposed <- results_transposed[, c(ncol(results_transposed), 1:(ncol(results_transposed)-1))]
rownames(results_transposed) <- NULL

#write.table (results_transposed , "results_table.csv", col.names = T, row.names=F, append= F, sep=',')

#### 13. ADDITIONAL RESULTS ####

## by subgroups of S=1

S1data <- analytic_data %>% 
  filter(S == 1)

# check if results are the same w/wo stabilization
S1data %<>%
  mutate(
    iosw_outcome= iosw * outcome,
    siosw_outcome = siosw * outcome,
    iptw_outcome = iptw * outcome,
    siptw_outcome = siptw * outcome,
    w_outcome = w * outcome,
    sw_outcome = sw * outcome
  )

S1data %>%
  group_by(trt) %>%
  summarize(
    iosw_Y = sum(iosw_outcome)/sum(iosw),
    siosw_Y = sum(siosw_outcome)/sum(siosw),
    iptw_Y = sum(iptw_outcome)/sum(iptw),
    siptw_Y = sum(siptw_outcome)/sum(siptw),
    w_Y = sum(w_outcome)/sum(w),
    sw_Y = sum(sw_outcome)/sum(sw)
  )

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

S0data <- analytic_data %>% 
  filter(S == 0)

# check if results are the same w/wo stabilization
S0data %<>%
  mutate(
    iosw_outcome= iosw_all * outcome,
    siosw_outcome = siosw_all * outcome,
    iptw_outcome = iptw * outcome,
    siptw_outcome = siptw * outcome,
    w_outcome = w_all * outcome,
    sw_outcome = sw_all * outcome
  )

S0data %>%
  group_by(trt) %>%
  summarize(
    iosw_Y = sum(iosw_outcome)/sum(iosw_all),
    siosw_Y = sum(siosw_outcome)/sum(siosw_all),
    iptw_Y = sum(iptw_outcome)/sum(iptw),
    siptw_Y = sum(siptw_outcome)/sum(siptw),
    w_Y = sum(w_outcome)/sum(w_all),
    sw_Y = sum(sw_outcome)/sum(sw_all)
  )

A <- S0data$trt
S <- S0data$S
Y <- S0data$outcome
iosw <- S0data$iosw_all
siosw <- S0data$siosw_all
iptw <- S0data$iptw
siptw <- S0data$siptw
w <- S0data$w_all
sw <- S0data$sw_all

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
A <- analytic_data$trt
S <- analytic_data$S
Y <- analytic_data$outcome
iosw <- analytic_data$iosw_all
siosw <- analytic_data$siosw_all
iptw <- analytic_data$iptw
siptw <- analytic_data$siptw
w <- analytic_data$w_all
sw <- analytic_data$sw_all

iosw_1 <- (sum(iosw)^-1) * sum(A*S*iosw*Y)
iosw_0 <- (sum(iosw)^-1) * sum((1-A)*S*iosw*Y)
iosw = iosw_1 - iosw_0

siosw_1 <- (sum(siosw)^-1) * sum(A*S*siosw*Y)
siosw_0 <- (sum(siosw)^-1) * sum((1-A)*S*siosw*Y)
siosw = siosw_1 - siosw_0

iptw_1 <- (sum(iptw)^-1) * sum(A*S*iptw*Y)
iptw_0 <- (sum(iptw)^-1) * sum((1-A)*S*iptw*Y)
iptw = iptw_1 - iptw_0

siptw_1 <- (sum(siptw)^-1) * sum(A*S*siptw*Y)
siptw_0 <- (sum(siptw)^-1) * sum((1-A)*S*siptw*Y)
siptw = siptw_1 - siptw_0

w_1 <- (sum(w)^-1) * sum(A*S*w*Y)
w_0 <- (sum(w)^-1) * sum((1-A)*S*w*Y)
w = w_1 - w_0

sw_1 <- (sum(sw)^-1) * sum(A*S*sw*Y)
sw_0 <- (sum(sw)^-1) * sum((1-A)*S*sw*Y)
sw = sw_1 - sw_0

# format additional results

additional_results <- data.frame(
  "Group and weight" = c('Proportion Events (A=1)', 'Proportion Events (A=0)', 'ATE'),  
  "S1 IOSW" = c(S1_iosw_1, S1_iosw_0, S1_iosw),
  "S1 SIOSW" = c(S1_siosw_1, S1_siosw_0, S1_siosw),
  "S1 IPTW" = c(S1_iptw_1, S1_iptw_0, S1_iptw),
  "S1 SIPTW" = c(S1_siptw_1, S1_siptw_0, S1_siptw),
  "S1 combo w" = c(S1_w_1, S1_w_0, S1_w),
  "S1 combo sw" = c(S1_sw_1, S1_sw_0, S1_sw),
  "S0 IOSW" = c(S0_iosw_1, S0_iosw_0, S0_iosw),
  "S0 SIOSW" = c(S0_siosw_1, S0_siosw_0, S0_siosw),
  "S0 IPTW" = c(S0_iptw_1, S0_iptw_0, S0_iptw),
  "S0 SIPTW" = c(S0_siptw_1, S0_siptw_0, S0_siptw),
  "S0 combo w" = c(S0_w_1, S0_w_0, S0_w),
  "S0 combo sw" = c(S0_sw_1, S0_sw_0, S0_sw),
  "all IOSW" = c(iosw_1, iosw_0, iosw),
  "all SIOSW" = c(siosw_1, siosw_0, siosw),
  "all IPTW" = c(iptw_1, iptw_0, iptw),
  "all SIPTW" = c(siptw_1, siptw_0, siptw),
  "all combo w" = c(w_1, w_0, w),
  "all combo sw" = c(sw_1, sw_0, sw)
)

additional_results %<>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

additional_results %<>%
  mutate(across(where(is.numeric), ~ .x*100))

additional_results_transposed <- as.data.frame(t(additional_results[,-1]))
colnames(additional_results_transposed) <- additional_results$Group.and.weight
additional_results_transposed$Group.and.weight <- rownames(additional_results_transposed)
additional_results_transposed <- additional_results_transposed[, c(ncol(additional_results_transposed), 1:(ncol(additional_results_transposed)-1))]
rownames(additional_results_transposed) <- NULL

#write.table (additional_results_transposed , "additional_results_table.csv", col.names = T, row.names=F, append= F, sep=',')



