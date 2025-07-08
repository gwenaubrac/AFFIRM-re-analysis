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
## ---------------------------

#### 1. LOAD PACKAGES ####

options(scipen = 999)

library(PSweight)
library(dplyr)
library(magrittr)
library(table1)
library(writexl)
library(tidysmd)
library(medicaldata)
library(cobalt)
library(ggplot2)
library(glmtoolbox)

#### 2. READ AND EXPLORE DATA ####

setwd("Z:/RPLATT/GWEN-OHRI")
data <- read.csv("AFFIRM_GWEN.csv")

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
  "other_thrombo_philia",
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

write.table (tab1 , "./results/table1_by_site_prior.csv", col.names = T, row.names=F, append= F, sep=',')

caption <- 'Baseline Characteristics of Cohort by Trial'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|unique_id"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1

write.table (tab1 , "./results/table1_by_trial_prior.csv", col.names = T, row.names=F, append= F, sep=',')


caption <- 'Baseline Characteristics of Cohort by Treatment'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|trt"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1

write.table (tab1 , "./results/table1_by_trt_prior.csv", col.names = T, row.names=F, append= F, sep=',')



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

write.table (tab2 , "./results/table1_by_trt_clean_covs.csv", col.names = T, row.names=F, append= F, sep=',')

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
    GA_at_start_cat == "16 to <20 weeks" |
    GA_at_start_cat == ">=20 weeks"
  )

analytic_data %<>%
  filter(!new_id %in% to_exclude$new_id)

# list of covariates adjusting for in models
analytic_covs_list <- c(
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
  "aspirin"
)

# covariates by single vs multi-site trials after restrictions
caption <- "Baseline Characteristics by Site After Restrictions"
tab3_formula <- as.formula(paste("~", paste(analytic_covs_list, collapse = "+"), "|single"))
tab3 <- table1(
  tab3_formula,
  data = analytic_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab3

write.table (tab3 , "./results/table1_by_site_restr.csv", col.names = T, row.names=F, append= F, sep=',')

# covariates by trial after applying exclusions
caption <- "Baseline Characteristcs by Trial After Restrictions"
tab3_formula <- as.formula(paste("~", paste(analytic_covs_list, collapse = "+"), "|unique_id"))
tab3 <- table1(
  tab3_formula,
  data = analytic_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab3
write.table (tab3 , "./results/table1_by_trial_restr.csv", col.names = T, row.names=F, append= F, sep=',')



#### 9. COMPUTE IOSW WEIGHTS AND CHECK BALANCE ####

# adapted from Dahabreh et al.
# thinking of S as "selection into multi-site trial", so:
# single site patients (single = 1 -> S = 0) -- this is our "target" population
# multi site patients (single = 0 -> S = 1) is our "study" population
# so our goal is to extend inference from multi-site to single-site patients

analytic_data %<>%
  mutate(S = if_else(single == 1, 0, 1))

analytic_covs_list <- c(
  "age__years_",
  "any_pe",
  "any_preterm",
  "any_sga",
  "bmi_above_30",
  "smoking_bin",
  "aspirin"
)

# originally we tried to incorporate treatment balancing directly in the IOSW weights (pa)
# but we struggled to achieve balance

# selection_model <- as.formula(paste("S ~", paste(analytic_covs_list, collapse = "+")))
# treatment_model <- as.formula(paste("trt ~", paste(analytic_covs_list, collapse = "+")))

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
S1data <- subset(analytic_data, S == 1)
summary(S1data[analytic_covs_list])

w_reg <- glm(selection_model, family = binomial(), data = analytic_data)
summary(w_reg)
ps <- predict(w_reg, newdata = analytic_data, type = "response")
analytic_data$ps <- ps

# w_reg2 <- glm(treatment_model, family = "binomial", data = S1data)
# summary(w_reg2)
# pa <- predict(w_reg2, newdata = analytic_data, type = "response")
# analytic_data$pa <- pa

# iosw = (analytic_data$trt * analytic_data$S * (1-ps))/(ps*pa) + ((1-analytic_data$trt) * analytic_data$S*(1-ps)/(ps*(1-pa)))
# analytic_data$iosw <- iosw

iosw = analytic_data$S * (1-ps)/ps
analytic_data$iosw <- iosw

# stabilize weights: multiply by marginal probability of selection
marg_s_model <- glm(S ~ 1, family = binomial(), data = analytic_data)
marg_ps <- predict(marg_s_model, newdata = analytic_data, type = "response")
analytic_data$marg_ps <- marg_ps
analytic_data %<>%
  mutate(siosw = iosw * marg_ps)


summary(analytic_data$ps)
#summary(analytic_data$pa)
summary(analytic_data$iosw)
sd(analytic_data$iosw)
summary(analytic_data$siosw)
sd(analytic_data$siosw)

analytic_data %>%
  group_by(S) %>% 
  summarize(
    N = sum(iosw)
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


ggsave("./results/p_selection_density_plot.png", plot = p, width = 8, height = 6, dpi = 300)

# compare weighted multi-site patients vs unweighted single-site patients
# want SMD < 0.1

# let's update iosw to 1 to patients with S = 0 (rather than 0 as is currently)
# so using iosw means looking at unweighted single-site (S = 0) + weighted multi-site (S = 1) together

covs <- analytic_data %>%
  select(all_of(analytic_covs_list))

analytic_data %<>%
  mutate(
    iosw = if_else(S == 0, 1, iosw),
    siosw = if_else(S == 0, 1, siosw)
  )

tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(iosw))

bal_site_iosw <- bal.tab(S ~ covs, data = analytic_data,
                      weights = 'iosw',
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

write.csv(bal_tab_site_iosw, "./results/bal_tab_site_iosw.csv", row.names = TRUE)

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

ggsave("./results/loveplot_site_iosw.png", plot = loveplot_site_iosw, width = 8, height = 6, dpi = 300, bg = "white")
dev.off()

# check SMD by treatment in weighted multi-site patients
wmulti_data <- analytic_data %>% filter(S == 1)
wmulti_covs <- wmulti_data %>%
  select(all_of(analytic_covs_list))

tidy_smd(wmulti_data, analytic_covs_list, .group = trt, .wts = c(iosw))

bal_trt_iosw <- bal.tab(trt ~ wmulti_covs, data = wmulti_data,
                      weights = 'iosw',
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

write.csv(bal_tab_trt_iosw, "./results/bal_tab_trt_iosw.csv", row.names = TRUE)

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
ggsave("./results/loveplot_trt_iosw.png", plot = loveplot_trt_iosw, width = 8, height = 6, dpi = 300, bg = "white")

# although iosw balances cov distribution between single and multi-site patients
# the weighs also introduce confouding between trt groups in weighted multi-site patients
# so we will compute iptw to balance by treatment in the weighted pseudopopulation

## IPTW IN WEIGHTED MULTI-SITE PATIENTS ##
library(survey)
wmulti_data <- analytic_data %>% filter(S == 1)

# specify "survey" design where "sampling probability" is iosw weight
design <- svydesign(ids = ~1, data = wmulti_data, weights = ~iosw)

# fit trt model to estimate probability of trt in iosw-weighted pseudopopulation
trt_model <- svyglm(trt ~ age__years_ + 
                          any_pe + 
                          any_preterm + 
                          any_sga + 
                          bmi_above_30 + 
                          smoking_bin + 
                          aspirin,
                          design = design,
                          family = binomial())

wmulti_data$pa <- predict(trt_model, type = "response")

wmulti_data$iptw <- with(wmulti_data, ifelse(trt == 1,
                                             1 / pa,
                                             1 / (1 - pa)))

# normalized weights
# p_trt <- weighted.mean(wmulti_data$trt, w = wmulti_data$iosw)
# 
# wmulti_data$siptw <- with(wmulti_data, ifelse(trt == 1,
#                                               p_trt / pa,
#                                               (1 - p_trt) / (1 - pa)))

# stabilized weights
marg_trt_model <- svyglm(trt ~ 1,
                    design = design,
                    family = binomial())

wmulti_data$marg_pa <- predict(marg_trt_model, type = "response")


wmulti_data$siptw <- with(wmulti_data, ifelse(trt == 1,
                                             marg_pa / pa,
                                             (1-marg_pa) / (1 - pa)))

summary(wmulti_data$iptw)
sd(wmulti_data$iptw)
summary(wmulti_data$siptw)
sd(wmulti_data$siptw)

analytic_data$combo_w <- NA
analytic_data$combo_sw <- NA
analytic_data$iptw <- NA
analytic_data$siptw <- NA
analytic_data$iptw[analytic_data$S == 1] <-wmulti_data$iptw
analytic_data$siptw[analytic_data$S == 1] <-wmulti_data$siptw
analytic_data$combo_w[analytic_data$S == 1] <-wmulti_data$iosw * wmulti_data$iptw
analytic_data$combo_sw[analytic_data$S == 1] <-wmulti_data$siosw * wmulti_data$siptw

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
                    family = binomial())

single_data$pa <- predict(trt_model_single, type = "response")

single_data$iptw <- with(single_data, ifelse(trt == 1,
                                             1 / pa,
                                             1 / (1 - pa)))

# normalize weights
# p_trt <- weighted.mean(single_data$trt, w = single_data$iosw)
# 
# single_data$siptw <- with(single_data, ifelse(trt == 1,
#                                               p_trt / pa,
#                                               (1 - p_trt) / (1 - pa)))
# 
# analytic_data$combo_w[analytic_data$S==0] <- single_data$siptw

# stabilize weights
marg_trt_model <- glm(trt ~ 1, data = single_data, family = binomial())
single_data$marg_pa <- predict(marg_trt_model, type = "response")


single_data$siptw <- with(single_data, ifelse(trt == 1,
                                              marg_pa / pa,
                                              (1-marg_pa) / (1 - pa)))

summary(single_data$iptw)
sd(single_data$iptw)
summary(single_data$siptw)
sd(single_data$siptw)

analytic_data$iptw[analytic_data$S == 0] <-single_data$iptw
analytic_data$siptw[analytic_data$S == 0] <-single_data$siptw
analytic_data$combo_w[analytic_data$S == 0] <-single_data$iosw * single_data$iptw
analytic_data$combo_sw[analytic_data$S == 0] <-single_data$siosw * single_data$siptw

# check trt balance with combined both weights
tidy_smd(analytic_data, analytic_covs_list, .group = trt, .wts = c(combo_w))

bal_trt_combo_w <- bal.tab(trt ~ covs, data = analytic_data,
                   weights = 'combo_w',
                   binary = "std", continuous = "std", 
                   stats = c('mean.diffs'),
                   s.d.denom = 'pooled',
                   un = TRUE, 
                   thresholds = c(m = 0.1)
) 

bal_tab_trt_combo_w <- bal_trt_combo_w$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_trt_combo_w)
row.names(bal_tab_trt_combo_w) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

write.csv(bal_tab_trt_combo_w, "./results/bal_tab_trt_combo_w.csv", row.names = TRUE)

loveplot_trt_combo_w <- love.plot(
  bal_trt_combo_w,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Treatment with IOSW + IPTW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_trt_combo_w <- loveplot_trt_combo_w + theme(
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

loveplot_trt_combo_w

ggsave("./results/loveplot_trt_combo_w.png", plot = loveplot_trt_combo_w, width = 8, height = 6, dpi = 300, bg = "white")

# check site balance with both weights
tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(combo_w))

bal_site_combo_w <- bal.tab(S ~ covs, data = analytic_data,
                         weights = 'combo_w',
                         binary = "std", continuous = "std", 
                         stats = c('mean.diffs'),
                         s.d.denom = 'pooled',
                         un = TRUE, 
                         thresholds = c(m = 0.1)
) 

bal_tab_site_combo_w <- bal_site_combo_w$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

# clean printed variable names
row.names(bal_tab_site_combo_w)
row.names(bal_tab_site_combo_w) <- c("Maternal age", "Any pre-eclampsia", "Any pre-term births", "Any SGA neonates", "BMI>30", "Smoking", "Aspirin Use")

write.csv(bal_tab_site_combo_w, "./results/bal_tab_site_combo_w.csv", row.names = TRUE)

loveplot_site_combo_w <- love.plot(
  bal_site_combo_w,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Balance by Site with IOSW + IPTW',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2,
  var.names = var_labels
)

loveplot_site_combo_w <- loveplot_site_combo_w + theme(
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

loveplot_site_combo_w
ggsave("./results/loveplot_site_combo_w.png", plot = loveplot_site_combo_w, width = 8, height = 6, dpi = 300, bg = "white")

# recap of all balance plots
library(gridExtra)

loveplot_grid <- grid.arrange(
  loveplot_site_iosw, 
  loveplot_trt_iosw, 
  loveplot_site_combo_w, 
  loveplot_trt_combo_w, 
  ncol = 2
  )

loveplot_grid
ggsave("./results/loveplot_grid.png", plot = loveplot_grid, width = 8, height = 6, dpi = 300, bg = "white")

p_weight_distribution <- ggplot(analytic_data, aes(x = combo_w, fill = as.factor(S))) +
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

ggsave("./results/p_weight_distribution.png", plot = p_weight_distribution, width = 8, height = 6, dpi = 300, bg = "white")

summary(analytic_data$combo_w)

analytic_data %>% 
  group_by(unique_id) %>% 
  summarize(
    combo_w = mean(combo_w)
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

# model probability of outcome in treated multi-site patients
outcome_model <- as.formula(paste("outcome ~", paste(analytic_covs_list, collapse = "+")))
outcome_model
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
summary(S1data_A1[analytic_covs_list])

OM1mod <- glm(formula = outcome_model, data = S1data_A1, family = "binomial" (link = "logit"))
summary(OM1mod)

# predict probability of outcome in whole population based on treated
p1 <- predict(OM1mod, newdata = analytic_data, type = "response")
analytic_data$p1 <- p1

# model probability of outcome in untreated multi-site patients
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
summary(S1data_A0[analytic_covs_list])

OM0mod <- glm(formula = outcome_model, data = S1data_A0, family = "binomial" (link = "logit"))

# predict probability of outcome in whole population based on untreated
p0 <- predict(OM0mod, newdata = analytic_data, type = "response")
analytic_data$p0 <- p0


S0sub <- subset(analytic_data, S == 0)
OM_1 <- mean(S0sub$p1)
OM_0 <- mean(S0sub$p0)
OM <- mean(S0sub$p1) - mean(S0sub$p0)
OM

round(OM_1, 4)*100
round(OM_0, 4)*100
round(OM, 4)*100

## C) INVERSE ODDS WEIGHTING

# unstabilized IOW
A <- analytic_data$trt
S <- analytic_data$S
w <- analytic_data$combo_w
Y <- analytic_data$outcome

# original version: divide by N0 (sample size in S = 0, or single-site patients) - A
IOW1a_1 <- (sum((1-S))^-1) * sum(A*S*w*Y)
IOW1a_0 <- (sum((1-S))^-1) * sum((1-A)*S*w*Y)

IOW1a = IOW1a_1 - IOW1a_0
IOW1a

round(IOW1a_1, 4)*100
round(IOW1a_0, 4)*100
round(IOW1a, 4)*100

# updated version: divide by sum of weights - B
IOW1b_1 <- (sum(w)^-1) * sum(A*S*w*Y)
IOW1b_0 <- (sum(w)^-1) * sum((1-A)*S*w*Y)

IOW1b = IOW1b_1 - IOW1b_0
IOW1b

round(IOW1b_1, 4)*100
round(IOW1b_0, 4)*100
round(IOW1b, 4)*100

# stabilized IOW
S0data <- subset(analytic_data, S == 0)
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
IOW1mod <- glm(formula = outcome ~ 1, data = S1data_A1, weights = S1data_A1$combo_w)
p1 <- predict(IOW1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
IOW0mod <- glm(formula = outcome ~ 1, data = S1data_A0, weights = S1data_A0$combo_w)
p0 <- predict(IOW0mod, newdata = S0data, type = "response")
IOW2_1 <- mean(p1)
IOW2_0 <- mean(p0)
IOW2 <- mean(p1) - mean(p0)
IOW2

round(IOW2_1, 4)*100
round(IOW2_0, 4)*100
round(IOW2, 4)*100

# normalize by dividing each weight by mean weight
analytic_data$w_norm <- with(analytic_data, ave(combo_w, trt, FUN = function(x) x / mean(x)))
w_norm <- analytic_data$w_norm

IOW1c_1 <- (sum((1-S))^-1) * sum(A*S*w_norm*Y)
IOW1c_0 <- (sum((1-S))^-1) * sum((1-A)*S*w_norm*Y)

IOW1c = IOW1c_1 - IOW1c_0
IOW1c

round(IOW1c_1, 4)*100
round(IOW1c_0, 4)*100
round(IOW1c, 4)*100

IOW1d_1 <- (sum(w_norm)^-1) * sum(A*S*w_norm*Y)
IOW1d_0 <- (sum(w_norm)^-1) * sum((1-A)*S*w_norm*Y)

IOW1d = IOW1d_1 - IOW1d_0
IOW1d

round(IOW1d_1, 4)*100
round(IOW1d_0, 4)*100
round(IOW1d, 4)*100

# stabilize using weighted mean (same as above)
# S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
# S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
# 
# IOW2_1 <- weighted.mean(S1data_A1$outcome, S1data_A1$combo_w)
# IOW2_0 <- weighted.mean(S1data_A0$outcome, S1data_A0$combo_w)
# IOW2 <- IOW2_1 - IOW2_0
# IOW2*100

# stabilize manually (same as above)
# w1 <- S1data_A1$iosw
# w0 <- S1data_A0$iosw
# y1 <- S1data_A1$outcome
# y0 <- S1data_A0$outcome
# 
# IOW2_1 <- sum(w1 * y1) / sum(w1)
# IOW2_0 <- sum(w0 * y0) / sum(w0)
# IOW2 <- IOW2_1 - IOW2_0
# IOW2



## D) DOUBLY ROBUST ESTIMATOR

# DR unstabilized
A <- analytic_data$trt
S <- analytic_data$S
Y <- analytic_data$outcome
p1 <- analytic_data$p1
p0 <- analytic_data$p0
w <- analytic_data$combo_w

# original version: divide by N0 (sample size in S = 0, or single-site patients) - A

DR1a_1 <- (sum((1-S))^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
DR1a_0 <- (sum((1-S))^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)

DR1a <- DR1a_1 - DR1a_0
DR1a

round(DR1a_1, 4)*100
round(DR1a_0, 4)*100
round(DR1a, 4)*100

# updated version: divide by sum of weights - B

DR1b_1 <- (sum(w)^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
DR1b_0 <- (sum(w)^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)

DR1b <- DR1b_1 - DR1b_0
DR1b

round(DR1b_1, 4)*100
round(DR1b_0, 4)*100
round(DR1b, 4)*100

# DR using stabilized IOW
sum1_DR2 <- sum(S*A*w*(Y-p1))
sum0_DR2 <- sum(S*(1-A)*w*(Y-p0))
norm1 <- (sum(S*A*w))^-1
norm0 <- (sum(S*(1-A)*w))^-1
DR2_1 <- norm1 * sum1_DR2 + (sum(1-S)^-1) * sum((1-S) * p1)
DR2_0 <- norm0 * sum0_DR2 + (sum(1-S)^-1) * sum((1-S) * p0)
DR2 <- DR2_1 - DR2_0
DR2

round(DR2_1, 4)
round(DR2_0, 4)
round(DR2, 4)

# DR using stabilized IOW w_norm directly

DR1c_1 <- (sum((1-S))^-1) * sum(S*A*w_norm*(Y-p1) + (1-S) * p1)
DR1c_0 <- (sum((1-S))^-1) * sum(S*(1-A)*w_norm*(Y-p0) + (1-S) * p0)

DR1c <- DR1c_1 - DR1c_0
DR1c

DR1d_1 <- (sum(w_norm)^-1) * sum(S*A*w_norm*(Y-p1) + (1-S) * p1)
DR1d_0 <- (sum(w_norm)^-1) * sum(S*(1-A)*w_norm*(Y-p0) + (1-S) * p0)

DR1d <- DR1d_1 - DR1d_0
DR1d

# DR using weighted regression
S0data <- subset(analytic_data, S == 0)
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
DR1mod <- glm(formula = outcome_model, data = S1data_A1, weights = S1data_A1$combo_w, family = "binomial" (link = "logit"))
p1 <- predict(DR1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
DR0mod <- glm(formula = outcome_model, data = S1data_A0, weights = S1data_A0$combo_w, family = "binomial" (link = "logit"))
p0 <- predict(DR0mod, newdata = S0data, type = "response")
DR3_1 <- mean(p1)
DR3_0 <- mean(p0)
DR3 <- DR3_1 - DR3_0
DR3

round(DR3_1, 4)*100
round(DR3_0, 4)*100
round(DR3, 4)*100

#### 11. FIT WEIGHTED GEE MODEL FOR OUTCOME ####  

# first we arrange data by id and grouping var (trial)
analytic_data %<>% 
  arrange(new_id, unique_id)

## A) CRUDE GEE

# crude from whole population

# GEE coefficients have same interpretation as in logistic model
# i.e., odds ratio
# see: https://library.virginia.edu/data/articles/getting-started-with-generalized-estimating-equations 

gee_crude_model <- glmgee(
  outcome ~ trt,
  id = unique_id, 
  data = analytic_data,
  family = binomial("logit"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  type = "bias-corrected" # use sample sample cov matrix
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
gee_weighted_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  weights = analytic_data_S1$combo_w
)

summary(gee_weighted_model)
exp(gee_weighted_model$coefficients[2]) # this is the OR

treated <- analytic_data_S1
treated$trt <- 1

untreated <- analytic_data_S1
untreated$trt <- 0

treated$pred <- predict(gee_weighted_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_weighted_model, newdata = untreated, type = "response")

GEE_IOW1_1 <- mean(treated$pred)
GEE_IOW1_0 <- mean(untreated$pred)
round(GEE_IOW1_1, 4)*100
round(GEE_IOW1_0, 4)*100
round(GEE_IOW1_1 - GEE_IOW1_0, 4)*100

rm(treated, untreated)

# now use identity to get the risk difference directly
# and apply Mancl & deRouen small sample bias variance correction after
gee_weighted_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "identity"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  weights = analytic_data_S1$combo_w
)

vcov_bias_corrected <- vcov(gee_weighted_model, type = "bias-corrected")
se_bias_corrected <- sqrt(diag(vcov_bias_corrected))

cbind(
  estimate = gee_weighted_model$coefficients[2]*100, 
  std_error = se_bias_corrected[2], 
  lower = gee_weighted_model$coefficients[2]*100 - se_bias_corrected[2]*1.95, 
  upper = gee_weighted_model$coefficients[2]*100 + se_bias_corrected[2]*1.95
  )

GEE_IOW1 <- coef(gee_weighted_model)[2]

# C) STABILZIED IOW GEE - DEFINITELY SOMETHING WRONG HERE

analytic_data_S1_norm <- analytic_data_S1
analytic_data_S1_norm$w_norm <- with(analytic_data_S1_norm, ave(iosw, trt, FUN = function(x) x / mean(x)))

# first use logit to get separate probability by treatment group
gee_weighted_model_norm <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  weights = analytic_data_S1_norm$w_norm
)

summary(gee_weighted_model_norm)
exp(gee_weighted_model_norm$coefficients[2]) # this is the OR

treated <- analytic_data_S1
treated$trt <- 1

untreated <- analytic_data_S1
untreated$trt <- 0

treated$pred <- predict(gee_weighted_model_norm, newdata = treated, type = "response")
untreated$pred <- predict(gee_weighted_model_norm, newdata = untreated, type = "response")

GEE_IOW2_1 <- mean(treated$pred)
GEE_IOW2_0 <- mean(untreated$pred)

rm(treated, untreated)

gee_weighted_model_norm <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1_norm,
  id = unique_id,
  family = binomial(link = "identity"),
  corstr = "exchangeable",
  weights = analytic_data_S1_norm$w_norm
)

vcov_bias_corrected_norm <- vcov(gee_weighted_model_norm, type = "bias-corrected")
se_bias_corrected_norm <- sqrt(diag(vcov_bias_corrected_norm))

cbind(
  Estimate = coef(gee_weighted_model_norm),
  `Std.error` = se_bias_corrected_norm
)

GEE_IOW2 <- coef(gee_weighted_model_norm)[2]

#### 12. FORMAT RESULTS ####

results <- data.frame(
  "measure" = c('Proportion Events (A=1)', 'Proportion Events (A=0)', 'ATE'),
  "Crude" = c(crude_1, crude_0, crude),
  "Outcome Model" = c(OM_1, OM_0, OM),
  "IOSW (sumN0)" = c(IOW1a_1, IOW1a_0, IOW1a),
  "IOSW (sumW)" = c(IOW1b_1, IOW1b_0, IOW1b),
  "Stabilized IOSW" = c(IOW2_1, IOW2_0, IOW2),
  "IOSW w_norm (sumN0)" = c(IOW1c_1, IOW1c_0, IOW1c),
  "IOSW w_norm (sumW)" = c(IOW1d_1, IOW1d_0, IOW1d),
  "Doubly Robust (sumN0)" = c(DR1a_1, DR1a_0, DR1a),
  "Doubly Robust (sumW)" = c(DR1b_1, DR1b_0, DR1b),
  "Doubly Robust (stabilized)" = c(DR2_1, DR2_0, DR2),
  "Doubly Robust (sumN0)" = c(DR1c_1, DR1c_0, DR1c),
  "Doubly Robust (sumW)" = c(DR1d_1, DR1d_0, DR1d),
  "Doubly Robust (weighted regression)" = c(DR3_1, DR3_0, DR3),
  "Crude (GEE)" = c(GEE_crude_1, GEE_crude_0, GEE_crude),
  "IOSW (GEE)"= c(GEE_IOW1_1, GEE_IOW1_0, GEE_IOW1),
  "Stabilized IOSW (GEE)"= c(GEE_IOW2_1, GEE_IOW2_0, GEE_IOW2)
)

results %<>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# transpose results data frame for easier reading
results_transposed <- as.data.frame(t(results[,-1]))
colnames(results_transposed) <- results$measure
results_transposed$Method <- rownames(results_transposed)
results_transposed <- results_transposed[, c(ncol(results_transposed), 1:(ncol(results_transposed)-1))]
rownames(results_transposed) <- NULL

