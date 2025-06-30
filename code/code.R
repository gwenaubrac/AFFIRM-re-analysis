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

# first we have to create some variables to help us define the outcome

# determine gestational age at PE onset using trial start and onset dates
data$active_trial_start <- as.Date(data$active_trial_start, format = "%d%b%Y")
data$pre_eclampsia_onset <- as.Date(data$pre_eclampsia_onset, format = "%d%b%Y")

data$PE_days_elapsed <- as.numeric(data$pre_eclampsia_onset - data$active_trial_start)
data$PE_weeks_elapsed <- round(data$PE_days_elapsed/7, 2)
data$GA_at_PE <- data$ga_at_start__weeks_ + data$PE_weeks_elapsed

# define early PE as PE occurring prior to 34 weeks GA
data %<>%
  mutate(early_PE = if_else(pre_eclampsia == 2 & GA_at_PE < 34, 1, 0))

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
# namely those with loss prior to 20 weeks or termination unrelated to outcome
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
    any = sum(
      (pregnancy_outcome == 2 & GA_at_delivery < 20 & !is.na(ga_at_start__weeks_) & !is.na(pregnancy_outcome_date)) |
        (pregnancy_outcome == 3 & GA_at_delivery < 20 & !is.na(ga_at_start__weeks_) & !is.na(pregnancy_outcome_date)) |
        (pregnancy_outcome == 4 & !is.na(ga_at_start__weeks_) & !is.na(pregnancy_outcome_date)) |
        (pregnancy_outcome == 5 & !is.na(ga_at_start__weeks_) & !is.na(pregnancy_outcome_date)) |
        (is.na(ga_at_start__weeks_) | is.na(pregnancy_outcome_date))
    ),
    any_test = sum(explained + unexplained + social + medical + missing)
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
    missing_early_PE = sum(is.na(early_PE) | early_PE == 99),
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


table(data$early_PE)
table(data$severe_pre_eclampsia)
table(data$severe_pre_eclampsia, data$early_PE) # some women have both severe and early PE
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

# write.table (tab1 , "tab1.csv", col.names = T, row.names=F, append= F, sep=',')

caption <- 'Baseline Characteristics of Cohort by Trial'
tab1_formula <- as.formula(paste("~", paste(tab1_covs, collapse = "+"), "|unique_id"))
tab1 <- table1(
  tab1_formula,
  data = tab1_data,
  overall = c(right = 'Total'),
  caption = caption
)
tab1




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
    misS_sga10 = sum(is.na(sga__10) | sga__10 == 99),
    miss_sga3 = sum(is.na(sga__3) | sga__3 == 99),
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

data %>% 
  summarize(
    misS_BMI = sum(is.na(bmi)),
    miss_race = sum(is.na(race) | race == ""),
    miss_paternity = sum(is.na(paternity) | paternity == 99),
    miss_loss_12 = sum(is.na(hx_late_loss_12_weeks) | hx_late_loss_12_weeks == 99),
    miss_loss_16 = sum(is.na(hx_late_loss_16_weeks) | hx_late_loss_16_weeks == 99),
    misS_sga10 = sum(is.na(sga__10) | sga__10 == 99),
    miss_sga3 = sum(is.na(sga__3) | sga__3 == 99),
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
data$smoking[is.na(data$smoking) == T] <- 1
data$vte[is.na(data$vte) == T] <- 1
data$hx_late_loss_20[is.na(data$hx_late_loss_20) == T] <- 0
data$hx_sga_5[is.na(data$hx_sga_5) == T] <- 1
data$hx_sga_3[is.na(data$hx_sga_3) == T] <- 1
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

data <- data %>%
  mutate(
    bmi = case_when(
      is.na(bmi) ~ NA,
      bmi < 18.5 ~ "<18.5",
      bmi >= 18.5 & bmi < 24.9 ~ "18.5-24.9",
      bmi >= 25 & bmi < 29.9 ~ "25-29.9",
      bmi >= 30 & bmi < 34.9 ~ "30-34.9",
      bmi >= 34.9 & bmi <39.9 ~ "34.9-39.9",
      bmi >=40 ~ ">=40"
    ),
    ga_at_start__weeks_ = case_when(
      ga_at_start__weeks_ < 10 ~ "<10 weeks",
      ga_at_start__weeks_ >= 10 & ga_at_start__weeks_ < 16 ~ "10 to <16 weeks",
      ga_at_start__weeks_ >= 16 & ga_at_start__weeks_< 20 ~ "16 to <20 weeks",
      ga_at_start__weeks_ >= 20 ~ ">=20 weeks"
    ),
    hx_late_loss_12_weeks = case_when(
      hx_late_loss_12_weeks == 0 ~ "0",
      hx_late_loss_12_weeks == 1 ~ "1",
      hx_late_loss_12_weeks >= 2 ~ "2+"
    ),
    hx_late_loss_16_weeks = case_when(
      hx_late_loss_16_weeks == 0 ~ "0",
      hx_late_loss_16_weeks == 1 ~ "1",
      hx_late_loss_16_weeks >= 2 ~ "2+"
    ),
    gravida = case_when(
      gravida == 0 ~ "0",
      gravida == 1 ~ "1",
      gravida == 2 ~ "2",
      gravida == 3 ~ "3",
      gravida >= 4 ~ "4+"
    ),
    hx_losses = case_when(
      hx_losses == 0 ~ "0",
      hx_losses == 1 ~ "1",
      hx_losses >= 2 ~ "2+"
    ),
    live_births = case_when(
      live_births == 0 ~ "0",
      live_births == 1 ~ "1",
      live_births >= 2 ~ "2+"
    ),
    hx_preterm_34 = case_when(
      hx_preterm_34 == 0 ~ "0",
      hx_preterm_34 >= 1 ~ "1+"
    ),
    hx_preterm_37 = case_when(
      hx_preterm_37 == 0 ~ "0",
      hx_preterm_37 >= 1 ~ "1+"
    ),
    hx_late_loss_20 = case_when(
      hx_late_loss_20 == 0 ~ "0",
      hx_late_loss_20 >= 1 ~ "1+"
    ),
    thrombophilia = case_when(
      fvl == 1 & pgm == 1 & anti_thrombin == 1 & protein_c == 1 & protein_s == 1 & apla == 1 & other_thrombo_philia == 1 ~ "none",
      fvl == 2 | pgm == 2 ~ "weak",
      protein_c == 2 | protein_s == 2 | other_thrombo_philia == 2 ~ "moderate",
      anti_thrombin == 2 | apla == 2 | fvl == 3 | pgm == 3 | thrombo_count > 1 ~ "strong"
    ),
    thrombophilia_new = case_when(
      thrombophilia == "none" ~ "no",
      thrombophilia == "weak" | thrombophilia == "moderate" | thrombophilia == "strong" ~ "yes"
    ),
    smoking_new = case_when(
      smoking == 1 | smoking == 3 ~ "no",
      smoking == 2 ~ "yes"
    )
  )

data %<>%
  mutate(across(c(
    bmi,
    race,
    paternity,
    hx_late_loss_12_weeks,
    hx_late_loss_16_weeks,
    hx_sga_10th,
    hx_pre_eclampsia,
    hx_abruption_delivery,
    thrombophilia,
    smoking,
    chronic_hyper_tension,
    diabetes,
    vte,
    family_history_of_vte,
    family_history_of_arterial_disea,
    gravida,
    hx_losses,
    hx_late_loss_20,
    live_births,
    hx_sga_5,
    hx_sga_3,
    hx_preterm_34,
    hx_preterm_37,
    hx_severe_pre_eclampsia,
    hx_early_pre_eclampsia,
    hx_abruption,
    aspirin,
    smoking_new,
    thrombophilia_new
  ), as.factor))

# remove 1 patient with missing age
data %<>% filter(!is.na(age__years_))

# get SMDs
# tidy_smd(data, tab1_covs, .group = single, na.rm = T)
# rm(caption, tab, tab_formula, tab_data)



#### 8. RESTRICTIONS PRIOR TO ANALYSES ####

# examine whether some patients are systematically not represented in single-site trials
# as will violate positivity in the weighting model
covs_list <- c(
  "ga_at_start__weeks_",
  "age__years_",
  "bmi",
  "hx_late_loss_12_weeks", 
  "hx_late_loss_16_weeks",
  "hx_sga_10th",
  "hx_pre_eclampsia",
  "hx_abruption_delivery",
  "thrombophilia",
  "smoking",
  "gravida",
  "hx_losses",
  "hx_late_loss_20",
  "live_births",
  "hx_sga_5",
  "hx_sga_3",
  "hx_preterm_34",
  "hx_preterm_37",
  "hx_early_pre_eclampsia",
  "hx_abruption",
  "aspirin"
)

# covariate distribution by trial
tab2_formula <- as.formula(paste("~", paste(covs_list, collapse = "+"), "|unique_id"))
tab2 <- table1(
  tab2_formula,
  data = data,
  overall = c(right = 'Total'),
  caption = caption
)

tab2

# exclusion 1: any pregnancy loss
# exclusion 2: any gravida >2
# exclusion 3: GA at active start >10 weeks
analytic_data <- data

analytic_data %<>%
  mutate(ga_at_start__weeks_ = as.factor(ga_at_start__weeks_))

table(analytic_data$hx_losses != 0)
table(analytic_data$hx_late_loss_12_weeks != 0)
table(analytic_data$hx_late_loss_16_weeks != 0)
table(analytic_data$hx_late_loss_20 != 0)
table(analytic_data$gravida != 2)
table(analytic_data$ga_at_start__weeks_ == "10 to <16 weeks")
table(analytic_data$ga_at_start__weeks_ == "16 to <20 weeks")
table(analytic_data$ga_at_start__weeks_ == ">=20 weeks")
table(analytic_data$unique_id, analytic_data$hx_preterm_34)

to_exclude <- analytic_data %>% 
  filter(
    hx_losses != 0 |
    hx_late_loss_12_weeks != 0 |
    hx_late_loss_16_weeks != 0 |
    hx_late_loss_20 != 0 |
    gravida != 2 |
    ga_at_start__weeks_ == "10 to <16 weeks" |
    ga_at_start__weeks_ == "16 to <20 weeks" |
    ga_at_start__weeks_ == ">=20 weeks"
  )

analytic_data %<>%
  filter(!new_id %in% to_exclude$new_id)


table(analytic_data$unique_id, analytic_data$hx_preterm_34)

analytic_covs_list <- c(
  #"ga_at_start__weeks_", # removed because excluded to keep only <16
  "age__years_",
  "hx_sga_10th",
  "hx_pre_eclampsia",
  #"hx_abruption_delivery", # removed bc extreme weights (low counts -> positivity)
  "thrombophilia_new",
  "smoking_new",
  "hx_sga_5",
  #"hx_sga_3", # removed bc causes extreme weights (low counts -> positivity)
  # "hx_preterm_34", removed bc extreme weights (low counts -> positivity)
  # "hx_preterm_37", # positivity violation by trt in S = 1 and S = 0
  "hx_early_pre_eclampsia"
  #"hx_abruption", # removed bc causes near perfect separation of outcome
  # "aspirin" # positivity violation by trt in S = 1 and S = 0
)

# drop empty factor levels for GA
analytic_data$ga_at_start__weeks_ <- droplevels(analytic_data$ga_at_start__weeks_)

# covariates by trial after applying exclusions
tab3_formula <- as.formula(paste("~", paste(analytic_covs_list, collapse = "+"), "|unique_id"))
tab3 <- table1(
  tab3_formula,
  data = analytic_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab3

#### 4. COMPUTE ATE ESTIMATORS ####

# compute different ATE estimators (adapted from Dahabreh et al.)
# thinking of S as "selection into multi-site trial", so:
# single site patients (single = 1 -> S = 0) -- this is our "target" population
# multi site patients (single = 0 -> S = 1) is our "study" population
# so we are extending inference from multi-site to single-site patients

analytic_data %<>%
  mutate(S = if_else(single == 1, 0, 1))

## A) CRUDE ATE (as a rate difference)

crude_1 <- sum(analytic_data$outcome[analytic_data$trt == 1])/sum(analytic_data$trt == 1)
crude_0 <- sum(analytic_data$outcome[analytic_data$trt == 0])/sum(analytic_data$trt == 0)
crude <- crude_1 - crude_0

round(crude_1, 3)
round(crude_0, 3)
round(crude, 3)

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

round(OM_1, 3)
round(OM_0, 3)
round(OM, 3)

## C) INVERSE ODDS WEIGHTING

selection_model <- as.formula(paste("S ~", paste(analytic_covs_list, collapse = "+")))
treatment_model <- as.formula(paste("trt ~", paste(analytic_covs_list, collapse = "+")))

# generate weights
S1data <- subset(analytic_data, S == 1)
summary(S1data[analytic_covs_list])

w_reg <- glm(selection_model, family = "binomial", data = analytic_data)
summary(w_reg)
ps <- predict(w_reg, newdata = analytic_data, type = "response")
analytic_data$ps <- ps

w_reg2 <- glm(treatment_model, family = "binomial", data = S1data)
summary(w_reg2)
pa <- predict(w_reg2, newdata = analytic_data, type = "response")
analytic_data$pa <- pa

w = (analytic_data$trt * analytic_data$S * (1-ps))/(ps*pa) + ((1-analytic_data$trt) * analytic_data$S*(1-ps)/(ps*(1-pa)))
analytic_data$w <- w

summary(analytic_data$ps)
summary(analytic_data$pa)
summary(analytic_data$w)

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
    mean_w = mean(w)
  ) %>% 
  ungroup()

test <- analytic_data %>% 
  filter(w > 30)

summary(test[analytic_covs_list])

ggplot(analytic_data, aes(x = ps, fill = as.factor(S))) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Selection probability",
    y = "Density",
    fill = "Multi-site trial",
    title = "Selection Probability Density Plot by Trial Type"
  ) +
  theme_minimal()

# unstabilized IOW
# this should be the same regardless of method used...
A <- analytic_data$trt
S <- analytic_data$S
w <- analytic_data$w
Y <- analytic_data$outcome

# original version: divide by N0 (sample size in S = 0, or single-site patients) - A

IOW1a_1 <- (sum((1-S))^-1) * sum(A*S*w*Y)
IOW1a_0 <- (sum((1-S))^-1) * sum((1-A)*S*w*Y)

IOW1a = IOW1a_1 - IOW1a_0
IOW1a

round(IOW1a_1, 3)
round(IOW1a_0, 3)
round(IOW1a, 3)

# updated version: divide by sum of weights - B

IOW1b_1 <- (sum(w)^-1) * sum(A*S*w*Y)
IOW1b_0 <- (sum(w)^-1) * sum((1-A)*S*w*Y)

IOW1b = IOW1b_1 - IOW1b_0
IOW1b

round(IOW1b_1, 3)
round(IOW1b_0, 3)
round(IOW1b, 3)

# stabilized IOW
S0data <- subset(analytic_data, S == 0)
S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
IOW1mod <- glm(formula = outcome ~ 1, data = S1data_A1, weights = w)
p1 <- predict(IOW1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
IOW0mod <- glm(formula = outcome ~ 1, data = S1data_A0, weights = w)
p0 <- predict(IOW0mod, newdata = S0data, type = "response")
IOW2_1 <- mean(p1)
IOW2_0 <- mean(p0)
IOW2 <- mean(p1) - mean(p0)
IOW2

round(IOW2_1, 3)
round(IOW2_0, 3)
round(IOW2, 3)

# stabilize by dividing each weight by mean weight
analytic_data$w_norm <- with(analytic_data, ave(w, trt, FUN = function(x) x / mean(x)))
w_norm <- analytic_data$w_norm

IOW1c_1 <- (sum((1-S))^-1) * sum(A*S*w_norm*Y)
IOW1c_0 <- (sum((1-S))^-1) * sum((1-A)*S*w_norm*Y)

IOW1c = IOW1c_1 - IOW1c_0
IOW1c

round(IOW1c_1, 3)
round(IOW1c_0, 3)
round(IOW1c, 3)

IOW1d_1 <- (sum(w_norm)^-1) * sum(A*S*w_norm*Y)
IOW1d_0 <- (sum(w_norm)^-1) * sum((1-A)*S*w_norm*Y)

IOW1d = IOW1d_1 - IOW1d_0
IOW1d

round(IOW1d_1, 3)
round(IOW1d_0, 3)
round(IOW1d, 3)

# stabilize using weighted mean (same as above)
# S1data_A1 <- subset(analytic_data, S == 1 & trt == 1)
# S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
# 
# IOW2_1 <- weighted.mean(S1data_A1$outcome, S1data_A1$w)
# IOW2_0 <- weighted.mean(S1data_A0$outcome, S1data_A0$w)
# IOW2 <- IOW2_1 - IOW2_0
# IOW2

# stabilize manually (same as above)
# w1 <- S1data_A1$w
# w0 <- S1data_A0$w
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
w <- analytic_data$w

# original version: divide by N0 (sample size in S = 0, or single-site patients) - A

DR1a_1 <- (sum((1-S))^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
DR1a_0 <- (sum((1-S))^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)

DR1a <- DR1a_1 - DR1a_0
DR1a

round(DR1a_1, 3)
round(DR1a_0, 3)
round(DR1a, 3)

# updated version: divide by sum of weights - B

DR1b_1 <- (sum(w)^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
DR1b_0 <- (sum(w)^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)

DR1b <- DR1b_1 - DR1b_0
DR1b

round(DR1b_1, 3)
round(DR1b_0, 3)
round(DR1b, 3)

# DR using stabilized IOW
sum1_DR2 <- sum(S*A*w*(Y-p1))
sum0_DR2 <- sum(S*(1-A)*w*(Y-p0))
norm1 <- (sum(S*A*w))^-1
norm0 <- (sum(S*(1-A)*w))^-1
DR2_1 <- norm1 * sum1_DR2 + (sum(1-S)^-1) * sum((1-S) * p1)
DR2_0 <- norm0 * sum0_DR2 + (sum(1-S)^-1) * sum((1-S) * p0)
DR2 <- DR2_1 - DR2_0
DR2

round(DR2_1, 3)
round(DR2_0, 3)
round(DR2, 3)

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
DR1mod <- glm(formula = outcome_model, data = S1data_A1, weights = w, family = "binomial" (link = "logit"))
p1 <- predict(DR1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(analytic_data, S == 1 & trt == 0)
DR0mod <- glm(formula = outcome_model, data = S1data_A0, weights = w, family = "binomial" (link = "logit"))
p0 <- predict(DR0mod, newdata = S0data, type = "response")
DR3_1 <- mean(p1)
DR3_0 <- mean(p0)
DR3 <- DR3_1 - DR3_0
DR3

round(DR3_1, 3)
round(DR3_0, 3)
round(DR3, 3)

#### 5. CHECK BALANCE AND UPDATE PS MODEL ACCORDINGLY ####

covs <- analytic_data %>%
  select(all_of(analytic_covs_list))

# compare weighted multi-site patients vs unweighted single-site patients
# want SMD < 0.1

# let's define "new_w" which gives 1 to patients with S = 0 (rather than 0 as is currently)
# so using new_w means looking at unweighted single-site (S = 0) + weighted mutli-site (S = 1) together

analytic_data %<>%
  mutate(
    new_w = if_else(S == 0, 1, w)
  )

tidy_smd(analytic_data, analytic_covs_list, .group = S, .wts = c(new_w))

bal_single <- bal.tab(S ~ covs, data = analytic_data,
                      weights = 'new_w',
                      binary = "std", continuous = "std", 
                      stats = c('mean.diffs'),
                      s.d.denom = 'pooled',
                      un = TRUE, 
                      thresholds = c(m = 0.1)
) 

bal_tab_single <- bal_single$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (Wsingle)' = Diff.Adj,
                'Balance' = M.Threshold)

bal_tab_single$Variable <- row.names(bal_tab_single)

# write_xlsx(bal_tab_single , "bal_tab_Wsingle.xlsx")

p <- love.plot(
  bal_single,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Covariate Balance with Inverse Odds Weights to Resemble Single-Center Trials',
  colors = c('#0d0887', '#fdbd22'),
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2
)

p <- p + theme(
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

p

#ggsave("loveplot_w.png", plot = p, width = 6, height = 6, units = "in", bg = 'white')
dev.off()



#### 6. FIT WEIGHTED GEE MODEL FOR OUTCOME ####  

# for glmgee use type="bias-corrected" to use Mancl & DeRouen (2001)
# adapted covariance estimator for small samples
# glmgee(outcome ~ trt, id = site, family = binomial("logit"), weights = w, type = "bias-corrected", corstr = "...")

# first we need to arrange data by id and grouping var (trial)
analytic_data %<>% 
  arrange(new_id, unique_id)

# and we are only interested in (weighted) multi-site
# so remove single-site patients (who all have weight of 0) for now
# because GEE model cannot handle non-positive weight

analytic_data_S1 <- analytic_data %>% 
  filter(S == 1)

## A) CRUDE GEE

# crude from whole population for now

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
round(GEE_crude, 3)
rm(treated, untreated)

## B) IOW GEE

# only in weighted multi-site (S = 1) patients
# because GEE model cannot handle weight of 0 for single-site patients

# first use logit to get separate probability by treatment group
gee_weighted_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "logit"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  weights = analytic_data_S1$w
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

rm(treated, untreated)

# now use identity to get the risk difference directly
# and apply Mancl & deRouen small sample bias variance correction after

gee_weighted_model <- glmgee(
  outcome ~ trt,
  data = analytic_data_S1,
  id = unique_id,
  family = binomial(link = "identity"),
  corstr = "exchangeable", # can try independence, exchangeable, unstructured
  weights = analytic_data_S1$w
)

vcov_bias_corrected <- vcov(gee_weighted_model, type = "bias-corrected")
se_bias_corrected <- sqrt(diag(vcov_bias_corrected))
cbind(Estimate = coef(gee_weighted_model), `Std.error` = se_bias_corrected)

GEE_IOW1 <- coef(gee_weighted_model)[2]

# C) STABILZIED IOW GEE - DEFINITELY SOMETHING WRONG HERE

analytic_data_S1_norm <- analytic_data_S1
analytic_data_S1_norm$w_norm <- with(analytic_data_S1_norm, ave(w, trt, FUN = function(x) x / mean(x)))

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

#### 7. FORMAT RESULTS ####

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
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# transpose results data frame for easier reading
results_transposed <- as.data.frame(t(results[,-1]))
colnames(results_transposed) <- results$measure
results_transposed$Method <- rownames(results_transposed)
results_transposed <- results_transposed[, c(ncol(results_transposed), 1:(ncol(results_transposed)-1))]
rownames(results_transposed) <- NULL

