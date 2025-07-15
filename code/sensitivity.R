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
## This code just defines step before and after running bootstrap function
## (to use proper sample and save results with corresponding name).
## I didn't copy paste the bootstrap function here directly, please refer
## to code provided separately.
##
## Sensitivity analyses:
## Repeat with switch S
## 1. subgroup with aspirin
## 2. subgroup without aspirin
## 3. restrict to GA at start <10 weeks
##
## ---------------------------

# set working directory to analysis, among:
# main_results
# target_multi_site
# aspirin
# no_aspirin
# restrict_GA_10w

setwd("Z:/RPLATT/GWEN-OHRI/aspirin")
data <- readRDS("../data_clean.R")

# 1/2. ASPIRIN

data %<>%
  filter(aspirin == 1)

data %<>%
  mutate(S = if_else(single == 1, 0, 1))

covs_list <- c(
  "age__years_",
  "any_pe",
  "any_sga",
  "bmi_above_30"
)

selection_model <- as.formula("S ~ 
                              ns(age__years_, df = 3) + 
                              any_pe + 
                              any_sga + 
                              bmi_above_30 
                              ")

treatment_model <- as.formula("trt ~ 
                              age__years_ + 
                              any_pe + 
                              any_sga + 
                              bmi_above_30
                              ")

outcome_model <- as.formula(paste("outcome ~", paste(covs_list, collapse = "+")))

setDT(data)

# run bootstrap function + get Ci
write.table (bs_ci_boot , "./results/bs_ci_aspirin1.csv", col.names = T, row.names=F, append= F, sep=',')
write.table (bs_ci_boot , "./results/bs_ci_aspirin0.csv", col.names = T, row.names=F, append= F, sep=',')

# 3. RESTRICT GA TO <10 WEEKS
# uncomment part in bootstrap code that says "UNCOMMENT FOR SENSITIVITY" 
# to add additional restriction
write.table (bs_ci_boot , "./results/bs_ci_restrict10GA.csv", col.names = T, row.names=F, append= F, sep=',')
