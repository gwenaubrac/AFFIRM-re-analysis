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

# replace preeclampsia by outcome
# replace 99 by NA

#### 1. LOAD PACKAGES ####

library(geepack)
library(PSweight)
library(dplyr)
library(magrittr)
library(table1)
library(writexl)
library(tidysmd)
library(medicaldata)
library(cobalt)
library(ggplot2)



#### 2. READ AND FORMAT DATA ####

setwd("my_path_res")

# let us define treatment (LMWH or not)
# and also create a variable for single-site trials
# given that the two trials in France are single-site trials

data %<>%
  mutate(
    single = if_else(Country == "France", 1, 0),
    trt = if_else(Allocation == "LMWH" | Allocation == "LMWH + Aspirin", 1, 0),
  )

table(data$single)
table(data$trt) # should have 480 treated and 483 untreated

# let us define the outcome, which is a composite of:
# early-onset PE (<34 weeks) or severe PE
# birth of SGA neonate (<5th percentile)
# placental abruption
# late unexplained pregnancy loss (>= 20 weeks)

data$PE_days_elapsed <- as.numeric(data$"Pre-eclamspia onset" - data$"Active trial start")
data$PE_weeks_elapsed <- data$PE_days_elapsed/7
data$GA_at_PE <- data$"GA at start (weeks)" + data$PE_weeks_elapsed

data %<>%
  mutate(early_PE = if_else("Pre-eclampsia" == 2 && GA_at_PE < 32, 1, 0))

data$del_days_elapsed <- as.numeric(data$"Pregnancy outcome date" - data$"Active trial start")
data$del_weeks_elapsed <- data$del_days_elapsed / 7
data$GA_at_delivery <- data$"GA at start (weeks)" + data$del_weeks_elapsed

data %<>%
  mutate(
    pregnancy_loss = if_else(
      "Pregnancy outcome" == 2 | # 2 = unexplained loss
        "Prenancy outcome" == 3 | # 3 = explained loss
        "Pregnancy outcome" == 5, # 5 = medical termination
      1, 
      0
    )
  )

table(data$pregnancy_loss)

data %<>%
  mutate(late_loss = if_else("pregnancy_loss" == 1 && GA_at_delivery >= 20, 1, 0))

table(data$late_loss)

data %<>%
  mutate(
    outcome = if_else(
      early_PE == 1 |
        "Severe pre-eclampsia" == 1 |
        "SGA<5" == 2 |
        "Abruption no delivery" == 2 |
        "Abruption with delivery" == 2 |
        late_loss == 1, 
      1, 
      0
    )
  )

table(data$outcome)

# let us exclude women with a competing event
# namely those with loss prior to 20 weeks or termination unrelated to outcome

data %<>%
  mutate(early_loss = if_else("pregnancy_loss" == 1 && GA_at_delivery < 20, 1, 0))

table(data$early_loss) # should be N=79
table(data$unrelated_termination) # should be N=7

data %<>%
  filter(early_loss == 0 & unrelated_termination == 0)

# let us define treatment as LMWH (whether with or without aspirin)
# and aspirin as an EMM

data %<>%
  mutate(
    aspirin = if_else(Allocation = 2 | Allocation == 3, 1, 0)
  )

table(data$aspirin)



#### 3. DESCRIBE STUDY POPULATION ####

data %<>%
  mutate(
    BMI = "Baseline weight (kg)" / ("Heigh (cm)"^2)
  )

# define variable for GA at LMWH start


cov_list <- c(
  "GA at start (weeks)", # might remove
  "", # added LMWH timing
  "Age (Years)",
  "Race",
  "BMI",
  "Paternity",
  "Hx Late loss 12 weeks",
  "Hx Late loss 16 weeks",
  "Hx SGA 10th",
  "Hx pre-eclampsia",
  "Hx abruption / delivery",
  "FVL",
  "PGM",
  "Antithrombin",
  "Protein C",
  "Protein S",
  "APLA",
  "Other thrombophilia",
  "Smoking",
  "Chronic hypertension",
  "Diabetes",
  "Systolic baseline",
  "Diastolic baseline",
  "VTE",
  "Family history of VTE",
  "Family history of arterial",
  "Gravida",
  "Hx losses",
  "Hx Late loss 20",
  "Live births",
  "Hx SGA 5",
  "Hx SGA 3",
  "Hx preterm 34",
  "Hx preterm 37",
  "Hx severe pre-eclamspsia",
  "Hx early pre-eclampsia",
  "Hx abruption",
  "Aspirin treatment"
  )

# need to wrap variable names in backticks since contain spaces
cov_list_backticked <- paste0("`", cov_list, "`")

summary(data[,cov_list])

tab_data <- data
caption <- 'Baseline Characteristics of Cohort'

# create table 1 by site
tab_formula <- as.formula(paste("~", paste(cov_list_backticked, collapse = "+"), "|single"))
tab <- table1(
  tab_formula,
  data = tab_data,
  overall = c(right = 'Total'),
  caption = caption
)

tab

# write.table (table1 , "tab1.csv", col.names = T, row.names=F, append= F, sep=',')

# get SMDs
tidy_smd(data, cov_list, .group = single)
rm(caption, tab, tab_formula, tab_data)



#### 4. COMPUTE ATE ESTIMATORS ####

# compute difference ATE estimators (adapted from Dahabreh et al.)
# single site (single = 1 -> S = 1) is our "study" population
# multi site (single = 0 -> S = 0) is our "target" population
# so we are extending inference from single-site to multi-site patients
# (single-site made to resemble multi-site)

## A) CRUDE ATE

crude_1 <- sum(data$outcome[data$trt == 1])
crude_0 <- sum(data$outcome[data$trt == 0])
crude <- crude_1 - crude_0


## B) OUTCOME MODELING

# model probability of outcome in treated single-site patients
outcome_model <- as.formula(paste("outcome ~", paste(cov_list_backticked, collapse = "+")))
S1data_A1 <- subset(data, single == 1 & trt == 1)
OM1mod <- glm(formula = outcome_model, data = S1data_A1)

# predict probability of outcome in whole population based on treated
p1 <- predict(OM1mod, newdata = data, type = "response")
data$p1 <- p1

# model probability of outcome in untreated single-site patients
S1data_A0 <- subset(data, single == 1 & trt == 0)
OM0mod <- glm(formula = outcome_model, data = S1data_A0)

# predict probability of outcome in whole population based on untreated
p0 <- predict(OM0mod, newdata = data, type = "response")
data$p0 <- p0

S0sub <- subset(data, single == 0)
OM_1 <- mean(S0sub$p1)
OM_0 <- mean(S0sub$p0)
OM <- mean(S0sub$p1) - mean(S0sub$p0)
OM



## C) INVERSE ODDS WEIGHTING

selection_model <- as.formula(paste("single ~", paste(cov_list_backticked, collapse = "+")))
treatment_model <- as.formula(paste("trt ~", paste(cov_list_backticked, collapse = "+")))

# generate weights
S1data <- subset(data, single == 1)
w_reg <- glm(selection_model, family = "binomial", data = data)
ps <- predict(w_reg, newdata = data, type = "response")

w_reg2 <- glm(treatment_model, family = "binomial", data = S1data)
pa <- predict(w_reg2, newdata = data, type = "response")

w = (data$trt * data$single * (1-ps))/(ps*pa) + ((1-data$trt) * data$single*(1-ps)/(ps*(1-pa)))
data$w <- w

# unstabilized IOW
# this should be the same regardless of method used...
A <- data$trt
S <- data$single
w <- data$w
Y <- data$preeclampsia
IOW1_1 <- (sum((1-S))^-1) * sum(A*S*w*Y)
IOW1_0 <- (sum((1-S))^-1) * sum((1-A)*S*w*Y)
IOW1 = IOW1_1 - IOW1_0
IOW1

# stabilized IOW
S0data <- subset(data, single == 0)
S1data_A1 <- subset(data, single == 1 & trt == 1)
IOW1mod <- glm(formula = preeclampsia ~ 1, data = S1data_A1, weights = w)
p1 <- predict(IOW1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(data, single == 1 & trt == 0)
IOW0mod <- glm(formula = preeclampsia ~ 1, data = S1data_A0, weights = w)
p0 <- predict(IOW0mod, newdata = S0data, type = "response")
IOW2_1 <- mean(p1)
IOW2_0 <- mean(p0)
IOW2 <- mean(p1) - mean(p0)
IOW2

# stabilize using weighted mean
# S0data <- subset(data, single == 0)
# S1data_A1 <- subset(data, single == 1 & trt == 1)
# S1data_A0 <- subset(data, single == 1 & trt == 0)
# 
# IOW2_1 <- weighted.mean(S1data_A1$preeclampsia, S1data_A1$w)
# IOW2_0 <- weighted.mean(S1data_A0$preeclampsia, S1data_A0$w)
# IOW2 <- IOW2_1 - IOW2_0
# IOW2

# stabilize manually
# w1 <- S1data_A1$w
# w0 <- S1data_A0$w
# y1 <- S1data_A1$preeclampsia
# y0 <- S1data_A0$preeclampsia
# 
# IOW2_1 <- sum(w1 * y1) / sum(w1)
# IOW2_0 <- sum(w0 * y0) / sum(w0)
# IOW2 <- IOW2_1 - IOW2_0
# IOW2



## D) DOUBLY ROBUST ESTIMATOR

# DR unstabilized
A <- data$trt
S <- data$single
Y <- data$preeclampsia
p1 <- data$p1
p0 <- data$p0
w <- data$w

DR1_1 <- (sum((1-S))^-1) * sum(S*A*w*(Y-p1) + (1-S) * p1)
DR1_0 <- (sum((1-S))^-1) * sum(S*(1-A)*w*(Y-p0) + (1-S) * p0)
DR1 <- DR1_1 - DR1_0
DR1

# DR using stabilized IOW
sum1_DR2 <- sum(S*A*w*(Y-p1))
sum0_DR2 <- sum(S*(1-A)*w*(Y-p0))
norm1 <- (sum(S*A*w))^-1
norm0 <- (sum(S*(1-A)*w))^-1
DR2_1 <- norm1 * sum1_DR2 + (sum(1-S)^-1) * sum((1-S) * p1)
DR2_0 <- norm0 * sum0_DR2 + (sum(1-S)^-1) * sum((1-S) * p0)
DR2 <- DR2_1 - DR2_0
DR2

# DR using weighted regression
S0data <- subset(data, single == 0)
S1data_A1 <- subset(data, single == 1 & trt == 1)
DR1mod <- glm(formula = outcome_model, data = S1data_A1, weights = w)
p1 <- predict(DR1mod, newdata = S0data, type = "response")
S1data_A0 <- subset(data, single == 1 & trt == 0)
DR0mod <- glm(formula = outcome_model, data = S1data_A0, weights = w)
p0 <- predict(DR0mod, newdata = S0data, type = "response")
DR3_1 <- mean(p1)
DR3_0 <- mean(p0)
DR3 <- DR3_1 - DR3_0
DR3



#### 5. CHECK BALANCE AND UPDATE PS MODEL ACCORDINGLY ####

covs <- data %>%
  select(cov_list)

# compare weighted multi-site patients vs unweighted single-site patients
# want SMD < 0.1
tidy_smd(data, cov_list, .group = single, .wts = c(w))

bal_single <- bal.tab(single ~ covs, data = data,
                      weights = 'w',
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

# first we need to arrange data by id and grouping var (site)
data <- data %>% 
  arrange(id, site)

## A) CRUDE GEE

# GEE coefficients have same interpretation as in logistic model
# i.e., odds ratio
# see: https://library.virginia.edu/data/articles/getting-started-with-generalized-estimating-equations 

gee_crude_model <- geeglm(
  outcome ~ trt,
  id = site,
  data = data,
  family = binomial("logit"),
  corstr = "independence"
)

summary(gee_crude_model)
exp(gee_crude_model$coefficients[2]) # this is the OR

# we got the OR but we want the outcome probability in trt1, trt0, and the ATE
# can follow steps from here:
# https://www.r-bloggers.com/2021/06/estimating-a-risk-difference-and-confidence-intervals-using-logistic-regression/

treated <- data
treated$trt <- 1

untreated <- data
untreated$trt <- 0

treated$pred <- predict(gee_crude_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_crude_model, newdata = untreated, type = "response")

GEE_crude_1 <- mean(treated$pred)
GEE_crude_0 <- mean(untreated$pred)

GEE_crude <- GEE_crude_1 - GEE_crude_0
rm(treated, untreated)


## B) IOW GEE
gee_weighted_model <- geeglm(
  outcome ~ trt,
  data = data,
  id = site,
  family = binomial("logit"),
  corstr = "independence", # can try independence, exchangeable, unstructured
  weights = data$w
)

summary(gee_weighted_model)
exp(gee_weighted_model$coefficients[2]) # this is the OR

treated <- data
treated$trt <- 1

untreated <- data
untreated$trt <- 0

treated$pred <- predict(gee_weighted_model, newdata = treated, type = "response")
untreated$pred <- predict(gee_weighted_model, newdata = untreated, type = "response")

GEE_IOW_1 <- mean(treated$pred)
GEE_IOW_0 <- mean(untreated$pred)

GEE_IOW <- GEE_IOW_1 - GEE_IOW_0
rm(treated, untreated)



# C) STABILZIED IOW GEE

S1data_A1 <- subset(data, single == 1 & trt == 1)
S1data_A0 <- subset(data, single == 1 & trt == 0)

GEE_IOW2_1 <- weighted.mean(S1data_A1$outcome, S1data_A1$w)
GEE_IOW2_0 <- weighted.mean(S1data_A0$outcome, S1data_A0$w)
GEE_IOW2 <- GEE_IOW2_1 - GEE_IOW2_0
GEE_IOW2

#### 7. FORMAT RESULTS ####

results <- data.frame(
  measure = c('Proportion Events (A=1)', 'Proportion Events (A=0)', 'ATE'),
  "Crude" = c(crude_1, crude_0, crude),
  "Outcome Model" = c(OM_1, OM_0, OM),
  "IOSW" = c(IOW1_1, IOW1_0, IOW1),
  "Stabilized IOSW" = c(IOW2_1, IOW2_0, IOW2),
  "Doubly Robust" = c(DR1_1, DR1_0, DR1),
  "Doubly Robust (stabilized)" = c(DR2_1, DR2_0, DR2),
  "Doubly Robust (weighted regression)" = c(DR3_1, DR3_0, DR3),
  "Crude (GEE)" = c(GEE_crude_1, GEE_crude_0, GEE_crude),
  "IOSW (GEE)"= c(GEE_IOW_1, GEE_IOW_0, GEE_IOW),
  "Stabilized IOSW (GEE)"= c(GEE_IOW2_1, GEE_IOW2_0, GEE_IOW2),
)
