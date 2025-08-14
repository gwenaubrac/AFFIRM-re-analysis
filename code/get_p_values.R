## GET P-VALUES FOR IRDS BETWEEN SINGLE AND WEIGHTED MULTI ##

setwd("Z:/RPLATT/GWEN-OHRI/main_results")
data <- readRDS("../data_clean.R")

data %<>% 
  mutate(S = if_else(single == 1, 0, 1))

res <- read.csv("bs_ci1.csv")


# 1. Wald p-value from CI width

# S1crude_prior, S0crude_prior
# S1crude, S0crude
# S1_OM, S0_OM
# S1_w, S0_w
# S1_DR1, S0_DR1
# S1_GEE_IOW1, S0_GEE_IOW1

S1name <- "S1_DR1"
S0name <- "S0_w"

S1 <- res$estimate[res$var_name == S1name]
S1_l <- res$lower_ci[res$var_name == S1name]
S1_u <- res$upper_ci[res$var_name == S1name]

S0 <- res$estimate[res$var_name == S0name]
S0_l <- res$lower_ci[res$var_name == S0name]
S0_u <- res$upper_ci[res$var_name == S0name]

# get standard error from CIs
S1_se <- (S1_u - S1_l) / (2*1.96)
S0_se <- (S0_u - S0_l) / (2*1.96)

diff <- S1 - S0
diff_se <- sqrt(S1_se^2 + S0_se^2)

z <- diff / diff_se
wald_p <- 2*pnorm(-abs(z))

diff
diff_se
z
wald_p

# 2. Calculate Cochran's Q and I2

# get estimates and std error for both "studies" where S1 and S0 are the "studies"
theta <- c(S1, S0)
se <- c(S1_se, S0_se)

w <- 1/(se^2)

theta_pooled <- sum(w*theta)/sum(w)
Q <- sum(w*(theta-theta_pooled)^2)

k <- length(theta)
df <- k-1
p_value <- pchisq(Q, df, lower.tail = FALSE)

I2 <- max(0, (Q-df)/Q) * 100

Q
p_value
I2

# cross check with output from R package
library(metafor)
meta_res <- rma.uni(yi = theta, sei = se, method = "FE") 
summary(meta_res)

# 3. Fisher p-value from raw proportions and N
# events in S1 = (prop events in S1 treated + prop events in S1 untreated) * sample size S1

S1_name_trt <- "S1crude_prior_1"
S1_name_untrt <- "S1crude_prior_0"
S1_N <- sum(data$S==1)

S1_prop_events <- (res$estimate[res$var_name == S1_name_trt]/100 + res$estimate[res$var_name == S1_name_untrt]/100)
S1_event <- S1_prop_events * S1_N
S1_no_event <- (1-S1_prop_events) * S1_N
S1_event
S1_no_event

S0_name_trt <- "S0crude_prior_1"
S0_name_untrt <- "S0crude_prior_0"
S0_N <- sum(data$S==0)

S0_prop_events <-(res$estimate[res$var_name == S0_name_trt]/100 + res$estimate[res$var_name == S0_name_untrt]/100) 
S0_event <- S0_prop_events * S0_N
S0_no_event <- (1-S0_prop_events) * S0_N
S0_event
S0_no_event

# calculate fisher's exact test
tbl <- matrix(c(S1_event, S1_no_event, S0_event, S0_no_event), nrow = 2, byrow = TRUE)
fisher_test <- fisher.test(tbl)
fisher_p <- fisher_test$p.value
fisher_p

# apply restrictions
to_exclude <- data %>%
  filter(
    hx_any_loss != 0 |
      gravida_cat != 2 |
      live_births_cat != 1 |
      any_vte == 1 |
      GA_at_start_cat == "16 to <20 weeks" |
      GA_at_start_cat == ">=20 weeks"
  )

data %<>%
  filter(!new_id %in% to_exclude$new_id)
