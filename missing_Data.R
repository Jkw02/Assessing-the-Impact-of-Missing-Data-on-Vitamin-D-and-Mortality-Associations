
# Load packages
library(tidyverse)
library(Amelia)
library(stats)
library(dplyr)
library(survival)
library(survminer)
library(broom)

# Load data
vitd_death <- read_dta("vitd_death.dta")

# Check missing vitamin D data 
vitd_death$miss_vitd <- is.na(vitd_death$vitd)
table(vitd_death$miss_vitd)

# Recode BMI
vitd_death$bmi_bcat <- factor(vitd_death$bmi_category, levels = c(2,0,1,3,4,5,6))

# Logistic regression for missing vitamin D
missmodel <- glm(miss_vitd ~ sex + age + factor(bmi_bcat) + diabetes + factor(skin) + factor(alcohol),
                 data = vitd_death, family = binomial())
tidy(missmodel, exponentiate = TRUE, conf.int = TRUE)

# Predictors of vitamin D levels
vitdmodel <- lm(vitd ~ sex + age + factor(bmi_bcat) + diabetes + factor(skin) + factor(alcohol),
                data = vitd_death)
tidy(vitdmodel, conf.int = TRUE)

# Check missing alcohol intake data
vitd_death$alcohol_recode <- vitd_death$alcohol
vitd_death$alcohol_recode[vitd_death$alcohol %in% c(3,7)] <- NA
vitd_death$miss_alc <- ifelse(vitd_death$alcohol %in% c(3,7), 1, 0)
table(vitd_death$miss_alc)

# Logistic regression for missing alcohol data
missmodel2 <- glm(miss_alc ~ sex + age + factor(bmi_bcat) + diabetes + factor(skin),
                  data = vitd_death, family = binomial())
tidy(missmodel2, exponentiate = TRUE, conf.int = TRUE)

# Indicator for current alcohol intake
vitd_death$current <- NA
vitd_death$current[vitd_death$alcohol %in% c(4,5,6)] <- 1
vitd_death$current[vitd_death$alcohol %in% c(1,2)] <- 0
table(vitd_death$current, useNA = "always")

missmodel3 <- glm(current ~ sex + age + factor(bmi_bcat) + factor(skin) + diabetes,
                  data = vitd_death, family = binomial())
tidy(missmodel3, exponentiate = TRUE, conf.int = TRUE)

#Create time variable for survival analysis
vitd_death$end_date <- as.Date(vitd_death$end_date, "1960-01-01")
vitd_death$time_var <- as.numeric(difftime(vitd_death$end_date, vitd_death$date_birth, units = "days")) / 365.25

# Complete case analysis
modelcoxph <- coxph(Surv(time = age, time2 = time_var, event = mortality) ~
                      vitd + sex + factor(bmi_bcat) + factor(alcohol),
                    data = vitd_death, ties = "breslow")
summary(modelcoxph)

# Missing indicator method
vitd_death$vitd_ind_cts <- vitd_death$vitd
vitd_death$vitd_ind_cts[is.na(vitd_death$vitd)] <- 0

modelcoxph2 <- coxph(Surv(time = age, time2 = time_var, event = mortality) ~
                       vitd_ind_cts + miss_vitd + sex + factor(bmi_bcat) + factor(alcohol),
                     data = vitd_death, ties = "breslow")
summary(modelcoxph2)

# Mean imputation 
vitd_death$vitd_meanimp <- vitd_death$vitd
vitd_death$vitd_meanimp[is.na(vitd_death$vitd_meanimp)] <- mean(vitd_death$vitd_meanimp, na.rm = TRUE)

modelcoxph3 <- coxph(Surv(time = age, time2 = time_var, event = mortality) ~
                       vitd_meanimp + sex + factor(bmi_bcat) + factor(alcohol),
                     data = vitd_death, ties = "breslow")
summary(modelcoxph3)

# Regression imputatioN
vitdmodel6 <- lm(vitd ~ sex + age + factor(bmi_bcat) + diabetes + factor(skin) + factor(alcohol),
                 data = vitd_death)
vitd_death$vitd_hat <- predict(vitdmodel6, newdata = vitd_death)
vitd_death$vitd_regressimp <- ifelse(is.na(vitd_death$vitd), vitd_death$vitd_hat, vitd_death$vitd)

modelcoxph4 <- coxph(Surv(time = age, time2 = time_var, event = mortality) ~
                       vitd_regressimp + sex + factor(bmi_bcat) + factor(alcohol),
                     data = vitd_death, ties = "breslow")
summary(modelcoxph4)

#Regression imputation including mortality
vitdmodel7 <- lm(vitd ~ sex + age + factor(bmi_bcat) + diabetes + factor(skin) + factor(alcohol) + mortality,
                 data = vitd_death)
vitd_death$vitd_hat <- predict(vitdmodel7, newdata = vitd_death)
vitd_death$vitd_regressimp <- ifelse(is.na(vitd_death$vitd), vitd_death$vitd_hat, vitd_death$vitd)

modelcoxph5 <- coxph(Surv(time = age, time2 = time_var, event = mortality) ~
                       vitd_regressimp + sex + factor(bmi_bcat) + factor(alcohol),
                     data = vitd_death, ties = "breslow")
summary(modelcoxph5)
