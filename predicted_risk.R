# Script for generating predicted risk in UKBB for SIGN

# Dependencies
library(data.table)
library(plyr)

# Load data
load(file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_complete_101119.RData')

# Function to convert linear predictor of logistic regression model to predicted risk
pred_logistic <- function(score,intercept){
  pred <- exp(intercept+score) / (1 + (exp(intercept+score)))
  return(pred)
}

# Recreate logistic regression models to obtain intercepts
intercept_charge <- glm(Prev_Atrial_fibrillation_or_flutter ~ age_5 + race_binary + ht_10 + wt_15 + sbp_20
                                               + dbp_10 + active_smoker + bp_med_combined + Prev_Diabetes_All
                                               + Prev_Heart_Failure_V2 + Prev_Myocardial_Infarction,family='binomial',data=phenos_complete)$coefficients[1]
intercept_sign <- glm(Prev_Atrial_fibrillation_or_flutter ~ age_5 + factor(sex,levels=c('female','male')) 
                                          + race_binary + active_smoker + Prev_Hypertension + Prev_Diabetes_All
                                          + Prev_Coronary_Artery_Disease + Prev_Stroke,family='binomial',data=phenos_complete)$coefficients[1]

# Create predicted risk values
names <- c('charge','charge_rw')
phenos_complete[,paste0(names,'_pred') := lapply(.SD,pred_logistic,intercept=(intercept_charge)),.SDcols=c('charge','charge_reweighted')]
phenos_complete[,sign_pred := lapply(.SD,pred_logistic,intercept=(intercept_sign)),.SDcols='sign']

# Multiply by 100
phenos_complete[,':='(charge_pred = charge_pred*100,charge_rw_pred = charge_rw_pred*100,sign_pred = sign_pred*100)]