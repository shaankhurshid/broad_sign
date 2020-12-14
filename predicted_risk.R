# Script for generating predicted risk in UKBB for SIGN

# Dependencies
library(data.table)
library(plyr)

# Load data
load(file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_complete_121120.RData')

# Function to convert linear predictor of logistic regression model to predicted risk
pred_logistic <- function(score,intercept){
  pred <- exp(intercept+score) / (1 + (exp(intercept+score)))
  return(pred)
}

# Recreate logistic regression models to obtain intercepts
intercept_charge <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ age_5 + race_binary + ht_10 + wt_15 + sbp_20
                                               + dbp_10 + active_smoker + bp_med_combined + prevalent_Diabetes_All
                                               + prevalent_Heart_Failure + prevalent_Myocardial_Infarction,family='binomial',data=phenos_complete)$coefficients[1]

intercept_sign <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ age_5 + factor(sex,levels=c('female','male')) 
                      + race_binary + active_smoker + prevalent_Hypertension + prevalent_Diabetes_All
                      + prevalent_Coronary_Artery_Disease + prevalent_Stroke, family='binomial',data=phenos_complete)$coefficients[1]

# Create predicted risk values
phenos_complete[,charge_rw_pred := lapply(.SD,pred_logistic,intercept=(intercept_charge)),.SDcols='charge_reweighted']
phenos_complete[,sign_pred := lapply(.SD,pred_logistic,intercept=(intercept_sign)),.SDcols='sign']

# Multiply by 100
phenos_complete[,':='(charge_rw_pred = charge_rw_pred*100,sign_pred = sign_pred*100)]