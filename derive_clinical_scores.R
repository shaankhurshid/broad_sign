# Script to develop and compare AF risk scores (CHARGE-AF versus SIGN)

# Dependencies
library(data.table)
library(plyr)
library(rms)

# Load pre-existing phenotype file from PREDICT-AF project
load(file="/Volumes/medpop_afib/skhurshid/PRS/phenos_082819.RData")

# Variable cleanup
## Fix BPs
phenos[,':='(sbp_20 = sbp0/20,dbp_10 = dbp0/10)]

## Binary smoking variable
phenos[,active_smoker := ifelse(tobacco_phenotype=="Current",1,0)]

###################################################################### DEVELOP CLINICAL SCORES

## Completeness variables
phenos[,charge_complete := ifelse(c(is.na(afagevisit0) | is.na(ht) | is.na (wt_all) |
                                      is.na(sbp_20) | is.na(dbp_10) |
                                      is.na(race_binary) | is.na(active_smoker)
                                    | is.na(bp_med_combined) | is.na(Prev_Diabetes_All) |
                                      is.na(Prev_Heart_Failure_V2) | is.na(Prev_Myocardial_Infarction)),0,1)]

phenos[,sign_complete := ifelse(c(is.na(age_5) | is.na(sex) | is.na (race_binary) |
                                    is.na(active_smoker) | is.na(Prev_Hypertension) |
                                    is.na(Prev_Diabetes_All) | is.na(Prev_Coronary_Artery_Disease) |
                                    is.na(Prev_Stroke)),0,1)]

# Restrict set for complete case analysis
phenos_complete <- phenos[c(sign_complete==1 & charge_complete==1)]

#################################### REWEIGHTED CHARGE (for prevalent AF)
## Cleanup

## Make original CHARGE score (using average AF risk as intercept)
phenos_complete[,charge := (avg_risk) + age_5*0.508 + (race_binary=="White")*0.465 + ht_10*0.248 + wt_15*0.115 + sbp_20*0.197
       + dbp_10*(-0.101) + active_smoker*0.359 + bp_med_combined*0.349 + Prev_Diabetes_All*0.237
       + Prev_Heart_Failure_V2*0.701 + Prev_Myocardial_Infarction*0.496]

## Model prevalent AF
prevalent_charge <- glm(Prev_Atrial_fibrillation_or_flutter ~ age_5 + race_binary + ht_10 + wt_15 + sbp_20
                                                                     + dbp_10 + active_smoker + bp_med_combined + Prev_Diabetes_All
                                                                     + Prev_Heart_Failure_V2 + Prev_Myocardial_Infarction,family='binomial',data=phenos_complete)

## Make reweighted CHARGE score
phenos_complete[,charge_reweighted := age_5*prevalent_charge$coefficients[2] + (race_binary=="White")*prevalent_charge$coefficients[3] 
       + ht_10*prevalent_charge$coefficients[4] + wt_15*prevalent_charge$coefficients[5] + sbp_20*prevalent_charge$coefficients[6]
       + dbp_10*prevalent_charge$coefficients[7] + active_smoker*prevalent_charge$coefficients[8] 
       + bp_med_combined*prevalent_charge$coefficients[9] + Prev_Diabetes_All*prevalent_charge$coefficients[10]
       + Prev_Heart_Failure_V2*prevalent_charge$coefficients[11] + Prev_Myocardial_Infarction*prevalent_charge$coefficients[12]]

#################################### SIGN SCORE (for prevalent AF)
## Model prevalent AF
prevalent_sign <- glm(Prev_Atrial_fibrillation_or_flutter ~ age_5 + factor(sex,levels=c('female','male')) 
                        + race_binary + active_smoker + Prev_Hypertension + Prev_Diabetes_All
                        + Prev_Coronary_Artery_Disease + Prev_Stroke, family='binomial',data=phenos_complete)

## Make SIGN score
phenos_complete[sign_complete==1,sign := age_5*prevalent_sign$coefficients[2] + (sex=='male')*prevalent_sign$coefficients[3] 
       + (race_binary=="White")*prevalent_sign$coefficients[4] + active_smoker*prevalent_sign$coefficients[5] 
       + Prev_Hypertension*prevalent_sign$coefficients[6] + Prev_Diabetes_All*prevalent_sign$coefficients[7]
       + Prev_Coronary_Artery_Disease*prevalent_sign$coefficients[8] + Prev_Stroke*prevalent_sign$coefficients[9]]

#################################### Standardized scores
phenos_complete[,':='(charge_std = (charge - mean(charge))/sd(charge),
                      charge_reweighted_std = (charge_reweighted - mean(charge_reweighted))/sd(charge_reweighted),
                      sign_std = (sign - mean(sign))/sd(sign))]

#################################### Compare score performance
# Bootstrap function (takes DT as input, and returns DT)
boot <- function(outcome,classifier,data,runs,size){
  setDF(data)
  out <- cbind(rep(NA,times=runs),rep(NA,times=runs))
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    out[i,1] <- lrm(sample[,outcome] ~ sample[,classifier],data=sample)$stats['R2']
    out[i,2] <- concordance(sample[,outcome] ~ sample[,classifier],data=sample)$concordance
    print(paste0('run ',i,' complete'))
  }
  setDT(data)
  return(out)
}


## Obtain concordances and R2 by bootstrapping
# CHARGE
charge_cstatr2 <- boot(outcome='Prev_Atrial_fibrillation_or_flutter',classifier='charge',data=phenos_complete,runs=200,size=nrow(phenos_complete))
charge_r2 <- c(mean(charge_cstatr2[,1]),mean(charge_cstatr2[,1])-1.96*sd(charge_cstatr2[,1]),mean(charge_cstatr2[,1])+1.96*sd(charge_cstatr2[,1]))
charge_cstat <- c(mean(charge_cstatr2[,2]),mean(charge_cstatr2[,2])-1.96*sd(charge_cstatr2[,2]),mean(charge_cstatr2[,2])+1.96*sd(charge_cstatr2[,2]))

# Reweighted CHARGE
charge_rw_cstatr2 <- boot(outcome='Prev_Atrial_fibrillation_or_flutter',classifier='charge_reweighted',data=phenos_complete,runs=200,size=nrow(phenos_complete))
charge_rw_r2 <- c(mean(charge_rw_cstatr2[,1]),mean(charge_rw_cstatr2[,1])-1.96*sd(charge_rw_cstatr2[,1]),mean(charge_rw_cstatr2[,1])+1.96*sd(charge_rw_cstatr2[,1]))
charge_rw_cstat <- c(mean(charge_rw_cstatr2[,2]),mean(charge_rw_cstatr2[,2])-1.96*sd(charge_rw_cstatr2[,2]),mean(charge_rw_cstatr2[,2])+1.96*sd(charge_rw_cstatr2[,2]))

# SIGN
sign_cstatr2 <- boot(outcome='Prev_Atrial_fibrillation_or_flutter',classifier='sign',data=phenos_complete,runs=200,size=nrow(phenos_complete))
sign_r2 <- c(mean(sign_cstatr2[,1]),mean(sign_cstatr2[,1])-1.96*sd(sign_cstatr2[,1]),mean(sign_cstatr2[,1])+1.96*sd(sign_cstatr2[,1]))
sign_cstat <- c(mean(sign_cstatr2[,2]),mean(sign_cstatr2[,2])-1.96*sd(sign_cstatr2[,2]),mean(sign_cstatr2[,2])+1.96*sd(sign_cstatr2[,2]))

# Obtain model fit and associations using standardized scores
## Model
mod_af_charge <- glm(Prev_Atrial_fibrillation_or_flutter ~ charge_std,family='binomial',data=phenos_complete) 
mod_af_charge_rw <- glm(Prev_Atrial_fibrillation_or_flutter ~ charge_reweighted_std,family='binomial',data=phenos_complete) 
mod_af_sign <- glm(Prev_Atrial_fibrillation_or_flutter ~ sign_std,family='binomial',data=phenos_complete) 
## Params
hr_charge <- with(mod_af_charge,c(exp(coefficients[2]),exp(confint(mod_af_charge,'charge_std')[1]),exp(confint(mod_af_charge,'charge_std')[2])))
hr_charge_rw <- with(mod_af_charge_rw,c(exp(coefficients[2]),exp(confint(mod_af_charge_rw,'charge_reweighted_std')[1]),exp(confint(mod_af_charge_rw,'charge_reweighted_std')[2])))
hr_sign <- with(mod_af_sign,c(exp(coefficients[2]),exp(confint(mod_af_sign,'sign_std')[1]),exp(confint(mod_af_sign,'sign_std')[2])))

# Save output
save(phenos,file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_101119.RData')
save(phenos_complete,file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_complete_101519.RData')
