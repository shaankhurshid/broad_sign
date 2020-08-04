# Script to develop and compare AF risk scores (CHARGE-AF versus SIGN)

# Dependencies
library(data.table)
library(plyr)
library(rms)

# Load pre-existing phenotype file from PREDICT-AF project
load(file="/Volumes/medpop_afib/skhurshid/accel/phenos_060520.RData")

# Variable cleanup
## Fix BPs
phenos[,':='(sbp_20 = SBP1/20,dbp_10 = DBP1/10)]

## Binary smoking variable
tob0 <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/tobacco_instance_0.csv')
setkey(tob0,sample_id); setkey(phenos,ID)
phenos[tob0,active_smoker := ifelse(is.na(i.value) | i.value==-3,NA,
                                    ifelse(i.value==2,1,0))]

###################################################################### DEVELOP CLINICAL SCORES
## Completeness variables
phenos[,charge_complete := ifelse(!is.na(charge),1,0)]

phenos[,sign_complete := ifelse(c(is.na(age_5) | is.na(sex) | is.na(race_binary) |
                                    is.na(active_smoker) | is.na(prevalent_Hypertension) |
                                    is.na(prevalent_Diabetes_All) | is.na(prevalent_Coronary_Artery_Disease) |
                                    is.na(prevalent_Stroke)),0,1)]

# Restrict set for complete case analysis
phenos_complete <- phenos[sign_complete==1 & charge_complete==1]

## Remove exclusions
withdrawals <- fread('/Volumes/medpop_afib/skhurshid/PRS/withdrawals_020420.csv') # UKBB withdrawals
phenos_complete <- phenos_complete[!(ID %in% withdrawals$V1)]

#################################### SIGN SCORE (for prevalent AF)
## Model prevalent AF
prevalent_sign <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ age_5 + factor(sex,levels=c('female','male')) 
                        + race_binary + active_smoker + prevalent_Hypertension + prevalent_Diabetes_All
                        + prevalent_Coronary_Artery_Disease + prevalent_Stroke, family='binomial',data=phenos_complete)

## Make SIGN score
phenos_complete[sign_complete==1,sign := age_5*prevalent_sign$coefficients[2] + (sex=='male')*prevalent_sign$coefficients[3] 
       + (race_binary=="White")*prevalent_sign$coefficients[4] + active_smoker*prevalent_sign$coefficients[5] 
       + prevalent_Hypertension*prevalent_sign$coefficients[6] + prevalent_Diabetes_All*prevalent_sign$coefficients[7]
       + prevalent_Coronary_Artery_Disease*prevalent_sign$coefficients[8] + prevalent_Stroke*prevalent_sign$coefficients[9]]

#################################### CHARGE-AF (for prevalent AF)
## Model prevalent AF
prevalent_charge <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ age_5 + race_binary + ht_10 + wt_15 + sbp_20
                        + dbp_10 + active_smoker + bp_med_combined + Prevalent_Diabetes_All
                        + prevalent_Heart_Failure + prevalent_Myocardial_Infarction,family='binomial',data=phenos_complete)

## Make reweighted CHARGE score
phenos_complete[,charge_reweighted := age_5*prevalent_charge$coefficients[2] + (race_binary=="White")*prevalent_charge$coefficients[3] 
                + ht_10*prevalent_charge$coefficients[4] + wt_15*prevalent_charge$coefficients[5] + sbp_20*prevalent_charge$coefficients[6]
                + dbp_10*prevalent_charge$coefficients[7] + active_smoker*prevalent_charge$coefficients[8] 
                + bp_med_combined*prevalent_charge$coefficients[9] + prevalent_Diabetes_All*prevalent_charge$coefficients[10]
                + prevalent_Heart_Failure*prevalent_charge$coefficients[11] + prevalent_Myocardial_Infarction*prevalent_charge$coefficients[12]]

#################################### Standardized scores
phenos_complete[,':='(charge_rw_std = (charge_reweighted - mean(charge_reweighted))/sd(charge_reweighted),
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
charge_cstatr2 <- boot(outcome='prevalent_Atrial_fibrillation_or_flutter_v2',classifier='charge_reweighted',data=phenos_complete,runs=200,size=nrow(phenos_complete))
charge_cstat <- c(as.numeric(concordance(prevalent_charge)[1]),as.numeric(concordance(prevalent_charge)[1])-1.96*sd(charge_cstatr2[,1]),as.numeric(concordance(prevalent_charge)[1])+1.96*sd(charge_cstatr2[,1]))
charge_r2 <- c(lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ charge,data=phenos_complete)$stats['R2'],
                  lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ charge,data=phenos_complete)$stats['R2']-1.96*sd(charge_cstatr2[,2]),
                  lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ charge,data=phenos_complete)$stats['R2']+1.96*sd(charge_cstatr2[,2]))

# SIGN
sign_cstatr2 <- boot(outcome='prevalent_Atrial_fibrillation_or_flutter_v2',classifier='sign',data=phenos_complete,runs=200,size=nrow(phenos_complete))
sign_cstat <- c(as.numeric(concordance(prevalent_sign)[1]),as.numeric(concordance(prevalent_sign)[1])-1.96*sd(sign_cstatr2[,1]),as.numeric(concordance(prevalent_sign)[1])+1.96*sd(sign_cstatr2[,1]))
sign_r2 <- c(lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ sign,data=phenos_complete)$stats['R2'],
                  lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ sign,data=phenos_complete)$stats['R2']-1.96*sd(sign_cstatr2[,2]),
                  lrm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ sign,data=phenos_complete)$stats['R2']+1.96*sd(sign_cstatr2[,2]))

# Obtain model fit and associations using standardized scores
## Model
mod_af_charge <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ charge_rw_std,family='binomial',data=phenos_complete) 
mod_af_sign <- glm(prevalent_Atrial_fibrillation_or_flutter_v2 ~ sign_std,family='binomial',data=phenos_complete) 
## Params
hr_charge <- with(mod_af_charge,c(exp(coefficients[2]),exp(confint(mod_af_charge,'charge_rw_std')[1]),exp(confint(mod_af_charge,'charge_rw_std')[2])))
hr_sign <- with(mod_af_sign,c(exp(coefficients[2]),exp(confint(mod_af_sign,'sign_std')[1]),exp(confint(mod_af_sign,'sign_std')[2])))
AIC(prevalent_charge); AIC(prevalent_sign)

# Save output
save(phenos,file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_080320.RData')
save(phenos_complete,file='/Volumes/medpop_afib/skhurshid/SIGN/phenos_complete_080320.RData')
