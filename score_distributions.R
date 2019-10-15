# Script to generate density plots for clinical scores in UKBB for SIGN derivation

# Dependencies
library(data.table)
library(ggplot2)
library(reshape2)

# Load data
load(file="/Volumes/medpop_afib/skhurshid/SIGN/phenos_complete_101119.RData")

# Create separate data.tables for cases/controls
prevalent_af <- phenos_complete[Prev_Atrial_fibrillation_or_flutter==1,]
noaf <- phenos_complete[Prev_Atrial_fibrillation_or_flutter==0,]

### RAW SCORE DENSITY PLOTS
###################################################################################### CHARGE
# Generate distribution
x <- list(v1=prevalent_af$charge,v2=noaf$charge)
data <- melt(x)

# Plot distribution
png(file='density_charge.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(8,16,1),expand=c(0,0.1),limits=c(8,16)) +
  scale_y_continuous(breaks=seq(0,0.60,0.05),expand=c(0,0),limits=c(0,0.60)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF Score',y='density')
dev.off()

############################################################################# CHARGE REWEIGHTED
# Generate distribution
x <- list(v1=prevalent_af$charge_reweighted,v2=noaf$charge_reweighted)
data <- melt(x)

# Plot distribution
png(file='density_charge_rw.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(8,21,1),expand=c(0,0.1),limits=c(8,21)) +
  scale_y_continuous(breaks=seq(0,0.45,0.05),expand=c(0,0),limits=c(0,0.45)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF Score (Reweighted)',y='density')
dev.off()

########################################################################## SIGN
# Generate distribution
x <- list(v1=prevalent_af$sign,v2=noaf$sign)
data <- melt(x)

# Plot distribution
png(file='density_sign.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(3,11,1),expand=c(0,0.1),limits=c(3,11)) +
  scale_y_continuous(breaks=seq(0,0.55,0.05),expand=c(0,0),limits=c(0,0.55)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='SIGN Score',y='density')
dev.off()

### PREDICTED PROBABILITY DENSITY PLOTS
###################################################################################### CHARGE
# Generate distribution
x <- list(v1=prevalent_af$charge_pred,v2=noaf$charge_pred)
data <- melt(x)

# Plot distribution
png(file='density_charge_pred.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(0,6,1),expand=c(0,0.05),limits=c(0,6)) +
  scale_y_continuous(breaks=seq(0,3.5,0.5),expand=c(0,0),limits=c(0,3.5)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF risk (CHARGE-AF, %)',y='density')
dev.off()

############################################################################# CHARGE REWEIGHTED
# Generate distribution
x <- list(v1=prevalent_af$charge_rw_pred,v2=noaf$charge_rw_pred)
data <- melt(x)

# Plot distribution
png(file='density_charge_rw_pred.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(0,20,5),expand=c(0,0.05),limits=c(0,20)) +
  scale_y_continuous(breaks=seq(0,0.9,0.10),expand=c(0,0),limits=c(0,0.9)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF risk (CHARGE-AF Reweighted, %)',y='density')
dev.off()

########################################################################## SIGN
# Generate distribution
x <- list(v1=prevalent_af$sign_pred,v2=noaf$sign_pred)
data <- melt(x)

# Plot distribution
png(file='density_sign_pred.png',height=540,width=760)
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55,bw='nrd0',adjust=1.4) +
  scale_x_continuous(breaks=seq(0,20,5),expand=c(0,0.05),limits=c(0,20)) +
  scale_y_continuous(breaks=seq(0,0.90,0.1),expand=c(0,0),limits=c(0,0.90)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('prevalent AF','no prevalent AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted AF risk (SIGN, %)',y='density')
dev.off()
