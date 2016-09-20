#Linear Models for determining predictions of performance on task with self reports 
library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(reshape)
library(multcomp)

alldat<-read.table('/home/rajchahal/552014/sasdata_allcombined.csv',sep=",",header=T)

#escore<-merge(AccScore, alldat,  all=T, by="Subject")
alldat<-within(alldat,{Subject<-factor(Subject)})
#Using RT and Err on incongruent trials flank control as DV's

library(MASS)
#for poor accuracy outliers could use robust regression (e.g., rlm in MASS) 

#RT incongruent Trials correlations: I_MeanRt and I_Error
#Neg Temperament
gm1 <- lm(T_NegativeTemperament ~ I_MeanRT, alldat)
gm1<-lm(T_NegativeTemperament~ I_Error, alldat)


#and Impulsivity
gm1 <- lm(T_Impulsivity ~ I_MeanRT, alldat)
gm1<-lm(T_Impulsivity~ I_Error, alldat)

#Propriety
gm1 <- lm(T_Propriety ~ I_MeanRT, alldat)
gm1<-lm(T_Propriety~ I_Error, alldat)

#Detachment
gm1 <- lm(T_Detachment ~ I_MeanRT, alldat)
gm1<-lm(T_Detachment~ I_Error, alldat)
##**Sig at .00744
cora<-cor(alldat$T_Detachment, alldat$I_Error, use="complete.obs")

#Entitlement
gm1 <- lm(T_Entitlement ~ I_MeanRT, alldat)
gm1<-lm(T_Entitlement~ I_Error, alldat)


#Self Harm
gm1 <- lm(T_SelfHarm ~ I_MeanRT, alldat)
gm1<-lm(T_SelfHarm~ I_Error, alldat)
#*Barely sig at .0765
cora<-cor(alldat$T_SelfHarm, alldat$I_Error, use="complete.obs")


#Aggression
gm1 <- lm(T_Aggression ~ I_MeanRT, alldat)
gm1<-lm(T_Aggression~ I_Error, alldat)
#**Very sig at .0003
cora<-cor(alldat$T_Aggression, alldat$I_Error, use="complete.obs")

#Dependency
gm1 <- lm(T_Dependency ~ I_MeanRT, alldat)
gm1<-lm(T_Dependency~ I_Error, alldat)
#*Barely sig at .054
cora<-cor(alldat$T_Dependency, alldat$I_Error, use="complete.obs")

#Mistrust
gm1 <- lm(T_Mistrust ~ I_MeanRT, alldat)
gm1<-lm(T_Mistrust~ I_Error, alldat)
#** Very Significant at .000259
cora<-cor(alldat$T_Mistrust, alldat$I_Error, use="complete.obs")

#Manipulativeness
gm1 <- lm(T_Manipulativeness ~ I_MeanRT, alldat)
gm1<-lm(T_Manipulativeness~ I_Error, alldat)
#*Very sig at .000417
cora<-cor(alldat$T_Manipulativeness, alldat$I_Error, use="complete.obs")

# Low Self Esteem
gm1 <- lm(T_LowSelfEsteem ~ I_MeanRT, alldat)
gm1<-lm(T_LowSelfEsteem~ I_Error, alldat)

#MPQ neg emotionality
gm1 <- lm(MPS_NET ~ I_MeanRT, alldat)
gm1<-lm(MPS_NET~ I_Error, alldat)
#*Barely sig at .0122
cora<-cor(alldat$MPS_NET, alldat$I_Error, use="complete.obs")

#MPQ alienation
gm1 <- lm(MPS_alT ~ I_MeanRT, alldat)
gm1<-lm(MPS_alT~ I_Error, alldat)
#*Barely sig at .08
cora<-cor(alldat$MPS_alT, alldat$I_Error, use="complete.obs")

#MPQ aggression
gm1<-lm(MPS_agT~I_MeanRT, alldat)
gm1<-lm(MPS_agT~I_Error, alldat)
#*Very sig at .00059
cora<-cor(alldat$MPS_agT, alldat$I_Error, use="complete.obs")

#MPQ stress reaction
gm1<-lm(MPS_srT~I_MeanRT, alldat)
gm1<-lm(MPS_srT~I_Error, alldat)

#mPQ alienation
gm1<-lm(MPS_alT~I_MeanRT, alldat)
gm1<-lm(MPS_alT~I_Error, alldat)
# sig at .08
cora<-cor(alldat$MPS_alT, alldat$I_Error, use="complete.obs")


#ATQ-ECTotal
gm1<-lm(ATQ_ECTotal~I_MeanRT, alldat)
#*Sig at .0476
cora<-cor(alldat$ATQ_ECTotal, alldat$I_MeanRT, use="complete.obs")
gm1<-lm(ATQ_ECTotal~I_Error, alldat)


#ATQ-InhibControl
gm1<-lm(ATQ_InhibControl~I_MeanRT, alldat)
gm1<-lm(ATQ_InhibControl~I_Error, alldat)

#ATQ-AttControl
gm1<-lm(ATQ_AttControl~I_MeanRT, alldat)
#*Sig at .0133
cora<-cor(alldat$ATQ_AttControl, alldat$I_MeanRT, use="complete.obs")

gm1<-lm(ATQ_AttControl~I_Error, alldat)

#correlations
cora<-cor(alldat$MPS_agT, alldat$I_Error, use="complete.obs")




cora<-cor(alldat$MPS_agT, alldat$I_Error, use="complete.obs")