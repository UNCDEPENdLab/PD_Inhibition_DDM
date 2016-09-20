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

#Drop eflanker outlier subjects
alldate<-subset(alldat, Subject!="98" & Subject!="39" & Subject!="50" & Subject!="61" & Subject!="96")

library(MASS)
#for poor accuracy outliers could use robust regression (e.g., rlm in MASS) 

#Sad Trials correlations
#Sad ACC and Neg Temperament
gm1 <- lm(T_NegativeTemperament ~ Sad_I_PctInaccurate, alldate)
#**Significant at .00915 # How to plot correlation matrix?
g<-ggplot(gm1, aes(x=Sad_I_PctInaccurate, y=T_NegativeTemperament))+
cora<-cor(alldate$T_NegativeTemperament, alldate$Sad_I_PctInaccurate, use="complete.obs")

  #Sad RT and Neg Temp
gm1<-lm(T_NegativeTemperament~ I_Sad, alldate)


#Sad ACC and Impulsivity
gm1 <- lm(T_Impulsivity ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Impulsivity~ I_Sad, alldate)

#Sad ACC and Propriety
gm1 <- lm(T_Propriety ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Propriety~ I_Sad, alldate)

#Sad ACC and Detachment
gm1 <- lm(T_Detachment ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Detachment~ I_Sad, alldate)

#Sad ACC and Entitlement
gm1 <- lm(T_Entitlement ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Entitlement~ I_Sad, alldate)
# * SIGNIFICANT AT .0471
cora<-cor(alldate$T_Entitlement, alldate$I_Sad, use="complete.obs")


#Sad ACC and Self Harm
gm1 <- lm(T_SelfHarm ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_SelfHarm~ I_Sad, alldate)

#Sad ACC and Aggression
gm1 <- lm(T_Aggression ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Aggression~ I_Sad, alldate)
  
#Sad ACC and Dependency
gm1 <- lm(T_Dependency ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_Dependency~ I_Sad, alldate)

#Sad ACC and Mistrust
gm1 <- lm(T_Mistrust ~ Sad_I_PctInaccurate, alldate)
#*** Significant at .0449
cora<-cor(alldate$T_Mistrust, alldate$Sad_I_PctInaccurate, use="complete.obs")

#Sad RT and Mistrust
gm1<-lm(T_Mistrust~ I_Sad, alldate)
#** Significant at .0423
cora<-cor(alldate$T_Mistrust, alldate$I_Sad, use="complete.obs")


#Sad ACC and Manipulativeness
gm1 <- lm(T_Manipulativeness ~ Sad_I_PctInaccurate, alldate)
#*Significant, barely at .0156
cora<-cor(alldate$T_Manipulativeness, alldate$Sad_I_PctInaccurate, use="complete.obs")

gm1<-lm(T_Manipulativeness~ I_Sad, alldate)

#Sad ACC and Low Self Esteem
gm1 <- lm(T_LowSelfEsteem ~ Sad_I_PctInaccurate, alldate)
gm1<-lm(T_LowSelfEsteem~ I_Sad, alldate)

#sad acc and MPQ neg emotionality
gm1 <- lm(MPS_NET ~ Sad_I_PctInaccurate, alldate)
#**SIGnificant at .00155
cora<-cor(alldate$MPS_NET, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_NET~ I_Sad, alldate)

#Sad acc and MPQ alienation
gm1 <- lm(MPS_alT ~ Sad_I_PctInaccurate, alldate)
#** Ver Sig at .00211
cora<-cor(alldate$MPS_alT, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_alT~ I_Sad, alldate)

#Sad acc and MPQ aggression
gm1<-lm(MPS_agT~Sad_I_PctInaccurate, alldate)
#*significant at .064
cora<-cor(alldate$MPS_agT, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_agT~I_Sad, alldate)

#Sad acc and MPQ stress reaction
gm1<-lm(MPS_srT~Sad_I_PctInaccurate, alldate)
#*** Very significant at .00664
cora<-cor(alldate$MPS_srT, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_srT~I_Sad, alldate)

#Sad acc and ATQ-ECTotal
gm1<-lm(ATQ_ECTotal~Sad_I_PctInaccurate, alldate)
#*Sig at .0188
cora<-cor(alldate$ATQ_ECTotal, alldate$Sad_I_PctInaccurate, use="complete.obs")

gm1<-lm(ATQ_ECTotal~I_Sad, alldate)




#Sad acc and ATQ-InhibControl
gm1<-lm(ATQ_InhibControl~Sad_I_PctInaccurate, alldate)
#*Barely sig at .0829
cora<-cor(alldate$ATQ_InhibControl, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(ATQ_InhibControl~I_Sad, alldate)

#Sad acc and ATQ-AttControl
gm1<-lm(ATQ_AttControl~Sad_I_PctInaccurate, alldate)
#*Sig at .0681
cora<-cor(alldate$ATQ_AttControl, alldate$Sad_I_PctInaccurate, use="complete.obs")
gm1<-lm(ATQ_AttControl~I_Sad, alldate)

#Ang Trials correlations
#Ang ACC and Neg Temperament
gm1 <- lm(T_NegativeTemperament ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_NegativeTemperament~ I_Ang, alldate)


#Ang ACC and Impulsivity
gm1 <- lm(T_Impulsivity ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Impulsivity~ I_Ang, alldate)

#Ang ACC and Propriety
gm1 <- lm(T_Propriety ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Propriety~ I_Ang, alldate)

#Ang ACC and Detachment
gm1 <- lm(T_Detachment ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Detachment~ I_Ang, alldate)

#Ang ACC and Entitlement
gm1 <- lm(T_Entitlement ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Entitlement~ I_Ang, alldate)
# * SIGNIFICANT AT .0531
cora<-cor(alldate$T_Entitlement, alldate$I_Ang, use="complete.obs")


#Ang ACC and Self Harm
gm1 <- lm(T_SelfHarm ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_SelfHarm~ I_Ang, alldate)

#Ang ACC and Aggression
gm1 <- lm(T_Aggression ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Aggression~ I_Ang, alldate)

#Ang ACC and Dependency
gm1 <- lm(T_Dependency ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Dependency~ I_Ang, alldate)

#Ang ACC and Mistrust
gm1 <- lm(T_Mistrust ~ Ang_I_PctInaccurate, alldate)
#*** Significant at .0641
cora<-cor(alldate$T_Mistrust, alldate$Ang_I_PctInaccurate, use="complete.obs")
gm1<-lm(T_Mistrust~ I_Ang, alldate)

#Ang ACC and Manipulativeness
gm1 <- lm(T_Manipulativeness ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Manipulativeness~ I_Ang, alldate)

#Ang ACC and Low Self Esteem
gm1 <- lm(T_LowSelfEsteem ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(T_LowSelfEsteem~ I_Ang, alldate)

#Ang acc and MPQ neg emotionality
gm1 <- lm(MPS_NET ~ Ang_I_PctInaccurate, alldate)
#*SIGnificant at .0942
cora<-cor(alldate$MPS_NET, alldate$Ang_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_NET~ I_Ang, alldate)

#Ang acc and MPQ alienation
gm1 <- lm(MPS_alT ~ Ang_I_PctInaccurate, alldate)
gm1<-lm(MPS_alT~ I_Ang, alldate)

#Ang acc and MPQ aggression
gm1<-lm(MPS_agT~Ang_I_PctInaccurate, alldate)
#*significant at .064
cora<-cor(alldate$MPS_agT, alldate$Ang_I_PctInaccurate, use="complete.obs")
gm1<-lm(MPS_agT~I_Sad, alldate)

#Ang acc and MPQ stress reaction
gm1<-lm(MPS_srT~Ang_I_PctInaccurate, alldate)
gm1<-lm(MPS_srT~I_Ang, alldate)

#Ang and MPQ alienation
gm1<-lm(MPS_alT~Ang_I_PctInaccurate, alldate)
gm1<-lm(MPS_alT~I_Ang, alldate)

#Ang acc and ATQ-ECTotal
gm1<-lm(ATQ_ECTotal~Ang_I_PctInaccurate, alldate)
gm1<-lm(ATQ_ECTotal~I_Ang, alldate)

#Ang acc and ATQ-InhibControl
gm1<-lm(ATQ_InhibControl~Ang_I_PctInaccurate, alldate)
gm1<-lm(ATQ_InhibControl~I_Ang, alldate)

#Ang acc and ATQ-AttControl
gm1<-lm(ATQ_AttControl~Ang_I_PctInaccurate, alldate)
gm1<-lm(ATQ_AttControl~I_Ang, alldate)

#Ang acc and eccentric perceps
gm1<-lm(T_EccentricPerceptions~Ang_I_PctInaccurate, alldate)
gm1<-lm(T_EccentricPerceptions~I_Ang, alldate)

#Ang acc and deviance
gm1<-lm(T_Deviance~Ang_I_PctInaccurate, alldate)
gm1<-lm(T_Deviance~I_Ang, alldate)


jj




















#Angry Trials correlations  
gm1 <- lm(T_NegativeTemperament ~ Ang_I_PctInaccurate, alldate)
#Not Significant
  
