#Linear Models for determining predictions of performance on task with self reports 
library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(reshape)
library(multcomp)

alldat<-read.table('/home/rajchahal/Desktop/WPIC/alldat.csv',sep=",",header=T)

#escore<-merge(AccScore, alldat,  all=T, by="Subject")
alldat<-within(alldat,{Subject<-factor(Subject)})
#Using RT and Err on incongruent trials flank control as DV's

library(MASS)
#for poor accuracy outliers could use robust regression (e.g., rlm in MASS) 
alldat$FlankAcc <- 1-(alldat$FlankAcc)
alldat$Flankpinac <- alldat$FlankAcc*100

#RT incongruent Trials correlations: I_MeanRt and I_Error
#Neg Temperament
gm1 <- lm(T_NegativeTemperament ~ I_MeanRT, alldat)
gm1<-lm(T_NegativeTemperament~ Flankpinac, alldat)


#and Impulsivity
gm1 <- lm(T_Impulsivity ~ I_MeanRT, alldat)
gm1<-lm(T_Impulsivity~ Flankpinac, alldat)

#Propriety
gm1 <- lm(T_Propriety ~ I_MeanRT, alldat)
gm1<-lm(T_Propriety~ Flankpinac, alldat)
#*Barely sig at .0613
cora<-cor(alldat$T_Propriety, alldat$Flankpinac, use="complete.obs")

#Detachment
gm1 <- lm(T_Detachment ~ I_MeanRT, alldat)
gm1<-lm(T_Detachment~ Flankpinac, alldat)
#*Sig at .00279
cora<-cor(alldat$T_Detachment, alldat$Flankpinac, use="complete.obs")


#Entitlement
gm1 <- lm(T_Entitlement ~ I_MeanRT, alldat)
gm1<-lm(T_Entitlement~ Flankpinac, alldat)


#Self Harm
gm1 <- lm(T_SelfHarm ~ I_MeanRT, alldat)
gm1<-lm(T_SelfHarm~ Flankpinac, alldat)
#*Sig at .0382
cora<-cor(alldat$T_SelfHarm, alldat$Flankpinac, use="complete.obs")


#Aggression
gm1 <- lm(T_Aggression ~ I_MeanRT, alldat)
gm1<-lm(T_Aggression~ Flankpinac, alldat)
#*Very sig at 5.46e-05
cora<-cor(alldat$T_Aggression, alldat$Flankpinac, use="complete.obs")

#Dependency
gm1 <- lm(T_Dependency ~ I_MeanRT, alldat)
gm1<-lm(T_Dependency~ Flankpinac, alldat)
#*Sig at .0773
cora<-cor(alldat$T_Dependency, alldat$Flankpinac, use="complete.obs")

#Mistrust
gm1 <- lm(T_Mistrust ~ I_MeanRT, alldat)
gm1<-lm(T_Mistrust~ Flankpinac, alldat)
#**Super sig at 5.66*10-5
cora<-cor(alldat$T_Mistrust, alldat$Flankpinac, use="complete.obs")

#Manipulativeness
gm1 <- lm(T_Manipulativeness ~ I_MeanRT, alldat)
gm1<-lm(T_Manipulativeness~ Flankpinac, alldat)
#*Very sig at .00146
cora<-cor(alldat$T_Manipulativeness, alldat$Flankpinac, use="complete.obs")

# Low Self Esteem
gm1 <- lm(T_LowSelfEsteem ~ I_MeanRT, alldat)
gm1<-lm(T_LowSelfEsteem~ Flankpinac, alldat)

#MPQ neg emotionality
gm1 <- lm(MPS_NET ~ I_MeanRT, alldat)
gm1<-lm(MPS_NET~ Flankpinac, alldat)
#*Barely sig at .0576
cora<-cor(alldat$MPS_NET, alldat$Flankpinac, use="complete.obs")

#MPQ alienation
gm1 <- lm(MPS_alT ~ I_MeanRT, alldat)
gm1<-lm(MPS_alT~ Flankpinac, alldat)


#MPQ aggression
gm1<-lm(MPS_agT~I_MeanRT, alldat)
gm1<-lm(MPS_agT~Flankpinac, alldat)
#*Very sig at .00435
cora<-cor(alldat$MPS_agT, alldat$Flankpinac, use="complete.obs")

#MPQ stress reaction
gm1<-lm(MPS_srT~I_MeanRT, alldat)
gm1<-lm(MPS_srT~Flankpinac, alldat)

#mPQ alienation
gm1<-lm(MPS_alT~I_MeanRT, alldat)
gm1<-lm(MPS_alT~Flankpinac, alldat)



#ATQ-ECTotal
gm1<-lm(ATQ_ECTotal~I_MeanRT, alldat)
#*Sig at .0476
cora<-cor(alldat$ATQ_ECTotal, alldat$I_MeanRT, use="complete.obs")
gm1<-lm(ATQ_ECTotal~Flankpinac, alldat)


#ATQ-InhibControl
gm1<-lm(ATQ_InhibControl~I_MeanRT, alldat)
gm1<-lm(ATQ_InhibControl~Flankpinac, alldat)

#ATQ-AttControl
gm1<-lm(ATQ_AttControl~I_MeanRT, alldat)
#*Sig at .0133
cora<-cor(alldat$ATQ_AttControl, alldat$I_MeanRT, use="complete.obs")

gm1<-lm(ATQ_AttControl~I_Error, alldat)

#correlations
cora<-cor(alldat$MPS_agT, alldat$I_Error, use="complete.obs")




cora<-cor(alldat$MPS_agT, alldat$I_Error, use="complete.obs")


cm<-lmerCellMeans(fm1)
d <- ggplot(alldat, aes(x=, y=Flankpinac, 
                    color=Block, 
                    )) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Block", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(d)