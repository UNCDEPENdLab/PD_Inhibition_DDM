library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
source('/Users/michael/Dropbox/Hallquist_K01/Products/PD Inhibition Manuscript 2013/PD_ Feb92015/Go Nogo/corwithtarget.R')
source('/home/rajchahal/Dropbox/Hallquist_K01/Products/PD Inhibition Manuscript 2013/PD_ Feb92015/Go Nogo/corwithtarget.R')

alldat<-read.table('/home/rajchahal/Dropbox/Hallquist_K01/Products/PD Inhibition Manuscript 2013/PD_ Feb92015/Go Nogo/alldat.csv',sep=",",header=T)
alldat<-within(alldat,{Subject<-factor(Subject)})
library(nFactors)
attach(alldat)
forP <- cbind(T_Propriety, T_RareVirtues, T_Deviance, T_BackDeviance,
              T_NegativeTemperament, T_Detachment, T_SelfHarm,T_Manipulativeness,T_Aggression,
              T_Dependency, T_LowSelfEsteem,T_Mistrust,T_Entitlement, T_EccentricPerceptions,
              T_PositiveTemperament, T_Exhibitionism, T_Disinhibition, T_Impulsivity, T_Workaholism,
              T_SuicideProneness,
              MPS_alT,MPS_agT,MPS_srT, MPS_wbT, MPS_spT, MPS_acT,
              MPS_scT, MPS_clT, MPS_haT, MPS_tdT, MPS_abT, MPS_PET,
              MPS_NET, MPS_COT, MPS_PGT, MPS_PCT, MPS_NGT, MPS_NLT,
              MPS_uvT,
              ATQ_ECTotal,ATQ_InhibControl,ATQ_AttControl, ATQ_AttShifting, ATQ_ActControl,
              zk10, zstai)
#PCA to determine number of factors
library(psych)

c<-na.omit(forP)
df <- data.frame(c)


df$T_NegativeTemperament <- NULL
df$T_PositiveTemperament<- NULL
df$T_Disinhibition<- NULL
df$T_SelfHarm <- NULL
df$T_SuicideProneness <- NULL
#After VIF analysis, dropping additional scales
df$MPS_PGT<-NULL
df$MPS_PET<-NULL
df$MPS_NET<-NULL
df$MPS_COT<-NULL
df$MPS_PCT<-NULL
df$MPS_NGT<-NULL
df$MPS_uvT<-NULL
df$T_RareVirtues<-NULL
df$T_BackDeviance<-NULL
df$T_Deviance<-NULL
df$ATQ_AttControl<-NULL
df$ATQ_AttShifting<-NULL
df$ATQ_InhibControl<-NULL
df$T_RareVirtues<-NULL
df$MPS_NLT<-NULL

#6 fa
fa.promax <- fa(df, nfactors=6, rotate="promax", fm="mle")
fa.diagram1<-fa.diagram(fa.promax$loadings, sort=TRUE, simple=TRUE, cex=.5,digits=1, e.size=.03, rsize=1, side="L" ) 
#take six and five factor solutions ->Lavaan as CFA to look at fit and Modification indeces
library(lavaan)
#specify the model for 6 ffactor sol
cfa.model<-'antagonism=~T_Detachment + MPS_scT+ T_Mistrust+T_Manipulativeness+MPS_alT+T_Aggression+MPS_agT
  neuroticism=~zstai+zk10+MPS_srT+T_Dependency+T_LowSelfEsteem
  constraint=~ T_Propriety+T_Impulsivity+MPS_clT+MPS_tdT
  extroversion=~ T_Exhibitionism+ T_Entitlement+MPS_spT+MPS_wbT
  achievement=~MPS_acT+T_Workaholism+ATQ_ECTotal
  openness=~MPS_abT+T_EccentricPerceptions'
# fit the model
fit <- cfa(cfa.model, data=df)
# display summary output
summary(fit, fit.measures=TRUE)
#got warning message that some indeces are negative


#5 fa
fa.promax <- fa(df, nfactors=5, rotate="promax", fm="mle")
fa.diagram1<-fa.diagram(fa.promax$loadings, sort=TRUE, simple=TRUE, cex=.5,digits=1, e.size=.03, rsize=1, side="L" ) 
#take  five factor solutions ->Lavaan as CFA to look at fit and Modification indeces
library(lavaan)
#specify the model for 5 ffactor sol
cfa.model2<-'antagonism=~T_Detachment + MPS_scT+ T_Mistrust+T_Manipulativeness+MPS_alT+T_Aggression+MPS_agT
  neuroticism=~zstai+zk10+MPS_srT+T_Dependency+T_LowSelfEsteem
  constraint=~ T_Propriety+T_Impulsivity+MPS_clT+MPS_tdT
  extro_open=~ T_Exhibitionism+ T_Entitlement+MPS_spT+MPS_wbT+MPS_abT+T_EccentricPerceptions
  achievement=~MPS_acT+T_Workaholism+ATQ_ECTotal'
 
# fit the model
fit <- cfa(cfa.model2, data=df)
# display summary output
summary(fit, fit.measures=TRUE)
standardizedSolution(fit, type="std.all")
#check MIS and go big ones first to check
mi=modindices(fit)
print(mi[which(mi$mi>10),])
#neuroticism=~T_workaholism Mi= 30.028
#T_aggression=~MPS_aggresion Mi=46.738
#T_Detachment~~MPS_scT mi=59.367
alldat$Agg_Total<-(alldat$T_Aggression+alldat$MPS_agT)/2


cfa.model3<-'antagonism=~T_Detachment + MPS_scT+ T_Mistrust+T_Manipulativeness+MPS_alT+Agg_Total
  neuroticism=~zstai+zk10+MPS_srT+T_Dependency+T_LowSelfEsteem
  constraint=~ T_Propriety+T_Impulsivity+MPS_clT+MPS_tdT
  extro_open=~ T_Exhibitionism+ T_Entitlement+MPS_spT+MPS_wbT+MPS_abT+T_EccentricPerceptions
  achievement=~MPS_acT+T_Workaholism+ATQ_ECTotal'
# fit the model
fit <- cfa(cfa.model3, data=df)
# display summary output
summary(fit, fit.measures=TRUE)
standardizedSolution(fit, type="std.all")
#check MIS and go big ones first to check
mi=modindices(fit)
print(mi[which(mi$mi>10),])
attach(alldat)
forP <- cbind(T_Propriety, T_RareVirtues, T_Deviance, T_BackDeviance,
              T_NegativeTemperament, T_Detachment, T_SelfHarm,T_Manipulativeness,
              T_Dependency, T_LowSelfEsteem,T_Mistrust,T_Entitlement, T_EccentricPerceptions,
              T_PositiveTemperament, T_Exhibitionism, T_Disinhibition, T_Impulsivity, T_Workaholism,
              T_SuicideProneness,Agg_Total,
              MPS_alT,MPS_srT, MPS_wbT, MPS_spT, MPS_acT,
              MPS_scT, MPS_clT, MPS_haT, MPS_tdT, MPS_abT, MPS_PET,
              MPS_NET, MPS_COT, MPS_PGT, MPS_PCT, MPS_NGT, MPS_NLT,
              MPS_uvT,
              ATQ_ECTotal,ATQ_InhibControl,ATQ_AttControl, ATQ_AttShifting, ATQ_ActControl,
              zk10, zstai)
c<-na.omit(forP)
df <- data.frame(c)
df$T_NegativeTemperament <- NULL
df$T_PositiveTemperament<- NULL
df$T_Disinhibition<- NULL
df$T_SelfHarm <- NULL
df$T_SuicideProneness <- NULL
#After VIF analysis, dropping additional scales
df$MPS_PGT<-NULL
df$MPS_PET<-NULL
df$MPS_NET<-NULL
df$MPS_COT<-NULL
df$MPS_PCT<-NULL
df$MPS_NGT<-NULL
df$MPS_uvT<-NULL
df$T_RareVirtues<-NULL
df$T_BackDeviance<-NULL
df$T_Deviance<-NULL
df$ATQ_AttControl<-NULL
df$ATQ_AttShifting<-NULL
df$ATQ_InhibControl<-NULL
df$T_RareVirtues<-NULL
df$MPS_NLT<-NULL
df$T_Aggression<-NULL
df$MPS_agT<-NULL


alldat$Agg_Total<-(alldat$T_Aggression+alldat$MPS_agT)/2


cfa.model3<-'antagonism=~T_Detachment + MPS_scT+ T_Mistrust+T_Manipulativeness+MPS_alT+Agg_Total
neuroticism=~zstai+zk10+MPS_srT+T_Dependency+T_LowSelfEsteem
constraint=~ T_Propriety+T_Impulsivity+MPS_clT+MPS_tdT
extro_open=~ T_Exhibitionism+ T_Entitlement+MPS_spT+MPS_wbT+MPS_abT+T_EccentricPerceptions
achievement=~MPS_acT+T_Workaholism+ATQ_ECTotal'
# fit the model
fit <- cfa(cfa.model3, data=df)
# display summary output
summary(fit, fit.measures=TRUE)
standardizedSolution(fit, type="std.all")
#check MIS and go big ones first to check
mi=modindices(fit)
print(mi[which(mi$mi>10),])
attach(alldat)
forP <- cbind(T_Propriety,
              T_Detachment,T_Manipulativeness,
              T_Dependency, T_LowSelfEsteem,T_Mistrust,T_Entitlement, T_EccentricPerceptions,
              T_Exhibitionism, T_Impulsivity, T_Workaholism,
              Agg_Total,
              MPS_alT,MPS_srT, MPS_wbT, MPS_spT, MPS_acT,
              MPS_scT, MPS_clT, MPS_haT, MPS_tdT, MPS_abT,
              MPS_NET, MPS_NLT,
              ATQ_ECTotal,
              zk10, zstai)
c<-na.omit(forP)
df <- data.frame(c)
df$T_NegativeTemperament <- NULL
df$T_PositiveTemperament<- NULL
df$T_Disinhibition<- NULL
df$T_SelfHarm <- NULL
df$T_SuicideProneness <- NULL
#After VIF analysis, dropping additional scales
df$MPS_PGT<-NULL
df$MPS_PET<-NULL
df$MPS_NET<-NULL
df$MPS_COT<-NULL
df$MPS_PCT<-NULL
df$MPS_NGT<-NULL
df$MPS_uvT<-NULL
df$T_RareVirtues<-NULL
df$T_BackDeviance<-NULL
df$T_Deviance<-NULL
df$ATQ_AttControl<-NULL
df$ATQ_AttShifting<-NULL
df$ATQ_InhibControl<-NULL
df$T_RareVirtues<-NULL
df$MPS_NLT<-NULL
df$T_Aggression<-NULL
df$MPS_agT<-NULL

#do means of the two scores to extract one score. might help the model fit

#got warning message that some indeces are negative

#read 5,6,7 in book on fit indeces and meaurement modeling
