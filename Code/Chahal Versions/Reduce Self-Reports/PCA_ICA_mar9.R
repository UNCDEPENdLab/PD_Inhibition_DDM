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
forP.pca<-princomp(na.omit(forP))
summary(forP.pca)
plot(forP.pca)
c<-na.omit(forP)
#dat.fa1<-factanal(c, factors=10, rotation="varimax")
#fa1<-fa(r=cor(na.omit(forP)), nfactors=10, rotate="varimax", SMC=FALSE, fm="minres")
cor(c)
pca1<-princomp(c, scores=T, cor=T)
summary(pca1)
#Eigen vals are close to 1 up until Comp. 11- then prop of variance is only 2%)
loadings(pca1)
#pca1$locadings
#scree plot of Eigen vals
plot(pca1)
screeplot(pca1, type="line", main="Scree Plot") #cut off at ~1 for var
#biplot
biplot(pca1)
#scores of components
pca1$scores
df <- data.frame(c)
form <- formula(paste("dummy ~ ",paste(strsplit(names(df)," "),collapse=" + ")))
df$dummy <- 1
##act control and suicide proneness are aliased (collinear)
m <- lm(form, df)
summary(m)
#car::vif(m) #won't run with aliasing

##take out ATQ total score since the sum of the ATQ subscales is equal to the total
##take out self-harm because it is the sum of SNAP low self-esteem and suicide proneness
df$ATQ_ECTotal <- NULL
df$T_SelfHarm <- NULL

form <- formula(paste("dummy ~ ",paste(strsplit(names(df)," "),collapse=" + ")))
df$dummy <- 1
##act control and suicide proneness are aliased (collinear)
m <- lm(form, df)

df$dummy <- NULL #take out constate column prior to factor analysis

fa1<-factanal(df, factors=8, rotation="promax", scores="regression")

library(psych)
##still a heywood case (likely due to collinearity, even if not perfect as above)
##I would try to take out any superordinate factor that represents sums of subscales (e.g., neg temperament, pos temperament).

#3/1/15
df$T_NegativeTemperament <- NULL
df$T_PositiveTemperament<- NULL
df$T_Disinhibition<- NULL
fa1<-factanal(df, factors=8, rotation="promax", scores="regression")

##Start here to check colinearitiess with any list
library(car)
c<-na.omit(forP)
df <- data.frame(c)
df$dummy <- 1
form <- formula(paste("dummy ~ ",paste(strsplit(names(df)," "),collapse=" + ")))
Mod<-lm(form,data=df)
vif.Mod<-vif(Mod)
#generally a vif above 5 is considered too high multilinearity/ should exclude or do stepwise regression
#identify aliasing culprits
# Choose a VIF cutoff under which a variable is retained (Zuur et al. 2010 
# MEE recommends 2)
library(fmsb)
library(MASS)
vif_func(df, thresh=5, trace=T)

# Look at the final model
keep.dat<-vif_func(df,thresh=5,trace=F)
form <- formula(paste("dummy ~ ",paste(keep.dat)," "),collapse=" + ")
form <- formula(paste("dummy ~ ",paste(strsplit(names(df)," "),collapse=" + ")))
form.in<-paste("dummy ~",paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=df)
summary(mod2)


#when dropping columns:
df$T_NegativeTemperament <- NULL
df$T_PositiveTemperament<- NULL
df$T_Disinhibition<- NULL
#df$ATQ_ECTotal <- NULL
df$T_SelfHarm <- NULL
df$T_SuicideProneness <- NULL
df$ATQ_ActControl <- NULL
#After VIF analysis, dropping additional scales
df$MPS_PGT<-NULL
df$MPS_PET<-NULL
df$MPS_NET<-NULL
df$MPS_COT<-NULL
df$MPS_PCT<-NULL
df$MPS_NGT<-NULL

##Mar 9, df2
df$MPS_uvT<-NULL
df$T_RareVirtues<-NULL
df$T_BackDeviance<-NULL
df$T_Deviance<-NULL
df$ATQ_AttControl<-NULL
df$ATQ_AttShifting<-NULL
df$ATQ_InhibControl<-NULL
#Mar 9, df3
df$T_RareVirtues<-NULL
df$MPS_NLT<-NULL
#Only use ATQ total(C) and drop subscale. pull out the big ones - should be in ML2 (from 6 factor diag)
#take out rare virtues and uvt


##From firsst pass: dont drop following
#After dropping all of the top values, factanal works without colinearity error
df$MPS_NLT<-NULL
df$ATQ_ECTotal<-NULL
df$T_SelfHarm<-NULL
df$T_Disinhibition<-NULL
df$MPS_srT<-NULL
df$T_Workaholism<-NULL
df$T_PositiveTemperament<-NULL
df$T_Detachment<-NULL

fa1<-factanal(df, factors=8, rotation="promax", scores="regression",method="mle")


##MAR 9 2015
#Fa in psych pacakage should give BIC to give best fit/number of factors
#RMSEA in ml based analysis- want to be .07or below.
fa.promax <- fa(df, nfactors=6, rotate="promax", fm="mle")
fa.diagram1<-fa.diagram(fa.promax$loadings, sort=TRUE, simple=TRUE, cex=.5,digits=1, e.size=.03, rsize=1, side="L" ) 
#dROPPING t_rarevirtues, t_deviance, t_backdeviance?, (MPS_uvT??)
###elbow plot of bic.. 6 is min or 5 is alt in terms of fit vs parimony tradeoff

#try with diff number of factors, get BIC and put in excel
#RMSEA goes down with more factors
#TLI goes up with more factors-want .9or above
#Set up SEM with Lavaan :modifcation indeces tell how much better model gets with allowed individual association
#then maybe use factor scores to give us a new heat map

# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation 
fit <- factanal(df, 8, rotation="promax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(df),cex=.7) # add variable names

##drop
fa.diagram1<-fa.diagram(fa.promax$loadings, sort=TRUE, simple=TRUE, digits=1, e.size=.05, rsize=1.5) 
fa.diagram2<-fa.diagram(fit$loadings, digits=1, e.size=.05,cut=.3, rsize=.1, side=2, cex=NULL)
##prune if loading onto mult factors weakly
dfcor<-cor(df)
fa.promax <- fa(df, nfactors=7, rotate="promax", fm="mle")
fa.diagram(fa.promax, simple=FALSE)
source("http://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")
#
install.packages("Rgraphviz")
library("Rgraphviz")
fa.rgraph(fa.promax$loadings,simple=TRUE,
         size=c(8,6), node.font=c("Helvetica", 14),
         edge.font=c("Helvetica", 10), 
         digits=1,main="Factor Analysis", rank.direction=c("LR"))
#
f3 <- fa(df,8,rotate="cluster")
fa.diagram(f3,side=3)

#
#fa5<-factanal(df, factors=4, rotation="promax", scores="regression")
structure.diagram(fa1, Phi=NULL,fy=NULL,labels=NULL,cut=.3,errors=FALSE,simple=TRUE,
                  regression=FALSE,lr=TRUE,Rx=NULL,Ry=NULL,digits=1,
                  main="Structural model")
#
structure.graph(fa1,  Phi = NULL,fy = NULL, out.file = NULL, labels = NULL, cut = 0.3, 
                errors=FALSE, simple=TRUE,regression=FALSE, size = c(8, 6), 
                node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), 
                rank.direction = c("RL"), digits = 1, 
                title = "Structural model")
ll<-structure.sem(fa1,  Phi = NULL, fy = NULL,out.file = NULL, labels = NULL,
              cut = 0.3, errors=TRUE, simple=TRUE,regression=FALSE)


fadaf <- fa(df, nfactors=4, rotate="oblimin", SMC=TRUE)
#still getting comp sing error
#when rotation=varimax, this code worked?.need to try promax
#because of error, going to try correlation matrix to figure out whhich scores not to include that are
#overlapping- for example subssscales. 
cor(c)
splom(cor(c))

corrs<-cor(c)
write.table(corrs, file="feb09_corrs.csv", col.names=T, sep=",")
#Michael's code
library(nFactors)

ev <- eigen(cor(c, use="pairwise.complete.obs")) # get eigenvalues  
ap.fa <- parallel(subject=nrow(c), var=ncol(c), rep=2000, quantile=.05, model="factors")
nS.fa <- nScree(x=ev$values, aparallel=ap.fa$eigen$qevpea, model="factors")
numRetain.fa <- nS.fa$Components$nparallel
#says to retain 10 which is same as PCA above
#plot(nS.fa)
library(psych)
f1 <- fa(x[,-1], fm="ml", nfactors=1, missing=TRUE, rotate="oblimin", SMC=TRUE) #oblique
f2 <- fa(x[,-1], fm="ml", nfactors=2, missing=TRUE, rotate="oblimin", SMC=TRUE)
f3 <- fa(x[,-1], fm="ml", nfactors=3, missing=TRUE, rotate="oblimin", SMC=TRUE)
f4 <- fa(c, fm="ml", nfactors=10, missing=TRUE, rotate="oblimin", SMC=TRUE)
#confused- didnt work
f4 <- fa(c, fm="ml", nfactors=10, missing=TRUE, rotate="oblimin", SMC=TRUE)
#The warnings and errors indicates that your matrix is singular, thus no solution exists to the optimization problem.

#This means you need to use a different method of factor analysis. Using fa() in package psych you have two alternatives to perform factor analysis given a singular matrix:
  
 # pa (Principal axis factor analysis)
#minres (Minimum residual factor analysis)
library(GPArotation)
z<-fa(r=cor(c), nfactors=4, rotate="varimax", SMC=FALSE, fm="minres")
loadings(z)
## Check here which vars to include in each factor then do SEM
#Another method
fit <- princomp(c, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# PCA Variable Factor Map 
library(FactoMineR)
result <- PCA(c) # graphs generated automatically

# Simple CFA Model
library(sem)
mydata.cov <- cov(c)
model.mydata <- specify.model() 
F1 ->  X1, lam1, NA
F1 ->  X2, lam2, NA 
F1 ->  X3, lam3, NA 
F2 ->  X4, lam4, NA 
F2 ->  X5, lam5, NA 
F2 ->  X6, lam6, NA 
X1 <-> X1, e1,   NA 
X2 <-> X2, e2,   NA 
X3 <-> X3, e3,   NA 
X4 <-> X4, e4,   NA 
X5 <-> X5, e5,   NA 
X6 <-> X6, e6,   NA 
F1 <-> F1, NA,    1 
F2 <-> F2, NA,    1 
F1 <-> F2, F1F2, NA
#confused on how to set which vars go where
