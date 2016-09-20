library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
source('/home/rajchahal/Desktop/WPIC/flanker_Control/corwithtarget.R')

alldat<-read.table('/home/rajchahal/Desktop/WPIC/alldat.csv',sep=",",header=T)
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
fa1<-factanal(c, factors=8, rotation="promax", scores="regression")
#still getting comp sing error
#when rotation=varimax, this code worked?. according to michael usse promax since this doesnt
#need orthonality
#because of error, going to try correlation matrix to figure out whhich scores not to include that are
#overlapping- for example subssscales. 
cor(c)
splom(cor(c))
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