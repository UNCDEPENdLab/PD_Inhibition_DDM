library(sas7bdat)
library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(Hmisc)
library(multcomp)
library(reshape)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
dat<-read.sas7bdat("/home/rajchahal/Dropbox/Hallquist_K01/Products/PD Inhibition Manuscript 2013/SAS Originals/recentprobes.sas7bdat")
# Variables of Interest: ProbeDisplay_RT, ProbeDisplay_ACC, Running_Block_ Running_Trial_?, Subject,Condition, 
# Condition= Nfam0, Nfam1, Nfam1RI, Nfam2, Y
#Y= positive trials- exclude eventually
#Nfam0= Negative Unfamiliar
# Nfam1= Negative Familiar (was in trial preceding)
# Nfam2= Negative Highly Familiary (was in two preceding trials)
# Nfam1RI= Negative Familiar Response Interferencer (was correct probe in preceding trial)

dat<-within(dat, {Subject<-factor(Subject)})
dat<-within(dat, {Condition<-factor(Condition,
                                           levels=c("Nfam0",
                                                    "Nfam1",
                                                    "Nfam2",
                                                    "Nfam1RI",
                                                    "Y"),
                                           labels=c("Neg Unfam",
                                                    "Neg Fam",
                                                    "Neg High Fam",
                                                    "Neg Fam RI",
                                                    "Pos"))})
#Subset out Positive trials and subset out RT=0 trials and Subset out the low acc people

#Drop RT = 0
dat<-subset(dat, ProbeDisplay_RT!="0")
dat<-subset(dat, ProbeDisplay_RT!="54")
#Will need to trim RT by condition also

#Drop Condition=Pos here?

#Group descriptives to find outliers on ACC using means
getDescriptivesByGroup <- function(df, outcomes, group, zscale=FALSE, tscale=FALSE) {
  require(plotrix)
  oneVar <- function(df, outcome) {
    stats <- c("mean", "sd", "std.error", "median", "min", "max")
    bigMatrix <- c()
    
    for (stat in stats) {
      if (zscale) {
        bigMatrix <- cbind(bigMatrix, tapply(as.vector(scale(df[[outcome]])), df[[group]], stat, na.rm=TRUE))  
      } else if (tscale) {
        scaleVec <- as.vector(scale(df[[outcome]]))
        scaleVec <- 50 + 10*scaleVec
        bigMatrix <- cbind(bigMatrix, tapply(scaleVec, df[[group]], stat, na.rm=TRUE))
      } else {
        bigMatrix <- cbind(bigMatrix, tapply(df[[outcome]], df[[group]], stat, na.rm=TRUE))
      }
      
      colnames(bigMatrix)[ncol(bigMatrix)] <- stat
    }
    
    bigMatrix <- as.data.frame(bigMatrix)
    bigMatrix$var <- outcome
    bigMatrix$group <- factor(rownames(bigMatrix))
    rownames(bigMatrix) <- NULL
    
    return(bigMatrix)    
  }
  
  allVars <- c()
  for (outcome in outcomes) {
    allVars <- rbind(allVars, oneVar(df, outcome))
  }
  allVars <- plyr::rename(allVars, c(group=group))
  return(allVars)
}
#end descriptives function
#type in vars of interest for descriptives table
desc<-getDescriptivesByGroup(dat, c("ProbeDisplay_ACC", "ProbeDisplay_RT"), "Subject")
accdesc<-subset(desc, var=="ProbeDisplay_ACC")
attach(accdesc)
boxplot(mean)
print(boxplot(mean))

#Outlier Acc Subs: .54, .52, .59, (overall). Maybe rerun this only using negative trials
dat<-subset(dat, Subject!="40" & Subject!="65" & Subject!="100")

rtdesc<-subset(desc, var=="ProbeDisplay_RT")
attach(rtdesc)
boxplot(mean)
print(boxplot(mean))
#doesnt look like there are outliers in RT

#subset out the pos trials
dat<-subset(dat, Condition!="Pos")


#RT Observed
g<-ggplot(dat, aes(x=Condition, y=ProbeDisplay_RT), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)
#Acc Obs
g<-ggplot(dat, aes(x=Condition, y=ProbeDisplay_ACC, ), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)
#RT Obs
i<-ggplot(dat, aes(x=Condition, y=ProbeDisplay_RT))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Recent Probes RT")+xlab("Condition\n")+ylab("Mean RTs")
i+ggtitle("Recent Probes RT")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))


#Plot Observed ACC in pretty graph
dd<-ggplot(dat, aes(x=Condition, y=ProbeDisplay_ACC))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer(palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Recent Probes Acc")+xlab("Condition\n")+ylab("Mean ACC (%)")
dd+ggtitle("Recent Probes ACC")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))

#+scale_fill_manual(values=cbPalette)

#General Linear Modelss
m2 <- glmer(ProbeDisplay_ACC ~ Condition+ (1|Subject), data=dat, family="binomial")
m<-lmerCellMeans(m2)

fm2<-lmer(ProbeDisplay_RT~Condition+(1|Subject), data=dat)
cm<-lmerCellMeans(fm2)


#Acc GLM


# Heatmap Code

library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
source('/home/rajchahal/Desktop/WPIC/flanker_Control/corwithtarget.R')

alldat<-read.table('/home/rajchahal/Desktop/WPIC/alldat.csv',sep=",",header=T)
#Vars: RP_Nfam1_ACC = Recent Probes Negative Familiar Acc
#Vars 2: RP_Nfam2_ACC = Recent Probes Negative Highly Familiar Acc
#Vars 3: RP_Nfam1RI_ACC = Recent Probes Familiar Resp Interference Acc
#Vars 4: RP_Nfam1_RT = Recent Probes Negative Familiar RT
#Vars 5: RP_Nfam2_RT = Recent Probes Negative Highly Familiar RT
#Vars 6: RP_Nfam1RI_RT = Recent Probes Familiar Resp Interference RT

alldat<-subset(alldat, Subject!="40" & Subject!="65" & Subject!="100")


alldat<-within(alldat,{Subject<-factor(Subject)})
#alldat$FlankAcc <- 1-(alldat$FlankAcc)
#alldat$Flankpinac <- alldat$FlankAcc*100
#alldat$Ang_I_PctInaccurate<-alldat$Ang_I_PctInaccurate*100
#alldat$Sad_I_PctInaccurate<-alldat$Sad_I_PctInaccurate*100
#dropped these scores manually from eflank pct inacc and rt on sheet
#alldat<-subset(alldat, Subject!="98" & Subject!="39" & Subject!="50" & Subject!="61" & Subject!="96")
traits <- c("T_Propriety", "T_RareVirtues", "T_Deviance", "T_BackDeviance",
            "T_NegativeTemperament", "T_Detachment", "T_SelfHarm","T_Manipulativeness","T_Aggression",
            "T_Dependency", "T_LowSelfEsteem","T_Mistrust","T_Entitlement", "T_EccentricPerceptions",
            "T_PositiveTemperament", "T_Exhibitionism", "T_Disinhibition", "T_Impulsivity", "T_Workaholism",
            "T_SuicideProneness",
            "MPS_alT","MPS_agT","MPS_srT", "MPS_wbT", "MPS_spT", "MPS_acT",
            "MPS_scT", "MPS_clT", "MPS_haT", "MPS_tdT", "MPS_abT", "MPS_PET",
            "MPS_NET", "MPS_COT", "MPS_PGT", "MPS_PCT", "MPS_NGT", "MPS_NLT",
            "MPS_uvT",
            "ATQ_ECTotal","ATQ_InhibControl","ATQ_AttControl", "ATQ_AttShifting", "ATQ_ActControl",
            "zk10", "zstai")
alldat$FamMinHiFamAcc<-(alldat$RP_Nfam1_ACC-alldat$RP_Nfam2_ACC) #Neg Fam- Neg High Fam (ACC)
alldat$FamMinFamRIAcc<-(alldat$RP_Nfam1_ACC-alldat$RP_Nfam1RI_ACC) #Neg Fam- Neg Fam RI (ACC)
alldat$HiFamMinFamRIAcc<-(alldat$RP_Nfam2_ACC-alldat$RP_Nfam1RI_ACC) #Hi Fam- Neg Fam RI (ACC)
alldat$FamMinHiFamRT<-(alldat$RP_Nfam1_RT-alldat$RP_Nfam2_RT) #Neg Fam- Neg High Fam (RT)
alldat$FamMinFamRIRT<-(alldat$RP_Nfam1_RT-alldat$RP_Nfam1RI_RT) #Neg Fam- Neg Fam RI (RT)
alldat$HiFamMinFamRIRT<-(alldat$RP_Nfam2_RT-alldat$RP_Nfam1RI_RT) #Hi Fam- Neg Fam RI (RT)

RecProb <- c("RP_Nfam1_ACC","RP_Nfam2_ACC", "RP_Nfam1RI_ACC","RP_Nfam1_RT", "RP_Nfam2_RT", "RP_Nfam1RI_RT")
RecProb <- c("FamMinHiFamAcc","FamMinFamRIAcc", "HiFamMinFamRIAcc","FamMinHiFamRT", "FamMinFamRIRT", "HiFamMinFamRIRT")

save(alldat,file='alldat.Rdata')
df<-alldat
load('alldat.Rdata')
cmat.list <- corwithtarget(df,target=RecProb,withvars=traits)
cmat.3dmat <- do.call(abind, list(cmat.list, along=0))
cmat.r <- melt( t(cmat.3dmat[,1,]),value.name='r' )
cmat.p <- melt( t(cmat.3dmat[,2,]),value.name='p' )
cmat <- merge(cmat.r,cmat.p, by=c("X1", "X2"))
names(cmat)[1:4] <- c('Trait','RecProb', 'r', 'p')
ge<-ggplot( cmat, aes(x=Trait,y=RecProb,fill=p,label=r) ) + theme_grey(base_size=22)+geom_tile() + scale_fill_continuous(high='#000000',low='#FF0000', limits=c(0,.10) ) + geom_text(color=I('white'), size=5)+theme(axis.title.x=element_text(face="bold",colour="Black", size=14),
                                                                                                                                                                                                                   axis.text.x=element_text(angle=45, vjust=1,size=10, hjust=1.0, 
                                                                                                                                                                                                                                            colour="black" ))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=10),
        axis.text.y=element_text( vjust=.5,size=10, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Personality Traits and Cognitive Inhibition (Recent Probes)")+xlab("Trait")+ylab("Recent Probes Performance")

t<-ge+scale_y_discrete(limits=rev((levels=c("FamMinHiFamAcc",
                                           "FamMinFamRIAcc", 
                                           "HiFamMinFamRIAcc", 
                                           "FamMinHiFamRT", 
                                           "FamMinFamRIRT", "HiFamMinFamRIRT"))),
                      labels=rev(c("Fam-HiFam ACC",
                                   "Fam-FamRI ACC",
                                   "HiFam-FamRI ACC",
                                   "Fam-HiFam RT",
                                   "Fam-FamRI RT",
                                   "HiFam-FamRI RT")))
  l<-t+ scale_x_discrete(labels=c("ATQ Act Control",
                                  "ATQ Attention Control", "ATQ Att Shifting",
                                  "ATQ Effortful Control", "ATQ Inhibitory Control",
                                  "MPQ Absorption", "MPQ Achievement",
                                  "MPQ Aggression", "MPQ Alienation",
                                  "MPQ Control", "MPQ Constraint", 
                                  "MPQ Harm Avoidance",
                                  "MPQ Neg Emotionality", "MPQ NegEmot-Agentic",
                                  "NegEmot-Alienated", "PosEmot-Communal","MPQ Pos Emotionality",
                                  "PosEmot-Agentic",
                                  "MPQ Social Closeness", "MPQ Social Potency",
                                  "MPQ Stress Reaction", "MPQ Traditionalism",
                                  "MPQ Unlikely Virtues", "MPQ Well Being",
                                  "SNAP Aggression", "SNAP Back Dev", "SNAP Dependency",
                                  "SNAP Detachment", "SNAP Deviance", "SNAP Disinhibition:",
                                  "SNAP Eccentric Perceptions", "SNAP Enitlement", "SNAP Exhibitionism",
                                  "SNAP Impulsivity",
                                  "SNAP Low Self-Esteem", "SNAP Manipulativeness",
                                  "SNAP Mistrust", "SNAP Neg Temp", "SNAP Pos Temp",
                                  "SNAP Propriety", "SNAP Rare Virtues",
                                  "SNAP Self-Harm", "SNAP Suicide Proneness",
                                  "SNAP Workaholism", "K10", "STAI")) 
  
  l+ theme(plot.title=element_text(lineheight=.8,face="bold", size=20))+coord_flip()




t<-ge+scale_y_discrete(limits=rev((levels=c("RP_Nfam1_ACC",
                                         "RP_Nfam2_ACC", 
                                         "RP_Nfam1RI_ACC", 
                                         "RP_Nfam1_RT", 
                                         "RP_Nfam2_RT", "RP_Nfam1RI_RT"))),
                    labels=rev(c("Neg. Fam. ACC",
                                 "Neg. High Fam. ACC",
                                 "Neg. Resp. Int. ACC",
                                 "Neg Fam. RT",
                                 "Neg High Fam. RT",
                                 "Neg Resp. Int. RT")))+
 l<-t+ scale_x_discrete(labels=c("ATQ Act Control",
                            "ATQ Attention Control", "ATQ Att Shifting",
                            "ATQ Effortful Control", "ATQ Inhibitory Control",
                            "MPQ Absorption", "MPQ Achievement",
                            "MPQ Aggression", "MPQ Alienation",
                            "MPQ Control", "MPQ Constraint", 
                            "MPQ Harm Avoidance",
                            "MPQ Neg Emotionality", "MPQ NegEmot-Agentic",
                            "NegEmot-Alienated", "PosEmot-Communal","MPQ Pos Emotionality",
                            "PosEmot-Agentic",
                            "MPQ Social Closeness", "MPQ Social Potency",
                            "MPQ Stress Reaction", "MPQ Traditionalism",
                            "MPQ Unlikely Virtues", "MPQ Well Being",
                            "SNAP Aggression", "SNAP Back Dev", "SNAP Dependency",
                            "SNAP Detachment", "SNAP Deviance", "SNAP Disinhibition:",
                            "SNAP Eccentric Perceptions", "SNAP Enitlement", "SNAP Exhibitionism",
                            "SNAP Impulsivity",
                            "SNAP Low Self-Esteem", "SNAP Manipulativeness",
                            "SNAP Mistrust", "SNAP Neg Temp", "SNAP Pos Temp",
                            "SNAP Propriety", "SNAP Rare Virtues",
                            "SNAP Self-Harm", "SNAP Suicide Proneness",
                            "SNAP Workaholism", "K10", "STAI")) +
  
 l+ theme(plot.title=element_text(lineheight=.8,face="bold", size=20))+coord_flip()

#factor Analysis
library(nFactor)
ev<-eigen(cor(cmat, use="pairwise.complete.obs"))
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
forP.pca<-princomp(na.omit(forP))
summary(forP.pca)
plot(forP.pca)



SAS Notes:Trait
  
# Start by subset out where trial does not include positive. Only want negative probe trials
#*note that only Y trials have 36 trials (the rest have 9)
#subject not in (40, 65, 100))); *near 50% in terms of accuracy;
#  rename Y = RP_Y_RT
#Nfam0 = RP_Nfam0_RT
#Nfam1 = RP_Nfam1_RT
#Nfam1RI = RP_Nfam1RI_RT
#Nfam2 = RP_Nfam2_RT;

#RP_N_RT = mean (of Nfam0 Nfam1 Nfam1RI Nfam2); *this is a mean of means, not agg of the trials;
#RP_familiar_RT = mean (of Nfam1 Nfam2); *Nfam1RI this is a mean of means, not agg of the trials;
#RP_RT = mean (of RP_N_RT Y);

#*CONSIDER SHIFTING SCORING TO SUM, NOT MEAN;
  #rename Y = RP_Y_ACC
#Nfam0 = RP_Nfam0_ACC
#Nfam1 = RP_Nfam1_ACC
#Nfam1RI = RP_Nfam1RI_ACC
#Nfam2 = RP_Nfam2_ACC;
#RP_N_ACC = mean(of Nfam0 Nfam1 Nfam1RI Nfam2);
#RP_ACC = mean(of Y RP_N_ACC

#*score error sums;
ods listing close;
proc tabulate data = diss.RecentProbes;
class subject condition;
var ProbeDisplay_ACC;
tables subject * ProbeDisplay_ACC * condition, sum;
ods output table=sumAccurate;

#ProbeDisplay_ACC_SUM by condition. add this column