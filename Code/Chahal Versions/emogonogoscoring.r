#Emo Go No Go
#Read in data
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
#GO no GO
egg<-read.sas7bdat("C:/Users/chahalr/Desktop/sas/gonogo.sas7bdat")
egg<-read.sas7bdat("/home/rajchahal/emogonogo/emotgonogo.sas7bdat")

egg<-within(egg, {Subject<-factor(Subject)})
egg<-within(egg, {Procedure_Trial_<-factor(Procedure_Trial_,
                                          levels=c("OneProc",
                                                   "ThreeProc",
                                                   "FiveProc"),
                                          labels=c("OneGo",
                                                   "ThreeGo",
                                                   "FiveGo"))})
egg<-within(egg, {TargetEmotion<-factor(TargetEmotion,
                                    levels=c("Happy", "Fear", "Anger","Neutral"),
                                    labels=c("Hap", "Fear", "Ang", "Neut"))})

egg<-within(egg, {EmotionPairs<-factor(EmotionPairs,
levels=c(1:8), labels=c("FearNeut", "NeutAng", "HapFear", "FearAng", "AngNeut", "FearHap", "AngFear", "NeutFear"))})
#Combine the emo pairs
#egg<-within(egg, {EmotionPairs<-factor(EmotionPairs,
 #                                      levels=c(1:8), labels=c("FearNeut", "AngNeut", "HapFear", "AngFear", "AngNeut", "HapFear", "AngFear", "FearNeut"))})

#recategorize vars
egg$EP2[egg$EmotionPairs%in%c('FearNeut', 'NeutFear')]<-1
egg$EP2[egg$EmotionPairs%in%c('AngNeut', 'NeutAng')]<-2
egg$EP2[egg$EmotionPairs%in%c('FearAng', 'AngFear')]<-3
egg$EP2[egg$EmotionPairs%in%c('HapFear', 'FearHap')]<-4

egg<-within(egg,{EP2<-factor(EP2, levels=c(1:4), labels=c("Fear-Neut", "Ang-Neut", "Fear-Ang", "Hap-Fear"))})

#Label where its go or no go
egg$GoType<-as.numeric(egg$Procedure_SubTrial=="NoGoProc")
egg<-within(egg,{GoType<-factor(GoType, levels=c(0:1), labels=c("Go", "NoGo"))})

#If GoDisplay_RT=0, no response. FailedResp=1, else FailedResp=0
#egg$failedResp <- as.numeric(gonogo$GoDisplay_RT==0)
#If NoGoDisplay_RESP="{space}" then falsealarm=1, else falsealarm=0
egg$falseAlarm<-as.numeric(egg$NoGoSlide_RESP=="{SPACE}")

#RT's: subset to each block type, then drop per cat on IQR's

#subset block type = OneGoTrial, then if GoDisplay_RT>3*IQR then RT=NA
UnoGo<-subset(egg, BlockType=="OneGo")
#Replace outlier values with NA
OneGo_RT_Trim<- ddply(UnoGo, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$GoSlide_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$GoSlide_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$GoSlide_RT_Trim <- subdf$GoSlide_RT
  subdf$GoSlide_RT_Trim[which(subdf$GoSlide_RT_Trim > upperCut)] <- NA
  return(subdf)
})
egg$GoSlide_RT_Trim<-1
egg[egg$BlockType=="OneGo",]<-OneGo_RT_Trim



#subset block type = ThreeGoTrial, then if GoDisplay_RT>3*IQR then RT=NA
ThreeGo<-subset(egg, BlockType=="ThreeGo")
#Replace outlier values with NA
ThreeGo_RT_Trim<- ddply(ThreeGo, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$GoDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$GoDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$GoDisplay_RT_Trim <- subdf$GoDisplay_RT
  subdf$GoDisplay_RT_Trim[which(subdf$GoDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})

gonogo[gonogo$BlockType=="ThreeGo",]<-ThreeGo_RT_Trim

#Subset blcok type=Fivego, then RT=NA if outlier
FiveGo<-subset(egg, BlockType=="FiveGo")
#Replace outlier values with NA
FiveGo_RT_Trim<- ddply(FiveGo, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$GoDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$GoDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$GoDisplay_RT_Trim <- subdf$GoDisplay_RT
  subdf$GoDisplay_RT_Trim[which(subdf$GoDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})

gonogo[gonogo$BlockType=="FiveGo",]<-FiveGo_RT_Trim


#Now that all outlier RTs have been replaced by NA...
#PLot observed values

#Rt Observed
g<-ggplot(egg, aes(x=TargetEmotion, y=GoSlide_RT, color=EP2, shape=Procedure_Trial_), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)

library(ggplot2)
g<-ggplot(egg, aes(x=TargetEmotion, y=GoSlide_RT, color=EP2, 
                      shape=Procedure_Trial_))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emotional GNG RT")+xlab("Emotion\n")+ylab("Mean RTs")
g+ggtitle("Emotional GNG RTs")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("BlockType")

#Save as 8x11 landscape

#In False alarms, it always no-Go so dont include block type and subset to only no-go trials
Stop<-subset(egg, GoType=="NoGo")

dd<-ggplot(egg, aes(x=TargetEmotion, y=falseAlarm, color=EP2, shape=Procedure_Trial_))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emotional GNG FA")+xlab("Emotion\n")+ylab("Mean FAs")
dd+ggtitle("Emotional GNG FAs")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("BlockType")




#Note on  models: will need to have two. one without hap fear and one with only hap fear
#emotion is a nested var of emotion pairs so change in model
#RT Model 1 without hap fear
egg_NoHapFear<-subset(egg, EP2!="Hap-Fear")
#here, when tried to include GoType- said error: require 2 or more levels
fm1<-lmer(GoSlide_RT~TargetEmotion*TargetEmotion/EP2*Procedure_Trial_+(1|Subject), data=egg_NoHapFear)
cm<-lmerCellMeans(fm1)
#Subset out nonsense vars in cm
#TO SUBSET OUT NONSENSE VARS, merge emo and emo pair and remove nonsense
cm$EmoEmoPair <- paste(cm[,1], cm[,3], sep=" ")
cm<-subset(cm, !(EmoEmoPair %in% c("Ang Fear-Neut", "Fear Ang-Neut", "Neut Fear-Ang")))

#PLOT RT Model
RTMod <- ggplot(cm, aes(x=TargetEmotion, y=GoSlide_RT, 
                    color=EP2, shape=Procedure_Trial_,
                    ymin=GoSlide_RT-se, ymax=GoSlide_RT+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
RTMod+ggtitle("Emotional GNG Model RT")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("BlockType")

#FA Model
#subset to NoGo Trials only since that when false alarm occurs
#egg_NoHapFear_nogo<-subset(egg_NoHapFear, GoType=="NoGo")
#block type should be included but errorss

mod <- gdata::drop.levels(egg_NoHapFear)
mod <- gdata::drop.levels(subset(mod, GoType=="NoGo"))

gm1 <- glmer(falseAlarm ~ TargetEmotion*TargetEmotion/EP2*Procedure_Trial_+ (1|Subject), data=mod, family="binomial")
m<-lmerCellMeans(gm1)

table(~BlockType + TargetEmotion + EP2, egg_NoHapFear)


m$EmoEmoPair <- factor(paste(m[,1], m[,2], sep=" "))
#Subset out nonsense vars
m<-subset(m, !(EmoEmoPair %in% c("Ang Fear-Neut", "Fear Ang-Neut", "Neut Fear-Ang")))
#Plot FA Model
#Aggregate across trials
m_aggTrials <- ddply(m, .(TargetEmotion, EP2), function(subdf) {
  subdf <- as.data.frame(lapply(subdf, function(col) {
    if (is.factor(col)) {
      head(col, n=1)
    } else {
      mean(col, na.rm=TRUE)
    }
  }))
  
  subdf$SubTrial <- NULL
  subdf
})
#PLOT Using AGG Scores from model
dh1_aggTrials <- ggplot(m_aggTrials, aes(x=TargetEmotion, y=falseAlarm, 
                                        color=EP2, 
                                        ymin=falseAlarm - se, ymax=falseAlarm + se)) +
  geom_pointrange(size=1, position=position_dodge(width=0.5))
#Convert form log to % false alarms
m_aggTrials$acchat <- with(m_aggTrials, { exp(falseAlarm)/(1+exp(falseAlarm))})
m_aggTrials$acchi <- with(m_aggTrials, { exp((falseAlarm+se))/(1+exp(falseAlarm+se))})
m_aggTrials$acclo <- with(m_aggTrials, { exp((falseAlarm-se))/(1+exp(falseAlarm-se))})

#Plot FA Model on percentage scale.
dh_aggTrials <- ggplot(m_aggTrials, aes(x=TargetEmotion, y=acchat, 
                                        color=EP2,
                                        ymin=acclo, ymax=acchi)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2") +
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(dh_aggTrials)


#then redo FA and RT models with including happy trials

#RT Model 2 with  only hap fear
egg_HapFear<-subset(egg, EP2=="Hap-Fear")
#here, when tried to include GoType- said error: require 2 or more levels

#Leave out Emotion Pairs in Model when its justs hap-fear
fm2<-lmer(GoSlide_RT~TargetEmotion*BlockType+(1|Subject), data=egg_HapFear)
cm<-lmerCellMeans(fm2)

#SNo need to subset out nonsense vars since EP not included in model

#PLOT RT Model
RTMod2 <- ggplot(cm, aes(x=TargetEmotion, y=GoSlide_RT, 
                        color=BlockType,
                        ymin=GoSlide_RT-se, ymax=GoSlide_RT+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())

#FA Model
#subset to NoGo Trials only since that when false alarm occurs
#egg_HapFear_nogo<-subset(egg_HapFear, GoType=="NoGo")
#*******block type not included but SHOULD
gm2 <- glmer(falseAlarm ~ TargetEmotion+ (1|Subject), data=egg_HapFear, family="binomial")
m<-lmerCellMeans(gm2)


#Plot FA Model
#Aggregate across trials
m_aggTrials <- ddply(m, .(TargetEmotion), function(subdf) {
  subdf <- as.data.frame(lapply(subdf, function(col) {
    if (is.factor(col)) {
      head(col, n=1)
    } else {
      mean(col, na.rm=TRUE)
    }
  }))
  
  subdf$SubTrial <- NULL
  subdf
})
#PLOT Using AGG Scores from model
dh2_aggTrials <- ggplot(m_aggTrials, aes(x=TargetEmotion, y=falseAlarm, 
                                         ymin=falseAlarm - se, ymax=falseAlarm + se)) +
  geom_pointrange(size=1, position=position_dodge(width=0.5))
#Convert form log to % false alarms
m_aggTrials$acchat <- with(m_aggTrials, { exp(falseAlarm)/(1+exp(falseAlarm))})
m_aggTrials$acchi <- with(m_aggTrials, { exp((falseAlarm+se))/(1+exp(falseAlarm+se))})
m_aggTrials$acclo <- with(m_aggTrials, { exp((falseAlarm-se))/(1+exp(falseAlarm-se))})

#Plot FA Model on percentage scale.
dh_aggTrials <- ggplot(m_aggTrials, aes(x=TargetEmotion, y=acchat, 
                                        ymin=acclo, ymax=acchi)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(dh_aggTrials)

#Now weve plotted both types of models on RT and FA, time to correlate to self report measures using heat map
#Most difficult to do was fear-ang emotion pair
#For DVS: Ang FA, Ang RT, Fear FA, fear RT

#For Mean RTS
#Group Descriptives
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

eggfear<-subset(egg,TargetEmotion=="Fear")
desc<-getDescriptivesByGroup(eggfear, c("GoSlide_RT"), "Subject")

alldat$FearRT<-desc$mean

#########################################################

library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
library(reshape)
library(ggplot2)
source('/home/rajchahal/Desktop/WPIC/flanker_Control/corwithtarget.R')

#read in all dat
alldat<-read.table('/home/rajchahal/Desktop/WPIC/alldat.csv',sep=",",header=T)
alldat<-within(alldat,{Subject<-factor(Subject)})

alldat<-subset(alldat, Subject!="99" & Subject!="7")

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

eggo <- c("fearang_FAs","angfear_FAs")

save(alldat,file='alldat.Rdata')
df<-alldat
load('alldat.Rdata')
cmat.list <- corwithtarget(df,target=eggo,withvars=traits)
cmat.3dmat <- do.call(abind, list(cmat.list, along=0))
cmat.r <- melt( t(cmat.3dmat[,1,]),value.name='r' )
cmat.p <- melt( t(cmat.3dmat[,2,]),value.name='p' )
cmat <- merge(cmat.r,cmat.p, by=c("X1", "X2"))
names(cmat)[1:4] <- c('traits','eggo', 'r', 'p')
ge<-ggplot( cmat, aes(x=traits,y=eggo,fill=p,label=r) ) + geom_tile() + scale_fill_continuous(high='#000000',low='#FF0000', limits=c(0,.10) ) + geom_text(color=I('white'),size=5)+theme(axis.title.x=element_text(face="bold",colour="Black", size=14),
                                                                                                                                                                                         axis.text.x=element_text(angle=90, vjust=.5,size=14, 
                                                                                                                                                                                                                  colour="black", ))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=10),
        axis.text.y=element_text( vjust=.5,size=10, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Personality Traits and Emo Go/No-Go Correlations")+xlab("Trait")+ylab("Response Inhbition Performance")
d<-ge+scale_y_discrete(limits=rev((levels=c("fearang_FAs",
                                            "angfear_FAs"))),
                       labels=rev(c("Fear(-Ang) FAs",
                                    "Ang(-Fear) FAs")))+
  scale_x_discrete(labels=c("ATQ Act Control",
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
d+theme(plot.title=element_text(lineheight=.8,face="bold", size=20))+coord_flip()

#Need to look at TRIM and VRIM sscores to know who to drop based on how they answred the questionnairs in final dataset

splom(cmat)
###

#  brew <- brewer.pal(name="RdBu",n=5)
#  #get color ramp
#  cc.brew <- colorRampPalette(brew)
#  #apply color ramp
#  cc <- rev(cc.brew(nrow(cmat.mat)))
#  
#  
#  #pdf("heatNPD.pdf", width=14, height=16)
#  heatmap.2(cmat.mat,symm=TRUE, col=cc, main="", Rowv=TRUE, Colv=TRUE, trace="none", density.info="none")
#  #title(main="Correlation Table (Ordered by Dendrogram)",font.main=1,outer=TRUE,line=-1)#,cex.main=2)
#  #dev.off()
#  



# 
# cora<-cor(alldate$T_NegativeTemperament, alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$T_Mistrust,            alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$T_Manipulativeness,    alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_NET,               alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_alT,               alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_agT,               alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_srT,               alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$ATQ_ECTotal,           alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$ATQ_InhibControl,      alldate$Sad_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$ATQ_AttControl,        alldate$Sad_I_PctInaccurate, use="complete.obs")
# 
# cora<-cor(alldate$T_Mistrust,            alldate$Ang_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_NET,               alldate$Ang_I_PctInaccurate, use="complete.obs")
# cora<-cor(alldate$MPS_agT,               alldate$Ang_I_PctInaccurate, use="complete.obs")
# 
# cora<-cor(alldate$T_Entitlement,         alldate$I_Sad,               use="complete.obs")
# cora<-cor(alldate$T_Mistrust,            alldate$I_Sad,               use="complete.obs")
# 
# cora<-cor(alldate$T_Entitlement,         alldate$I_Ang,               use="complete.obs")

