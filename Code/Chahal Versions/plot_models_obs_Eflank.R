library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(reshape)
library(multcomp)
library(hmisc)

eflank<-read.table('/home/rajchahal/552014/sasdata_emotflanker.csv',sep=",",header=T)

eflank<-within(eflank,{Subject<-factor(Subject)})
eflank<-within(eflank,{emotion<-factor(emotion,levels=1:4, labels=c("Sad", "Ang", "Happy", "Neut"))})
eflank<-within(eflank,{EmotionPairs<-factor(EmotionPairs, levels=1:6, labels=c("Ang-Sad", "Hap-Sad", "Neut-Sad", "Ang-Neut", "Ang-Hap", "Hap-Neut"))})
eflank<-within(eflank,{congruent<-factor(congruent,levels=0:1, labels=c("incongruent", "congruent"))})


#Replace outlier values of subject trial means with NA in new col- FlankerDisplay_RT_Trim
eflank <- ddply(eflank, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$FlankerDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$FlankerDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$FlankerDisplay_RT_Trim <- subdf$FlankerDisplay_RT
  subdf$FlankerDisplay_RT_Trim[which(subdf$FlankerDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})

#Dropped Subjects: 96, 39, 50, 61, 98
eflank<-subset(eflank, Subject!="98" & Subject!="39" & Subject!="50" & Subject!="61" & Subject!="96")

RTflank<-subset(eflank, FlankerDisplay_RT_Trim>0)
RTflank_noAS<-subset(RTflank, EmotionPairs!="Ang-Sad")
RTflank_noASI<-subset(RTflank_noAS, congruent=="incongruent")
eflank_noAS<-subset(eflank, EmotionPairs!="Ang-Sad")
eflank_noASI<-subset(eflank_noAS, congruent=="incongruent")

#RT Lmer Model
fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank)
#Using new lmerCellMeans update:
cm<-lmerCellMeans(fm1)
#TO SUBSET OUT NONSENSE VARS, merge emo and emo pair and remove nonsense
cm$EmoEmoPair <- paste(cm[,1], cm[,4], sep=" ")
cm<-subset(cm, !(EmoEmoPair %in% c("Happy Ang-Sad", "Neut Ang-Sad", "Neut Hap-Sad",
                                   "Ang Hap-Sad", "Neut Hap-Sad", "Ang Neut-Sad",
                                   "Happy Neut-Sad", "Sad Ang-Neut", "Happy Ang-Neut",
                                   "Sad Ang-Hap", "Neut Ang-Hap", "Sad Hap-Neut", 
                                   "Ang Hap-Neut")))

m_aggTrials <- ddply(cm, .(emotion, congruent, EmotionPairs), function(subdf) {
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

#Plot RT Model
d <- ggplot(cm, aes(x=emotion, y=FlankerDisplay_RT_Trim, 
                                        color=EmotionPairs, shape=congruent,
                                        ymin=FlankerDisplay_RT_Trim-se, ymax=FlankerDispslay_RT_Trim+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(dh_aggTrials)








#ACC GLM Model
gm1 <- glmer(FlankerDisplay_ACC ~ emotion*emotion/EmotionPairs*congruent + SubTrial  + (1+SubTrial | Subject), eflank, family="binomial")
m<-lmerCellMeans(gm1)
m$EmoEmoPair <- factor(paste(m[,1], m[,4], sep=" "))
#Subset out nonsense vars
m<-subset(m, !(EmoEmoPair %in% c("Happy Ang-Sad", "Neut Ang-Sad", "Neut Hap-Sad",
                                 "Ang Hap-Sad", "Neut Hap-Sad", "Ang Neut-Sad",
                                 "Happy Neut-Sad", "Sad Ang-Neut", "Happy Ang-Neut",
                                 "Sad Ang-Hap", "Neut Ang-Hap", "Sad Hap-Neut", 
                                 "Ang Hap-Neut")))

#Plot Acc Model
#Aggregate across trials
m_aggTrials <- ddply(m, .(emotion, congruent, EmotionPairs), function(subdf) {
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

#Plot using aggregate scores from model
dh_aggTrials <- ggplot(m_aggTrials, aes(x=emotion, y=FlankerDisplay_ACC, 
                                        color=EmotionPairs, shape=congruent,
                                        ymin=FlankerDisplay_ACC - se, ymax=FlankerDisplay_ACC + se)) +
  geom_pointrange(size=1, position=position_dodge(width=0.5))

#Convert from log to percentage of correct
m_aggTrials$acchat <- with(m_aggTrials, { exp(FlankerDisplay_ACC)/(1+exp(FlankerDisplay_ACC))})
m_aggTrials$acchi <- with(m_aggTrials, { exp((FlankerDisplay_ACC+se))/(1+exp(FlankerDisplay_ACC+se))})
m_aggTrials$acclo <- with(m_aggTrials, { exp((FlankerDisplay_ACC-se))/(1+exp(FlankerDisplay_ACC-se))})

#Plot ACC Model on percentage scale.
dh_aggTrials <- ggplot(m_aggTrials, aes(x=emotion, y=acchat, 
                                        color=EmotionPairs, shape=congruent,
                                        ymin=acclo, ymax=acchi)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2") +ylim(0.85, .98)+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(dh_aggTrials)


# Observed Plots

#Acc Observed
g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_ACC, color=EmotionPairs, 
                      shape=congruent))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emotional Flanker Accuracy")+xlab("Emotion\n")+ylab("Mean Accuracy(%)\n")
g+ggtitle("Emotional Flanker Accuracy")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Congruence")


#RT Observed
eflank$RTse<-std.error(eflank$FlankerDisplay_RT_Trim)
g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_RT_Trim, color=EmotionPairs, 
                      shape=congruent))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(angle=90, vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emotional Flanker Reaction Time")+xlab("Emotion\n")+ylab("Mean Reaction Time (ms)\n")
g+ggtitle("Emotional Flanker Reaction Time")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))+scale_shape("Congruence")

