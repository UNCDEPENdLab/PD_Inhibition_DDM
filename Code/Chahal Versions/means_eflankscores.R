library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(reshape)
library(multcomp)

eflank<-read.table('C:/Users/chahalr/Desktop/552014/sasdata_emotflanker.csv',sep=",",header=T)
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
#Mean scores on acc 
AccScore<-ddply(eflank, .(Subject, EmotionPairs, emotion,congruent), 
                     summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccEmoEmoPair<-cast(AccEmoEmoPair, Subject~EmotionPairs+emotion+congruent, mean)
#colnames(AccEmoEmoPair)<-c("Subject", "Acc_Sad_HapSad", "Acc_Hap_HapSad", "Acc_Sad_NeutSad", "Acc_Neut_NeutSad",
#                           "Acc_Ang_AngNeut", "Acc_Neut_AngNeut", "Acc_Ang_AngHap",
#                           "Acc_Hap_AngHap", "Acc_Happy_HapNeut", "Acc_Neut_HapNeut")

#Mean scores RT
RTScore<-ddply(eflank, .(Subject, EmotionPairs, emotion,congruent), 
                    summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoEmoPair<-cast(RTEmoEmoPair, Subject~EmotionPairs+emotion+congruent, mean)
#colnames(AccEmoEmoPair)<-c("Subject", "Acc_Sad_HapSad", "Acc_Hap_HapSad", "Acc_Sad_NeutSad", "Acc_Neut_NeutSad",
#                           "Acc_Ang_AngNeut", "Acc_Neut_AngNeut", "Acc_Ang_AngHap",
#                           "Acc_Hap_AngHap", "Acc_Happy_HapNeut", "Acc_Neut_HapNeut")

#Merge Data Frames
AccScore$mean_RT<-(RTScore$mean_RT)

#Merge two data frames of means for acc and RT
#escore<-Reduce(function(...)merge(...,all=T, by="Subject"),
#                        list(RTEmoEmoPair, AccEmoEmoPair))
#escore2<-merge(AccScore, RTScore,  all=T, by="Subject")


#Convert to PCt  Inaccurate
library(MASS)
#for poor accuracy outliers could use robust regression (e.g., rlm in MASS) 
gm1 <- lm(T_NegativeTemperament ~ AccSad_I, EFlank_SelfReport)

AccScore$pctinac <- 1-(AccScore$mean_ACC)
AccScore$pctinac <- AccScore$pctinac*100

#Merge with self reports
alldat<-read.table('/home/rajchahal/552014/sasdata_allcombined.csv',sep=",",header=T)
escore<-merge(AccScore, alldat,  all=T, by="Subject")


gm1 <- lm(T_NegativeTemperament ~ pctinac, EFlank_SelfReport)

summary(gm1)