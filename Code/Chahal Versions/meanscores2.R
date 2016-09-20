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

RTflank<-subset(eflank, FlankerDisplay_RT_Trim>0)
RTflank_noAS<-subset(RTflank, EmotionPairs!="Ang-Sad")
RTflank_noASI<-subset(RTflank_noAS, congruent=="incongruent")
eflank_noAS<-subset(eflank, EmotionPairs!="Ang-Sad")
eflank_noASI<-subset(eflank_noAS, congruent=="incongruent")

##################### Mean Scores Without Angry-Sad block ############################

#get means of accuracy for each subject on trials across emotions (exclude ang-sad block)

AccEmoPairs<-ddply(eflank, .(Subject, EmotionPairs), 
                   summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
  AccEmoPairs<-cast(AccEmoPairs, Subject~EmotionPairs, mean)
  colnames(AccEmoPairs) <- c("Subject","AccHapSad", "AccNeutSad","AccAngNeut", "AccAngHap", "AccHapNeut" )

AccEmo<-ddply(eflank, .(Subject, emotion), 
              summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
  AccEmo<-cast(AccEmo, Subject~emotion, mean)
  colnames(AccEmo) <- c("Subject","AccSad", "AccAng","AccHap", "AccNeut" )

AccCong<-ddply(eflank_noAS, .(Subject, congruent), 
               summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
  AccCong<-cast(AccCong, Subject~congruent, mean)
  colnames(AccCong) <- c("Subject","AccIncong", "AccCong")

#ACC by emotion pair and emotion on no ang-sad trials
AccEmoEmoPair<-ddply(eflank, .(Subject, EmotionPairs, emotion,congruent), summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccEmoEmoPair<-cast(AccEmoEmoPair, Subject~EmotionPairs+emotion, mean)
colnames(AccEmoEmoPair)<-c("Subject", "Acc_Sad_HapSad", "Acc_Hap_HapSad", "Acc_Sad_NeutSad", "Acc_Neut_NeutSad",
                          "Acc_Ang_AngNeut", "Acc_Neut_AngNeut", "Acc_Ang_AngHap",
                          "Acc_Hap_AngHap", "Acc_Happy_HapNeut", "Acc_Neut_HapNeut") 

#Acc by emotions on no ang-sad trials on only incongruent trials
AccEmoIncong<-ddply(eflank_noASI, .(Subject, emotion),
                    summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
  AccEmoIncong<-cast(AccEmoIncong, Subject~emotion, mean)
  colnames(AccEmoIncong) <- c("Subject","AccSad_I", "AccAng_I", "AccHappy_I", "AccNeut_I")

#Acc by emotion pair on no ang-sad trials, only incongruent trials
AccEmoPairIncong<-ddply(eflank_noASI, .(Subject, EmotionPairs), 
                     summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
  AccEmoPairIncong<-cast(AccEmoPairIncong, Subject~EmotionPairs, mean)
  colnames(AccEmoPairIncong)<-c("Subject", "AccHapSad_I", "AccNeutSad_I", "AccAngNeut_I", "AccAngHap_I", "AccHapNeut_I") 

#Acc by emotion pair and emotion on no ang-sad, only incongruent trials
AccEmoEmoPairI<-ddply(eflank_noASI, .(Subject, EmotionPairs, emotion), 
                     summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccEmoEmoPairI<-cast(AccEmoEmoPairI, Subject~EmotionPairs+emotion, mean)
colnames(AccEmoEmoPairI)<-c("Subject", "Acc_Sad_HapSad_I", "Acc_Hap_HapSad_I", "Acc_Sad_NeutSad_I", "Acc_Neut_NeutSad_I",
                           "Acc_Ang_AngNeut_I", "Acc_Neut_AngNeut_I", "Acc_Ang_AngHap_I",
                           "Acc_Hap_AngHap_I", "Acc_Happy_HapNeut_I", "Acc_Neut_HapNeut_I") 



#Get means for each subject on reaction time, excluding ang-sad block
RTEmoPairs<-ddply(RTflank_noAS, .(Subject, EmotionPairs), 
                   summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoPairs<-cast(RTEmoPairs, Subject~EmotionPairs, mean)
colnames(RTEmoPairs) <- c("Subject","RTHapSad", "RTNeutSad","RTAngNeut", "RTAngHap", "RTHapNeut" )

RTEmo<-ddply(RTflank_noAS, .(Subject, emotion), 
              summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
  RTEmo<-cast(RTEmo, Subject~emotion, mean)
  colnames(RTEmo) <- c("Subject","RTSad", "RTAng","RTHap", "RTNeut" )

RTCong<-ddply(RTflank_noAS, .(Subject, congruent), 
               summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
  RTCong<-cast(RTCong, Subject~congruent, mean)
  colnames(RTCong) <- c("Subject","RTIncong", "RTCong")

#RT by emotion pair and emotion on no ang-sad trials
RTEmoEmoPair<-ddply(RTflank_noAS, .(Subject, EmotionPairs, emotion), 
                          summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoEmoPair<-cast(RTEmoEmoPair, Subject~EmotionPairs+emotion, mean)
colnames(RTEmoEmoPair)<-c("Subject", "RT_Sad_HapSad", "RT_Hap_HapSad", "RT_Sad_NeutSad", "RT_Neut_NeutSad",
                                "RT_Ang_AngNeut", "RT_Neut_AngNeut", "RT_Ang_AngHap",
                                "RT_Hap_AngHap", "RT_Happy_HapNeut", "RT_Neut_HapNeut") 


#RT by emotions on no ang-sad trials on only incongruent trials
RTEmoIncong<-ddply(RTflank_noASI, .(Subject, emotion),
                    summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoIncong<-cast(RTEmoIncong, Subject~emotion, mean)
colnames(RTEmoIncong) <- c("Subject","RTSad_I", "RTAng_I", "RTHappy_I", "RTNeut_I")

#RT by emotion pair on no ang-sad trials, only incongruent trials
RTEmoPairIncong<-ddply(RTflank_noASI, .(Subject, EmotionPairs), 
                        summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoPairIncong<-cast(RTEmoPairIncong, Subject~EmotionPairs, mean)
colnames(RTEmoPairIncong)<-c("Subject", "RTHapSad_I", "RTNeutSad_I", "RTAngNeut_I", "RTAngHap_I", "RTHapNeut_I") 

#RT by emotion pair and emotion on no ang-sad trials, only incongruent
RTEmoEmoPairIncong<-ddply(RTflank_noASI, .(Subject, EmotionPairs, emotion), 
                       summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
RTEmoEmoPairIncong<-cast(RTEmoEmoPairIncong, Subject~EmotionPairs+emotion, mean)
colnames(RTEmoEmoPairIncong)<-c("Subject", "RT_Sad_HapSad_I", "RT_Hap_HapSad_I", "RT_Sad_NeutSad_I", "RT_Neut_NeutSad_I",
                                "RT_Ang_AngNeut_I", "RT_Neut_AngNeut_I", "RT_Ang_AngHap_I",
                                "RT_Hap_AngHap_I", "RT_Happy_HapNeut_I", "RT_Neut_HapNeut_I") 


######### MERGE ALL MEAN SCORES
AEflankAllMeans<-Reduce(function(...)merge(...,all=T, by.x="Subject"),
                        list(RTEmoPairIncong, RTEmoIncong,
                             RTCong, RTEmo, RTEmoPairs, RTEmoEmoPair, RTEmoEmoPairIncong,
                             AccEmoPairIncong,AccEmoIncong, AccEmoPairs,
                             AccEmo, AccCong, AccEmoEmoPair,
                             AccEmoEmoPairI))
#write.table(AEflankAllMeans, file="AllEmoFlankMeans.csv", col.names=T, sep=",")

#EflankAllMeans<-merge(RTEmoPairIncong, RTEmoIncong,  all=T, by.x="Subject")


######### Merge all mean scores with self report scores (not sure why sub #3 is missing)
selfreport<-read.table('/home/rajchahal/552014/May7/selfreports.csv',sep=",",header=T)
selfreport<-within(selfreport,{Subject<-factor(Subject)})
#drop same subs from self reports
seflreport<-subset(selfreport, Subject!="98" & Subject!="39" & Subject!="50" & Subject!="61" & Subject!="96")
#Merge by sub to mean scores table
EFlank_SelfReport<-Reduce(function(...)merge(...,all=T, by.x="Subject"),
                        list(selfreport, AEflankAllMeans))

write.table(EFlank_SelfReport, file="eflank_SelfReport.csv", col.names=T, sep=",")



################ Use EFlank_SelfReport for linear models #########
ex:
negemotionality~%inacc on ang trials, % inacc on sad trials + (1|Subject)
#(Should calculate % inaccurate first)
library(MASS)
#for poor accuracy outliers could use robust regression (e.g., rlm in MASS) 
gm1 <- lm(T_NegativeTemperament ~ AccSad_I, EFlank_SelfReport)

EFlank_SelfReport$pctinac <- EFlank_SelfReport$AccSad_I * 100

gm1 <- lm(T_NegativeTemperament ~ pctinac, EFlank_SelfReport)

summary(gm1)


#new lmer cell means model works for both directions of nesting: emotion/EmoPair and EmoPair/emotion
#then subset out the nonsense predictions
#ggplot(cm)
#NExt Steop: plot the lmer cell means (predicted) and plot the observed data( already done):
# See how well the predicted plots match the observed plots: they should be simialar for good models
