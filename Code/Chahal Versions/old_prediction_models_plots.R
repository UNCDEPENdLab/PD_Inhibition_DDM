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
summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+emotion/EmotionPairs*congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))

fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank)
#summary(fm1)
#conmat <- do.call(rbind, 
#                  getSimpleEffectsMatrix(
#                    names(fixef(fm1)), 
#                    iv.slice="emotion",
#                    iv.test="EmotionPairs", 
#                    levels.slice=levels(eflank$emotion),
#                    levels.test=levels(eflank$EmotionPairs)
#                  )
#)

#summary(glht(fm1, linfct=conmat))

cm<-lmerCellMeans(fm1)

#TO SUBSET OUT NONSENSE VARS, merge emo and emo pair and remove nonsense
cm$EmoEmoPair <- paste(cm[,1], cm[,4], sep=" ")

#cm<-subset(cm, EmoEmoPair!="Happy Ang-Sad")
#To subset out multiple Emo Pairs at one time, use this:
cm<-subset(cm, !(EmoEmoPair %in% c("Happy Ang-Sad", "Neut Ang-Sad", "Neut Hap-Sad",
                                   "Ang Hap-Sad", "Neut Hap-Sad", "Ang Neut-Sad",
                                   "Happy Neut-Sad", "Sad Ang-Neut", "Happy Ang-Neut",
                                   "Sad Ang-Hap", "Neut Ang-Hap", "Sad Hap-Neut", 
                                   "Ang Hap-Neut")))

#Now we can plot the model for RT on all trials
dh<-ggplot(cm, aes(x=emotion, y=FlankerDisplay_RT_Trim, 
                      color=EmotionPairs, shape=congruent))+ stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=1)
b<-h+geom_bar(stat="identity", position="dodge")+coord_cartesian(ylim=c(500,800))
d<-b+facet_grid(EmotionPairs~., scale="free_y")                           
z<-h+geom_point()


#Trying another model for RT
fm2<-lmer(FlankerDisplay_RT_Trim~emotion*EmotionPairs/emotion*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank)
#summary(fm1)
#conmat <- do.call(rbind, 
#                  getSimpleEffectsMatrix(
#                    names(fixef(fm1)), 
#                    iv.slice="emotion",
#                    iv.test="EmotionPairs", 
#                    levels.slice=levels(eflank$emotion),
#                    levels.test=levels(eflank$EmotionPairs)
#                  )
#)

#summary(glht(fm1, linfct=conmat))

cm2<-lmerCellMeans(fm2)

#TO SUBSET OUT NONSENSE VARS, merge emo and emo pair and remove nonsense
cm2$EmoEmoPair <- paste(cm2[,1], cm2[,2], sep=" ")

#cm<-subset(cm, EmoEmoPair!="Happy Ang-Sad")
#To subset out multiple Emo Pairs at one time, use this:
cm2<-subset(cm2, !(EmoEmoPair %in% c("Happy Ang-Sad", "Neut Ang-Sad", "Neut Hap-Sad",
                                   "Ang Hap-Sad", "Neut Hap-Sad", "Ang Neut-Sad",
                                   "Happy Neut-Sad", "Sad Ang-Neut", "Happy Ang-Neut",
                                   "Sad Ang-Hap", "Neut Ang-Hap", "Sad Hap-Neut", 
                                   "Ang Hap-Neut")))

#Need to subset cm before plotting. Subset out the parts that are nonsense predictions
#NOw we can Plot

#h<-ggplot(cm2, aes(x=emotion, y=FlankerDisplay_RT_Trim, 
#                  fill=EmotionPairs))+facet_wrap(~congruent) 
#b<-h+geom_bar(stat="identity", position="dodge")+coord_cartesian(ylim=c(500,780))
#d<-b+facet_grid(EmotionPairs~., scale="free_y")                           
#z<-h+geom_point()


##Original plot. not cm or fm1
#h<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_RT_Trim,
#                      fill=EmotionPairs))+stat_summary(
#                        fun.data="mean_cl_boot", 
#                        position=position_dodge(width=1.25), 
#                        geom="bar", width=0.5)+facet_wrap(~congruent)


#Or bar graph
#RTmeanseflank<-ddply(eflank, .(Subject, congruent, emotion, EmotionPairs), 
#                     summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))

#n<-ggplot(RTmeanseflank, aes(x=emotion, y=mean_RT, fill=EmotionPairs))+facet_wrap(~congruent)
#q<-n+geom_bar(stat="identity", position="dodge")+coord_cartesian(ylim=c(700,1350))



#Linear model for ACC on all trials
gm1 <- glmer(FlankerDisplay_ACC ~ emotion*emotion/EmotionPairs*congruent + SubTrial  + (1+SubTrial | Subject), eflank, family="binomial")
gm0 <- glmer(FlankerDisplay_ACC ~ emotion*congruent + SubTrial  + (1+SubTrial | Subject), eflank, family="binomial")
#For model-bulding. 
anova(gm0, gm1)


#POsner conflict score?

#gm1diag<-glm.diag(gm1)
#glm.diag.plots(gm1, gm1diag)

m<-lmerCellMeans(gm1)
m$EmoEmoPair <- factor(paste(m[,1], m[,4], sep=" "))
#Subset out nonsense vars
m<-subset(m, !(EmoEmoPair %in% c("Happy Ang-Sad", "Neut Ang-Sad", "Neut Hap-Sad",
                                   "Ang Hap-Sad", "Neut Hap-Sad", "Ang Neut-Sad",
                                   "Happy Neut-Sad", "Sad Ang-Neut", "Happy Ang-Neut",
                                   "Sad Ang-Hap", "Neut Ang-Hap", "Sad Hap-Neut", 
                                   "Ang Hap-Neut")))

#Plot Acc Model

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

dh_aggTrials <- ggplot(m_aggTrials, aes(x=emotion, y=FlankerDisplay_ACC, 
                               color=EmotionPairs, shape=congruent,
                               ymin=FlankerDisplay_ACC - se, ymax=FlankerDisplay_ACC + se)) +
  geom_pointrange(size=1, position=position_dodge(width=0.5))

m_aggTrials$acchat <- with(m_aggTrials, { exp(FlankerDisplay_ACC)/(1+exp(FlankerDisplay_ACC))})
m_aggTrials$acchi <- with(m_aggTrials, { exp((FlankerDisplay_ACC+se))/(1+exp(FlankerDisplay_ACC+se))})
m_aggTrials$acclo <- with(m_aggTrials, { exp((FlankerDisplay_ACC-se))/(1+exp(FlankerDisplay_ACC-se))})

dh_aggTrials <- ggplot(m_aggTrials, aes(x=emotion, y=acchat, 
                                        color=EmotionPairs, shape=congruent,
                                        ymin=acclo, ymax=acchi)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Emotion Pair", palette="Set2") +ylim(0.85, .98)+theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
       axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
       axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+theme(panel.grid.minor=element_blank())
print(dh_aggTrials)

#Observed ACC Plot
m_agg <- ddply(eflank, .(emotion, congruent, EmotionPairs), function(subdf) {
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
g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_ACC, color=EmotionPairs, 
                      shape=congruent))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=1.5)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=18),
 axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Flanker Accuracy")+xlab("Emotion")+ylab("Accuracy")
  


# Observed RT Plot
eflank$RTse<-std.error(eflank$FlankerDisplay_RT_Trim)
g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_RT_Trim, color=EmotionPairs, 
                      shape=congruent, ymin=FlankerDisplay_RT_Trim-RTse, 
                      ymax=FlankerDisplay_RT_Trim+RTse))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=1.5)+
  scale_color_brewer("Emotion Pair", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Emo Flanker Reaction Time")+xlab("Emotion")+ylab("Reaction Time (ms)")






#Observed Acc's plot
eflank$EmoEmoPair <- factor(paste(m[,1], m[,4], sep=" "))
eflank$se<-ddply(eflank, .(Subject, emotion, EmotionPairs, congruent), summarize, se=std.error(FlankerDisplay_ACC, na.rm=TRUE))
AccEmo<-cast(AccEmo, Subject~emotion, mean)
colnames(AccEmo) <- c("Subject","AccSad", "AccAng","AccHap", "AccNeut" )

g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_ACC, color=EmotionPairs, shape=congruent)+
                   stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=1)
accd <- ggplot(eflank, aes(x=emotion, y=FlankerDisplay_ACC, 
                                                  color=EmotionPairs, shape=congruent,
                                                  ymin=FlankerDisplay_ACC - std.error(FlankerDisplay_ACC), ymax=FlankerDisplay_ACC + std.error(FlankerDisplay_ACC)) +geom_pointrange(size=1, position=position_dodge(width=0.5))
          
#g<-ggplot(eflank, aes(x=emotion, y=FlankerDisplay_ACC, fill=EmotionPairs))+
#  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=.8), geom="bar", width=.5)+facet_wrap(~congruent)

