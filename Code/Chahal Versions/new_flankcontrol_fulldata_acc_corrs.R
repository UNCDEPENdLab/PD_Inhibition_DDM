flank<-read.sas7bdat('home/rajchahal/Desktop/WPIC/flanker_Control/flankerscores.sas7bdat',sep=",",header=T)
conflank<-read.table('/home/rajchahal/Desktop/WPIC/flanker_Control/flanker_control.csv',sep=",",header=T)
#TrialSlide_ACC
#TrialSlide_RT
#Incongruent: 0=cong, 1=incong
#CongruentBlock: 0=incong, 1=cong
conflank<-within(conflank,{Subject<-factor(Subject)})
conflank<-within(conflank,{Trial<-factor(Incongruent,levels=0:1, labels=c("congruent", "incongruent"))})
conflank<-within(conflank,{Block<-factor(CongruentBlock, levels=c("1", "0"), labels=c("Mostly Congruent", "Mostly Incongruent"))})


#Replace outlier values of subject trial means with NA in new col- FlankerDisplay_RT_Trim
eflank <- ddply(eflank, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$FlankerDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$FlankerDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$FlankerDisplay_RT_Trim <- subdf$FlankerDisplay_RT
  subdf$FlankerDisplay_RT_Trim[which(subdf$FlankerDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})

cm<-lmerCellMeans(fm1)
d <- ggplot(cm, aes(x=Target, y=Error, 
                    color=Block, 
                    ymin=Error-se, ymax=Error+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Block", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(d)