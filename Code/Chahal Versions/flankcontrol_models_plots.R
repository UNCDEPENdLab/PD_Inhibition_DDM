using alldata frame
flanker
flank<-read.table('/home/rajchahal/Desktop/WPIC/flanker_Control/flanker_control.csv',sep=",",header=T)
flank<-within(flank,{Subject<-factor(Subject)})
flank<-within(flank,{Target<-factor(Target, labels=c("Congruent", "Incongruent"))})
flank<-within(flank,{Block<-factor(Block, labels=c("Mostly Congruent", "Mostly Incongruent"))})

#RT Observed
g<-ggplot(conflank, aes(x=Trial, y=TrialSlide_RT, color=Block))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2.5)+
  scale_color_brewer("Block", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text( vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Flanker Reaction Times")+xlab("\nTrial Congruency")+ylab("Mean Reaction Time (ms)\n")
g+ggtitle("Flanker Reaction Time")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))

#RT Model
fm1<-lmer(TrialSlide_RT~Trial*Block+(1|Subject), data=conflank)

#Doing omnibusss testss ussing car anova function on lmer get chi square vals for main effx
#Using new lmerCellMeans update:
cm<-lmerCellMeans(fm1)
d <- ggplot(cm, aes(x=Target, y=RT, 
                    color=Block, 
                    ymin=RT-se, ymax=RT+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Block", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(d)


#Errors observed
conflank$TrialSlide_ACC <- 1-(alldat$FlankAcc)
alldat$Flankpinac <- alldat$FlankAcc*100

g<-ggplot(conflank, aes(x=Trial, y=TrialSlide_ACC, color=Block))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Block", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90,vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Flanker Accuracy")+xlab("\nTrial Congruency")+ylab("Mean Accuracy (%)\n")
g+ggtitle("Flanker Accuracy")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))
#Err model 
fm1<-lmer(TrialSlide_ACC~Trial*Block+(1|Subject), data=conflank)
#Using new lmerCellMeans update:
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

#according to plots, the models are good fits



