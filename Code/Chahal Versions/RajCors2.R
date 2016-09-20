
library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
source('/home/rajchahal/Desktop/WPIC/flanker_Control/corwithtarget.R')

#This dataset has added flanker control acc means column
alldat<-read.table('/home/rajchahal/Desktop/WPIC/alldat.csv',sep=",",header=T)
alldat<-within(alldat,{Subject<-factor(Subject)})
alldat$FlankAcc <- 1-(alldat$FlankAcc)
alldat$Flankpinac <- alldat$FlankAcc*100
alldat$Ang_I_PctInaccurate<-alldat$Ang_I_PctInaccurate*100
alldat$Sad_I_PctInaccurate<-alldat$Sad_I_PctInaccurate*100
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
flank <- c("Sad_I_PctInaccurate","Ang_I_PctInaccurate", "I_Sad","I_Ang", "Flankpinac", "I_MeanRT")

save(alldat,file='alldat.Rdata')
df<-alldat
load('alldat.Rdata')
cmat.list <- corwithtarget(df,target=flank,withvars=traits)
cmat.3dmat <- do.call(abind, list(cmat.list, along=0))
cmat.r <- melt( t(cmat.3dmat[,1,]),value.name='r' )
cmat.p <- melt( t(cmat.3dmat[,2,]),value.name='p' )
cmat <- merge(cmat.r,cmat.p, by=c("X1", "X2"))
names(cmat)[1:4] <- c('Trait','Flank', 'r', 'p')
ge<-ggplot( cmat, aes(x=Trait,y=Flank,fill=p,label=r) ) + theme_grey(base_size=22)+geom_tile() + scale_fill_continuous(high='#000000',low='#FF0000', limits=c(0,.10) ) + geom_text(color=I('white'), size=5)+theme(axis.title.x=element_text(face="bold",colour="Black", size=14),
      axis.text.x=element_text(angle=45, vjust=1,size=10, hjust=1.0, 
                               colour="black" ))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=10),
        axis.text.y=element_text( vjust=.5,size=10, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Personality Traits and Interference Control Correlations")+xlab("Trait")+ylab("Interference Control Performance")
ge+scale_y_discrete(limits=rev((levels=c("Flankpinac",
                                         "I_MeanRT", 
                                         "Sad_I_PctInaccurate", 
                                         "Ang_I_PctInaccurate", 
                                         "I_Sad", "I_Ang"))),
                    labels=rev(c("Flank Errors",
                             "Mean Flank RT",
                             "Mean E-Flank Sad Errors",
                             "Mean E-Flank Angry Errors",
                             "Mean E-Flank Sad RT",
                             "Mean E-Flank Angry RT")))+
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
                            "SNAP Workaholism", "K10", "STAI")) +
  
  theme(plot.title=element_text(lineheight=.8,face="bold", size=20))+coord_flip()
   

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
