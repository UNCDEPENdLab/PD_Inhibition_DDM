library(sas7bdat)
library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(Hmisc)
library(reshape)
#GO no GO
gonogo<-read.sas7bdat("C:/Users/chahalr/Desktop/sas/gonogo.sas7bdat")
gonogo<-read.sas7bdat("/home/rajchahal/Desktop/gonogo/gonogo.sas7bdat")

write.table(gonogo, "gonogo.csv", col.names=T, sep=",")
#head(gonogo)
gonogo<-read.table("C:/Users/chahalr/Documents/gonogo.csv", header=T, sep=",")
#These next two steps are not needed if using the CSV. Only needed for SAS
gonogo<-within(gonogo, {Subject<-factor(Subject)})
gonogo<-within(gonogo, {BlockType<-factor(BlockType,
                                          levels=c("OneGoTrial",
                                                   "ThreeGoTrial",
                                                   "FiveGoTrial",
                                                   "SevenGoTrial"),
                                          labels=c("OneGo",
                                                   "ThreeGo",
                                                   "FiveGo",
                                                   "SevenGo"))})

#If GoDisplay_RT=0, no response. FailedResp=1, else FailedResp=0
gonogo$failedResp <- as.numeric(gonogo$GoDisplay_RT==0)
#If NoGoDisplay_RESP="{space}" then falsealarm=1, else falsealarm=0
gonogo$falseAlarm<-as.numeric(gonogo$NoGoDisplay_RESP=="{SPACE}")
#RT's: subset to each block type, then drop per cat on IQR's

#subset block type = OneGoTrial, then if GoDisplay_RT>3*IQR then RT=NA
OneGo<-subset(gonogo, BlockType=="OneGo")
#Replace outlier values with NA
OneGo_RT_Trim<- ddply(OneGo, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$GoDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$GoDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$GoDisplay_RT_Trim <- subdf$GoDisplay_RT
  subdf$GoDisplay_RT_Trim[which(subdf$GoDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})
gonogo$GoDisplay_RT_Trim<-1
gonogo[gonogo$BlockType=="OneGo",]<-OneGo_RT_Trim


#subset block type = ThreeGoTrial, then if GoDisplay_RT>3*IQR then RT=NA
ThreeGo<-subset(gonogo, BlockType=="ThreeGo")
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
FiveGo<-subset(gonogo, BlockType=="FiveGo")
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


#Subset to block type= seven go, then drop outliers
SevenGo<-subset(gonogo, BlockType=="SevenGo")
#Replace outlier values with NA
SevenGo_RT_Trim<- ddply(SevenGo, .(Subject), function(subdf) {
  pct75 <- quantile(subdf$GoDisplay_RT, probs=0.75, na.rm=T)
  iqr <- IQR(subdf$GoDisplay_RT, na.rm=T)
  upperCut <- pct75 + 3*iqr
  subdf$GoDisplay_RT_Trim <- subdf$GoDisplay_RT
  subdf$GoDisplay_RT_Trim[which(subdf$GoDisplay_RT_Trim > upperCut)] <- NA
  return(subdf)
})

gonogo[gonogo$BlockType=="SevenGo",]<-SevenGo_RT_Trim

#Exclude trials where RT=0, don't count
gonogo<-subset(gonogo, GoDisplay_RT_Trim!="0")

#Now that all outlier RTs have been replaced by NA...
#Group DESCRIPTIVES
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

desc2<-getDescriptivesByGroup(gonogo, c("GoDisplay_RT_Trim", "failedResp", "falseAlarm"), "Subject")

#Check who was an outlier
j<-subset(desc, var=="falseAlarm")
a<-boxplot(j$mean)
print(a) #two outliers at 64% subject-7 and 59%, subject-99 percent false alarms%
#drop em:

gonogo<-subset(gonogo, Subject!="99" & Subject!="7")
#PLot observed values

#Rt Observed
g<-ggplot(gonogo, aes(x=BlockType, y=GoDisplay_RT_Trim), colour='red', size=3)+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5),
               size=2.5)

z<-ggplot(gonogo, aes(x=BlockType,y=GoDisplay_RT_Trim))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2.5)+
  scale_color_brewer("Block", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text( vjust=.5,size=20, 
                                  colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90, vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Go/No-Go Reaction Times")+xlab("\nBlockType")+ylab("Mean Reaction Time (ms)\n")
z+ggtitle("Go/No-Go Reaction Time")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))

#RT Model
fm1<-lmer(GoDisplay_RT_Trim~BlockType*Trial+(1|Subject), data=gonogo)
cm<-lmerCellMeans(fm1)
d <- ggplot(cm, aes(x=BlockType, y=GoDisplay_RT_Trim, 
                    ymin=GoDisplay_RT_Trim-se, ymax=GoDisplay_RT_Trim+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(d)
#Model plot matches well.
  
#FRS observed ; NEed to do with dataset that doesnt exclude RT=0 since thats what failed response is
pp<-ggplot(gonogo, aes(x=BlockType,y=failedResp))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Block", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90,vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Failed Responses Go/No-Go")+xlab("\nBlock Type")+ylab("Mean % of Failed Responses\n")
pp+ggtitle("Go/ No-Go Failed Responses")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))
#FR Model
fm2<-glmer(failedResp~BlockType+(1+Trial | Subject), data=gonogo)
#Using new lmerCellMeans update:
cm2<-lmerCellMeans(fm2)
e <- ggplot(cm2, aes(x=BlockType, y=failedResp, 
                    ymin=failedResp-se, ymax=failedResp+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Block", palette="Set2")+
   theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(e)
##Model without Trial in it looks most like observed failed resp plot

#False Alarms Observed
tt<-ggplot(gonogo, aes(x=BlockType,y=falseAlarm))+
  stat_summary(fun.data="mean_cl_boot", position=position_dodge(width=0.5), size=2)+
  scale_color_brewer("Block", palette="Set2")+theme_bw()+
  theme(axis.title.x=element_text(face="bold",colour="Black", size=24),
        axis.text.x=element_text(vjust=.5,size=20, 
                                 colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=24),
        axis.text.y=element_text(angle=90,vjust=.5,size=20, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="False Alarms (Errors) Go/No-Go")+xlab("\nBlock Type")+ylab("Mean % false alarms\n")
tt+ggtitle("Go/ No-Go False Alarms")+theme(plot.title=element_text(lineheight=.8,face="bold", size=30))

#FAs Model
fm2<-glmer(falseAlarm~BlockType+(1+Trial | Subject), data=gonogo)
#Using new lmerCellMeans update:
cm2<-lmerCellMeans(fm2)
e <- ggplot(cm2, aes(x=BlockType, y=falseAlarm, 
                     ymin=falseAlarm-se, ymax=falseAlarm+se)) +
  geom_pointrange(size=1.5, position=position_dodge(width=.5)) +
  scale_color_brewer("Block", palette="Set2")+
  theme(axis.title.x=element_text(face="bold",colour="#990000", size=18),
        axis.text.x=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(axis.title.y=element_text(face="bold",colour="#990000", size=18),
        axis.text.y=element_text(angle=90, vjust=.5,size=16, colour="black"))+
  theme(panel.grid.minor=element_blank())
print(e)

#fit well
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

traits <- c("T_Propriety", "T_Detachment", "T_SelfHarm","T_Manipulativeness","T_Aggression", "T_Dependency", "T_LowSelfEsteem","MPS_alT","MPS_agT","MPS_srT", "MPS_wbT", "MPS_spT", "MPS_acT",
"MPS_scT", "MPS_clT", "MPS_haT", "MPS_tdT", "MPS_abT", "MPS_PET", "MPS_NET", "MPS_COT",
"MPS_uvT", "MPS_trT", "ATQ_ECTotal","ATQ_InhibControl","ATQ_AttControl","T_Mistrust","T_Entitlement")

gono <- c("MeanGoRT_1","sumFalseAlarms_1")

save(alldat,file='alldat.Rdata')
df<-alldat
load('alldat.Rdata')
cmat.list <- corwithtarget(df,target=gono,withvars=traits)
cmat.3dmat <- do.call(abind, list(cmat.list, along=0))
cmat.r <- melt( t(cmat.3dmat[,1,]),value.name='r' )
cmat.p <- melt( t(cmat.3dmat[,2,]),value.name='p' )
cmat <- merge(cmat.r,cmat.p, by=c("X1", "X2"))
names(cmat)[1:4] <- c('traits','gono', 'r', 'p')
ge<-ggplot( cmat, aes(x=traits,y=gono,fill=p,label=r) ) + geom_tile() + scale_fill_continuous(high='#000000',low='#FF0000', limits=c(0,.10) ) + geom_text(color=I('white'))+theme(axis.title.x=element_text(face="bold",colour="Black", size=18),
                                                                                                                                                                                  axis.text.x=element_text(angle=90, vjust=.5,size=12, 
                                                                                                                                                                                                           colour="black", ))+
  theme(axis.title.y=element_text(face="bold",colour="Black", size=18),
        axis.text.y=element_text( vjust=.5,size=12, colour="black"))+
  theme(panel.grid.minor=element_blank())+
  labs(title="Personality Traits and Go/No-Go Correlations")+xlab("Trait")+ylab("Response Inhbition Performance")
ge+scale_y_discrete(limits=rev((levels=c("MeanGoRT_1",
                                         "sumFalseAlarms_1"))),
                    labels=rev(c("OneGo RTs",
                                 "OneGo False Alarms")))+
  scale_x_discrete(labels=c("ATQ Attention Control",
                            "ATQ Effortful Control", "ATQ Inhibitory Control",
                            "MPQ Absorption", "MPQ Achievement",
                            "MPQ Aggression", "MPQ Alienation",
                            "MPQ Control", "MPQ Constraint", 
                            "MPQ Harm Avoidance",
                            "MPQ Neg Emotionality", "MPQ Pos Emotionality",
                            "MPQ Social Closeness", "MPQ Social Potency",
                            "MPQ Stress Reaction", "MPQ Traditionalism",
                            "MPQ TRIN", "MPQ Unlikely Virtues", "MPQ Well Being",
                            "SNAP Aggression", "SNAP Dependency",
                            "SNAP Detachment", "SNAP Enitlement",
                            "SNAP Low Self-Esteem", "SNAP Manipulativeness",
                            "SNAP Mistrust",
                            "SNAP Propriety", "SNAP Self-Harm"))


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
