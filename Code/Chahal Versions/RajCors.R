
library(RColorBrewer)
library(gplots)
library(abind)
library(Hmisc)
source('/home/rajchahal/Desktop/WPIC/flanker_Control/corwithtarget.R')

traits <- c("T_Propriety", "T_Detachment", "T_SelfHarm", "T_NegativeTemperament","T_Manipulativeness","T_Aggression", "T_Dependency", "T_LowSelfEsteem", "MPS_NET","MPS_alT","MPS_agT","MPS_srT","ATQ_ECTotal","ATQ_InhibControl","ATQ_AttControl","T_Mistrust","T_Entitlement")
flank <- c("Sad_I_PctInaccurate","Ang_I_PctInaccurate", "I_Sad","I_Ang", "Flankpinac")

df <- load('cmat.Rdata')
cmat.list <- corwithtarget(df,target=flank,withvars=traits)
cmat.3dmat <- do.call(abind, list(cmat.list, along=0))
cmat.r <- melt( t(cmat.3dmat[,1,]),value.name='r' )
cmat.p <- melt( t(cmat.3dmat[,2,]),value.name='p' )
cmat <- merge(cmat.r,cmat.p)
names(cmat)[1:2] <- c('Trait','Flank')
ggplot( cmat, aes(x=Trait,y=Flank,fill=p,label=r) ) + geom_tile() + scale_fill_continuous(high='#000000',low='#FF0000', limits=c(0,.05) ) + geom_text(color=I('white'))

save(cmat,file='cmat.Rdata')

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
