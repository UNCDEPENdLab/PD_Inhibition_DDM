setwd("~/Box Sync/DEPENd/Projects/PD_Inhibition_DDM/Data/SAS Originals/")
library(sas7bdat)
library(plyr)
library(dplyr)
library(scales)
library(psych)
library(retimes)
Flanker<- read.sas7bdat("flanker.sas7bdat")
#------------------------------------------------------------------------------------------------------------
##format to HDDM specifications
colnames(Flanker)[colnames(Flanker)=="TrialSlide_ACC"] <- "response"
colnames(Flanker)[colnames(Flanker)=="TrialSlide_RT"] <- "rt"

#Flanker <- Flanker[complete.cases(Flanker[,c(5, 7)]),]    #drop cases where either rt or accuracy are NA
#------------------------------------------------------------------------------------------------------------

#To rename all the blocks, such as C1,I1,I2,C2
#library(dplyr)
Flanker <- Flanker %>% group_by(Subject) %>% mutate(Block = sapply(1:length(Running), function(row) {
  if (is.nan(.[[row, "CongruentTrials_Cycle"]])) {
    return(paste0("I", .[[row, "IncongruentTrials_Cycle"]]))
  } else {
    return(paste0("C", .[[row, "CongruentTrials_Cycle"]]))
  }
})) %>% ungroup()
#------------------------------------------------------------------------------------------------------------
##Drop unneeded variables
Flanker <- Flanker %>% select(-Clock_StartTimeOfDay, -SessionDate, -SessionTime, -Central, -CongruentTrials_Cycle, -CongruentTrials_Sample, 
                              -Correct, -Flanker, -IncongruentTrials_Cycle, -IncongruentTrials_Sample, -Procedure, -Running, -TrialSlide_CRESP,
                              -TrialSlide_DurationError, TrialSlide_OnsetDelay, -TrialSlide_RESP, -TrialSlide_RTTime, -computer)

Flanker <- data.frame(Flanker)
#------------------------------------------------------------------------------------------------------------
##creates new column Flanker$Stim that indicates direction of central (target) stimulus and Flankers
Flanker <- Flanker %>% mutate(CongruentStim=recode(CongruentTrials, `1`="CLFL", `2`="CLFR", `3`="CRFR", `4`="CRFL"),
                              IncongruentStim=recode(IncongruentTrials, `1`="CLFL", `2`="CLFR", `3`="CRFR", `4`="CRFL"))
##CLFR = central stim left, flankers right; CRFL = central stim right, flankers left etc

Flanker$Stim <- apply(Flanker[,c("CongruentStim", "IncongruentStim")], 1, function(row) {
  ifelse(is.na(row["CongruentStim"]), row["IncongruentStim"], row["CongruentStim"])
})

Flanker[,8] <- Flanker[,8]*.001

Flanker <- Flanker %>% select(-CongruentTrials, -IncongruentTrials, -CongruentStim, -IncongruentStim)
head(Flanker)
summary(Flanker$rt)      #includes min, max, and number of NA's

#To get the mean accuracy and rt (across trials) for each participant
#library(scales)
Flanker.sum <- cbind(aggregate(response ~ Subject, FUN = mean, data = Flanker), aggregate(rt ~ Subject, FUN = mean, data = Flanker))
Flanker.sum <- Flanker.sum[,c(1,2,4)]

#quality check: how many trials have subjects completed?
ddply(Flanker,.(Subject),nrow)   
#quality check: how many trials in each block? 
##when we write a function for this, we want to make sure "expected trials" and "expected trials per block" are arguments in our function
##use stopifnot()
#count(Flanker, c("Subject", "Block"))
#quality check: I & C trials by subject and by block, check for .3 and .7 in percentage row
#Flanker.Blocks <- count(Flanker, c("Subject", "Block", "Incongruent"))
#Flanker.Blocks$percentage <- Flanker.Blocks$freq / 40

#10/18/2016
#Create a new df and get the accuracy for each participant by each block 
#Flanker.Blocks_accuracy <- count(Flanker, c("Subject", "Block"))
#n <- 40
#Flanker.Blocks_accuracy$accuracy <-  aggregate(Flanker$response,list(rep(1:(nrow(Flanker)%/%n+1),each=n,len=nrow(Flanker))),mean, na.rm = TRUE)[-1];
#Flanker.Blocks_accuracy$accuracy$x <- percent(Flanker.Blocks_accuracy$accuracy$x)

#10/20/2016
#Get the accuracy of blocks by each participant
Flanker.sum1 <- Flanker[grep("C1", Flanker$Block), ]
Flanker.sum3 <- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$C1_Accuracy <- Flanker.sum3$response

Flanker.sum1 <- Flanker[grep("C2", Flanker$Block), ]
Flanker.sum3 <- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$C2_Accuracy <- Flanker.sum3$response

Flanker.sum1 <- Flanker[grep("I1", Flanker$Block), ]
Flanker.sum3 <- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$I1_Accuracy <- Flanker.sum3$response

Flanker.sum1 <- Flanker[grep("I2", Flanker$Block), ]
Flanker.sum3 <- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$I2_Accuracy <- Flanker.sum3$response

#Get the accuracy of Incongruent and congruent
Flanker.sum1 <- Flanker[grep("0", Flanker$Incongruent), ]
Flanker.sum3<- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$Incongruent_Accuracy <- Flanker.sum3$response

Flanker.sum1 <- Flanker[grep("1", Flanker$Incongruent), ]
Flanker.sum3 <- aggregate(response ~ Subject, FUN = mean, data = Flanker.sum1)
Flanker.sum$Congruent_Accuracy <- Flanker.sum3$response

#11/1/2016
#library(psych)
#Exgaussian mean
#library(retimes)
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Total_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("C1", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$C1_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("C2", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$C2_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("I1", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$I1_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("I2", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$I2_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("0", Flanker$Incongruent), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Incongruent_EM <- Flanker.sum2$mu + Flanker.sum2$tau

Flanker.sum1 <- Flanker[grep("1", Flanker$Incongruent), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Congruent_EM <- Flanker.sum2$mu + Flanker.sum2$tau

#Exgaussian Standard Deviation
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Total_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("C1", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$C1_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("C2", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$C2_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("I1", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$I1_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("I2", Flanker$Block), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$I2_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("0", Flanker$Incongruent), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Incogruent_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

Flanker.sum1 <- Flanker[grep("1", Flanker$Incongruent), ]
Flanker.sum2 <- aggregate(response ~ Subject, FUN = mexgauss, data = Flanker.sum1)
Flanker.sum2 <- data.frame(Flanker.sum2$response[c(1:111, 0),])
Flanker.sum$Congruent_ESD <- sqrt(Flanker.sum2$sigma^2 + Flanker.sum2$tau^2)

# 11/7/2016
#Boxplot for subject 1
Flanker.sum3

#participant1 <- Flanker %>% select(-response, -TrialSlide_OnsetDelay, -TrialSlide_OnsetTime, -CongruentBlock, -Stim)
participant1<-Flanker[1:160,]
boxplot(participant1[,6]~participant1[,2],main="Sub 1", xlab="Time(rt)", ylim=c(0,1), las=1, horizontal=TRUE)

#loop for making boxplots 
par(mfrow = c(1,1))
#pdf("Boxplots by condition.pdf", width=16, height=10)
for(i in seq(from=1, to=17760, by=160)){
  boxplot(Flanker[i:(i+159),6]~Flanker[i:(i+159),2],main=paste0("Subj_",Flanker[i, 1]), xlab="Time(rt)", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)
  
  i <- i+159
}
#dev.off()

#11/15/2016
Flanker.sum1 <- Flanker[grep("0", Flanker$Incongruent), ]
boxplot(Flanker.sum1[,6]~Flanker.sum1[,1],xlab="Time(rt)", main = "Congruent", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

Flanker.sum1 <- Flanker[grep("1", Flanker$Incongruent), ]
boxplot(Flanker.sum1[,6]~Flanker.sum1[,1],xlab="Time(rt)", main = "Incongruent", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

#Flanker$rt1 <-winsor(Flanker$rt, trim = 0.02, na.rm = TRUE)
#Flanker$rt <-winsor(Flanker$rt, trim = 0.02, na.rm = TRUE)

###winsorising based on MH's winsorize_vec function defined above

#The original data
bycondition1 <- split(Flanker, Flanker$Incongruent)
congruent <- bycondition1[["0"]]
incongruent <- bycondition1[["1"]]

congruent_CC <- congruent[grep("1", congruent$CongruentBlock), ]
boxplot(congruent_CC$rt~congruent_CC$Subject,xlab="Time(rt)", main = "CC_normal", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

congruent_CI <- congruent[grep("0", congruent$CongruentBlock), ]
boxplot(congruent_CI$rt~congruent_CI$Subject,xlab="Time(rt)", main = "CI_normal", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

incongruent_II <- incongruent[grep("1", incongruent$CongruentBlock), ]
boxplot(incongruent_II$rt~incongruent_II$Subject,xlab="Time(rt)", main = "II_normal", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

incongruent_IC <- congruent[grep("0", incongruent$CongruentBlock), ]
boxplot(incongruent_IC$rt~incongruent_IC$Subject,xlab="Time(rt)", main = "IC_normal", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

#Without outliers 
# boxplot(congruent_CC$rt~congruent_CC$Subject,xlab="Time(rt)", main = "CC_no outliers", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE, outline = FALSE)
# boxplot(congruent_CI$rt~congruent_CI$Subject,xlab="Time(rt)", main = "CI_no outliers", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE, outline = FALSE)
# boxplot(incongruent_II$rt~incongruent_II$Subject,xlab="Time(rt)", main = "II_no outliers", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE, outline = FALSE)
#boxplot(incongruent_IC$rt~incongruent_IC$Subject,xlab="Time(rt)", main = "IC_no outliers", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE, outline = FALSE)


#The data after winsorized
rtecdf <- ecdf(Flanker$rt)
q98 <- quantile(Flanker$rt, c(.98), na.rm=TRUE)
q02 <- quantile(Flanker$rt, c(.02), na.rm=TRUE)

winsorize_vec <- function(v, cutoff) {
  #detect whether cutoff is high or low
  q50 <- quantile(v, .50, na.rm=TRUE)
  if (cutoff > q50) { high <- TRUE } else { high <- FALSE }
  if (high) {
    v[v >= cutoff] <- max(v[v < cutoff], na.rm=TRUE)    
  } else {
    v[v <= cutoff] <- min(v[v > cutoff], na.rm=TRUE)    
  }
  return(v)
}

Flanker<- Flanker %>% group_by(Subject) %>% mutate(rt_98wins = winsorize_vec(rt, q98)) %>% ungroup()

#changed <- filter(congruent, rt_98wins != rt) %>% select(rt_98wins, rt) %>% print(n=100)
#length(na.omit(congruent$rt))
bycondition1 <- split(Flanker, Flanker$Incongruent)
congruent <- bycondition1[["0"]]
incongruent <- bycondition1[["1"]]

congruent_CC <- congruent[grep("1", congruent$CongruentBlock), ]
boxplot(congruent_CC$rt_98wins~congruent_CC$Subject,xlab="Time(rt)", main = "CC_winsorized", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

congruent_CI <- congruent[grep("0", congruent$CongruentBlock), ]
boxplot(congruent_CI$rt_98wins~congruent_CI$Subject,xlab="Time(rt)", main = "CI_winsorized", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

incongruent_II <- incongruent[grep("1", incongruent$CongruentBlock), ]
boxplot(incongruent_II$rt_98wins~incongruent_II$Subject,xlab="Time(rt)", main = "II_winsorized", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)

incongruent_IC <- congruent[grep("0", incongruent$CongruentBlock), ]
boxplot(incongruent_IC$rt_98wins~incongruent_IC$Subject,xlab="Time(rt)", main = "IC_winsorized", ylim=c(0,1), las=1, horizontal=TRUE, plot=TRUE)



