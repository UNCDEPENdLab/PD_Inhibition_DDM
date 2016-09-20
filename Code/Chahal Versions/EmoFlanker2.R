library(ez, ggplot2)
library(sas7bdat)
library(lme4)
library(plotrix)
library(plyr)
library(reshape)
library(multcomp)


eflank<-read.sas7bdat('/emotflanker.sas7bdat',sep=",",header=T)
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



#Group Descriptives
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

desc<-getDescriptivesByGroup(eflank, c("FlankerDisplay_ACC", "FlankerDisplay_RT_Trim"), "Subject")

#To check descriptives on incongruent trials only:
eflankincong<-subset(eflank, congruent=="incongruent")
descincong<-getDescriptivesByGroup(eflankincong, c("FlankerDisplay_ACC", "FlankerDisplay_RT_Trim"), "Subject")

#Plot Accuracy on All Trials to find outliers
accdesc<-subset(desc, var=="FlankerDisplay_ACC")
attach(accdesc)
boxplot(mean)
print(boxplot(mean))

#Plot Accuracy on Incongruent trials to find outliers
accdescincong<-subset(descincong, var=="FlankerDisplay_ACC")
attach(accdescincong)
boxplot(mean)
print(boxplot(mean))


#Dropping subjects: from analysis: In orig manuscript, subs with 30& or more inacc reponses were excluded.
#High percentage of inaccurate responses (59% or more) were excluded
#Dropping subjects 61 and 98
#additional Dropped subjects by Raj: 
#Dropped Subjects: 96, 39, 50, 61, 98
eflank<-subset(eflank, Subject!="98" & Subject!="39" & Subject!="50" & Subject!="61" & Subject!="96")
#Other Descriptives: Means within conditions
#RT means
RTmeanseflank<-ddply(eflank, .(Subject, congruent, emotion, EmotionPairs), 
                   summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
#ACC means
ACCeflank<-ddply(eflank, .(Subject, congruent, emotion, EmotionPairs), 
                   summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
ACCeflank
AccEmoPairs<-ddply(eflank, .(Subject, EmotionPairs), 
                      summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccEmoPairs
AccEmo<-ddply(eflank, .(Subject, emotion), 
                   summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccEmo
AccCong<-ddply(eflank, .(Subject, congruent), 
              summarize, mean_ACC=mean(FlankerDisplay_ACC, na.rm=TRUE))
AccCong


#Subsets for further analysis:

RTflank<-subset(eflank, FlankerDisplay_RT_Trim>0)
RTflank_noAS<-subset(RTflank, EmotionPairs!="Ang-Sad")
RTflank_noASI<-subset(RTflank_noAS, congruent=="incongruent")
eflank_noAS<-subset(eflank, EmotionPairs!="Ang-Sad")
eflank_noASI<-subset(eflank_noAS, congruent=="incongruent")

#####################PLOTS#####################################
#ggplot rt by  emo pair and emo
sd(FlankerDisplay_RT_Trim, na.rm=TRUE)
unique(eflank$Subject)
sqrt(114)
#Standard Error= Standard Deviation Divided by Square Root of Sample Size (191.487/10.677)
se<-17.93443
ggplot(eflank, aes(x=emotion, y=FlankerDisplay_RT_Trim, ymin=FlankerDisplay_RT_Trim - se, 
                   ymax=FlankerDisplay_RT_Trim + se, color=EmotionPairs,))+geom_boxplot()

#EZ Plot RT by Emo
#All Trials
rtplot1=ezPlot(data=RTflank,dv=.(FlankerDisplay_RT_Trim), 
               wid=.(Subject),  
               between=NULL,
               within=emotion,
               within_full=emotion,
               x=.(emotion),
               do_lines=TRUE, 
               do_bars=TRUE, 
               x_lab='emo', 
               y_lab='RT',
)
print(rtplot1)

#RT on correct responses by emotion and congruency
rtplot1=ezPlot(data=RTflank[RTflank$FlankerDisplay_ACC==1,],dv=.(FlankerDisplay_RT_Trim), 
              wid=.(Subject),  
              between=NULL,
              within=.(emotion, congruent),
              within_full=emotion,
              split=congruent,
              x=.(emotion),
              do_lines=TRUE, 
              do_bars=TRUE, 
              x_lab='Emotion', 
              y_lab='RT',
              split_lab='congruency'
)
rtplot1=rtplot1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(rtplot1)

#RT Means by Emotion Pair and Congruency on correct trials
rtplot2=ezPlot(data=RTflank,dv=.(FlankerDisplay_RT_Trim), 
               wid=.(Subject),  
               between=NULL,
               within=.(EmotionPairs, congruent),
               within_full=EmotionPairs,
               split=congruent,
               x=.(EmotionPairs),
               do_lines=TRUE, 
               do_bars=TRUE, 
               x_lab='Emotion Pairs', 
               y_lab='RT',
               split_lab='congruency'
)
rtplot2=rtplot2+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(rtplot2)

#RT Means by Emotion Pair and Congruency on correct trials without Angry-Sad block
RTflank_noAS<-subset(RTflank, EmotionPairs!="Ang-Sad")
#str(RTflank)
rtplot3=ezPlot(data=RTflank_noAS,dv=.(FlankerDisplay_RT_Trim), 
               wid=.(Subject),  
               between=NULL,
               within=.(EmotionPairs, congruent),
               within_full=EmotionPairs,
               split=congruent,
               x=.(EmotionPairs),
               do_lines=TRUE, 
               do_bars=TRUE, 
               x_lab='Emotion Pairs', 
               y_lab='RT',
               split_lab='congruency'
)
rtplot3=rtplot3+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(rtplot3)

#RT Means by Emotion and congruency without Angry-Sad block
RTflank_noAS<-subset(RTflank, EmotionPairs!="Ang-Sad")
rtplot4=ezPlot(data=RTflank_noAS[RTflank_noAS$FlankerDisplay_ACC==1,],dv=.(FlankerDisplay_RT_Trim), 
               wid=.(Subject),  
               between=NULL,
               within=.(emotion, congruent),
               within_full=emotion,
               split=congruent,
               x=.(emotion),
               do_lines=TRUE, 
               do_bars=TRUE, 
               x_lab='Emotion', 
               y_lab='RT',
               split_lab='congruency'
)
rtplot4=rtplot4+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(rtplot4)

#EZ Plot Acc by Emo on all trials
#percentage of correct trials per subject on conditions

#Accuracy Means by emotion/congruency
accplot1=ezPlot(data=eflank,dv=.(FlankerDisplay_ACC), 
                wid=.(Subject),  
                between=NULL,
                within=.(emotion, congruent),
                within_full=emotion,
                split=congruent,
                x=.(emotion),
                do_lines=TRUE, 
                do_bars=TRUE, 
                x_lab='Emotion', 
                y_lab='Accuracy',
                split_lab='congruency'
)
accplot1=accplot1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(accplot1)

#Accuracy Means by Emotion Pairs and Congruency
accplot2=ezPlot(data=eflank,dv=.(FlankerDisplay_ACC), 
                wid=.(Subject),  
                between=NULL,
                within=.(EmotionPairs, congruent),
                within_full=EmotionPairs,
                split=congruent,
                x=.(EmotionPairs),
                do_lines=TRUE, 
                do_bars=TRUE, 
                x_lab='Emotion Pairs', 
                y_lab='Accuracy',
                split_lab='congruency'
)
accplot2=accplot2+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(accplot2)

#Accuracy Means by Emotion Pairs and Congruency without Angry-Sad block
eflank_noAS<-subset(eflank, EmotionPairs!="Ang-Sad")
accplot3=ezPlot(data=eflank_noAS,dv=.(FlankerDisplay_ACC), 
                wid=.(Subject),  
                between=NULL,
                within=.(EmotionPairs, congruent),
                within_full=EmotionPairs,
                split=congruent,
                x=.(EmotionPairs),
                do_lines=TRUE, 
                do_bars=TRUE, 
                x_lab='Emotion Pairs', 
                y_lab='Accuracy',
                split_lab='congruency'
)
accplot3=accplot3+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(accplot3)

#Acc means by emotion and cong withuot angry-sad block
accplot4=ezPlot(data=eflank_noAS,dv=.(FlankerDisplay_ACC), 
                wid=.(Subject),  
                between=NULL,
                within=.(emotion, congruent),
                within_full=emotion,
                split=congruent,
                x=.(emotion),
                do_lines=TRUE, 
                do_bars=TRUE, 
                x_lab='Emotion', 
                y_lab='Accuracy',
                split_lab='congruency'
)
accplot4=accplot4+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
print(accplot4)

##Plot RT by Emo and Emo Pair (And possibly congruency)

RTmeanseflank<-ddply(eflank, .(Subject, congruent, emotion, EmotionPairs), 
                     summarize, mean_RT=mean(FlankerDisplay_RT_Trim, na.rm=TRUE))
ggplot(RTmeanseflank, aes(x=emotion, y=mean_RT)+ geom_point())










##Plot Mean ACC by Emo and Emo pair(and possibly congruency)



##########################All Linear Models##############################

#Linear model for RT on all trials
summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+emotion/EmotionPairs*congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))

cm <- lmerCellMeans(fm1)
plot(cm)

#new lmer cell means model works for both directions of nesting: emotion/EmoPair and EmoPair/emotion
#then subset out the nonsense predictions
#ggplot(cm)
#NExt Steop: plot the lmer cell means (predicted) and plot the observed data( already done):
  # See how well the predicted plots match the observed plots: they should be simialar for good models
lmerCellMeans(fm1)
#Linear model for ACC on all trials
gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + SubTrial + congruent + (1+SubTrial | Subject), eflank, family="binomial")
summary(gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + Trial + congruent + (1+Trial | Subject), eflank, family="binomial"))
gm1diag<-glm.diag(gm1)
glm.diag.plots(gm1, gm1diag)

##RERRUN THIS^^^^^^^^^

#No angry-sad blocks, linear models

#ezMixed for RT without ang-sad block
ezrt<-ezMixed(data=eflank_noAS, dv=,(FlankerDisplay_RT_Trim), random=.(Subject, Trial), 
              fixed=.(congruent, emotion, EmotionPairs), family=gaussian)
print(ezrt$summary)

# Linear Model LMER of RT without ang-sad block
summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion/EmotionPairs+congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank_noAS))

#Linear Model LMER of RT without ang-sad block on incongruent trials only
eflank_noASIncong<-subset(eflank_noAS, congruent=="incongruent")
summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion/EmotionPairs+(1|Subject)+Trial+(1+Trial|Subject), data=eflank_noASIncong))

#Accuracy models without ang-sad and incongruent only: Ezmixed, glmer: 

#ezMixed for ACC without ang-sad block
ezacc<-ezMixed(data=eflank_noAS, dv=,(FlankerDisplay_ACC), random=.(Subject, Trial), 
               fixed=.(congruent, emotion, EmotionPairs), family=binomial)
print(ezacc$summary)

#Linear model GLM for ACC without ang-sad block
summary(gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + Trial + congruent + (1+Trial | Subject), eflank_noAS, family="binomial"))

#Linear model GLM for ACC without ang-sad block on incongruent trials only
summary(gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + Trial + (1+Trial | Subject), eflank_noASIncong, family="binomial"))


#couldnt figure out... 
#Function for OneWayAnova
oneWayANOVA <- function(dataset, outcome, group) {
  require(gplots)
  aovformula <- as.formula(paste(outcome, "~", group))
  plotmeans(aovformula, data=dataset)
  
  cat("\n\n-----\nGroup means and SDs\n\n")
  meanSd <- tapply(dataset[[outcome]], dataset[[group]], mean, na.rm=TRUE)
  meanSd <- rbind(meanSd, tapply(dataset[[outcome]], dataset[[group]], sd, na.rm=TRUE)) 
  row.names(meanSd) <- c("Mean", "SD")
  print(meanSd)
  
  cat("\n\n-----\nOne-way ANOVA\n")
  print(summary(anovaModel <- aov(aovformula, dataset)))
  
  #  cat("\n\n-----Planned contrasts\n")
  #  cat("\nAll with class 1\n")
  #  #setup planned contrasts here (fit.contrasts requires one to test fewer contrasts than the num factor levels)
  #  plannedContrasts <- rbind( "1 vs 2" = c(-1, 1, 0, 0, 0),
  #      "1 vs 3" = c(-1, 0, 1, 0, 0),
  #      "1 vs 4" = c(-1, 0, 0, 1, 0),
  #      "1 vs 5" = c(-1, 0, 0, 0, 1))
  #  
  #  print(fit.contrast(anovaModel, group, plannedContrasts))
  
  #  cat("\n3v5, 2v4\n")  	
  #  plannedContrasts <- rbind( "3 vs 5" = c(0, 0, -1, 0, 1),
  #      "2 vs 4" = c(0, -1, 0, 1, 0))
  #  
  #  print(fit.contrast(anovaModel, group, plannedContrasts))
  
  cat("\n\n-----\nTukey's HSD Post Hocs\n")
  print(TukeyHSD(anovaModel, ordered=FALSE))
}
oneWayANOVA ("eflank", "FlankerDisplay_RT_Trim", "emotion")

#potentially add ability to get empirical bayes trajectories per person using coef?
lmerCellMeans <- function(lmerObj, divide=NULL, n.divide=3, divide.prefix=TRUE, n.cont=30, fixat0=NULL,
                          yoked=NULL) { 
  #print cell means for lmer by expanding level combinations and multiplying against fixed effects  
  predNames <- attr(terms(lmerObj), "term.labels")
  whichME <- attr(terms(lmerObj), "order")
  predNames <- predNames[whichME==1]
  
  #sort yoked predictors to the end
  if (!is.null(yoked)) {
    yoked.tf <- as.numeric(predNames %in% sapply(yoked, "[[", "source"))
    predNames <- predNames[rev(order(yoked.tf))]
  }
  
  predData <- list()
  #divide into categorical and continuous predictors, determine which continuous predictors to discretize
  for (f in predNames) {
    if (f %in% fixat0) { predData[[f]] <- 0 #compute model prediction when this term is 0
                         #} else if (attr(terms(lmerObj), "dataClasses")[f] == "factor") {
    } else if (attr(attr(lmerObj@frame, "terms"), "dataClasses")[f] == "factor") {
      predData[[f]] <- levels(lmerObj@frame[[f]])  
    } else {
      if (f %in% divide) {
        #divide into -1 SD, M, + 1 SD; or -2SD, -1SD, M, +1SD, +2SD
        fsd <- sd(lmerObj@frame[[f]], na.rm=TRUE)
        fm <- mean(lmerObj@frame[[f]], na.rm=TRUE)
        predData[[f]] <- if (n.divide==3) { c(fm-fsd, fm, fm+fsd)
        } else { c(fm-fsd*2, fm-fsd, fm, fm+fsd, fm+fsd*2) }
      } else if (!is.null(yoked) && any(pmatch <- which(f %in% sapply(yoked, "[[", "dest") == TRUE))) {
        next #do nothing -- must propagate after expand.grid
      }	else {
        if (!is.null(names(n.cont))) {
          #Named vector specifying number of points to predict for each IV
          if (is.na(n.cont[f])) stop("Cannot locate number of continuous pred points for: ", f)
          predData[[f]] <- seq(min(lmerObj@frame[[f]], na.rm=TRUE), max(lmerObj@frame[[f]], na.rm=TRUE), length=n.cont[f])          
        } else {
          #treat as truly continuous predictor and compute models estimates across the range of observed values
          predData[[f]] <- seq(min(lmerObj@frame[[f]], na.rm=TRUE), max(lmerObj@frame[[f]], na.rm=TRUE), length=n.cont)
        }
      }
    }
  }
  
  #dependent variable
  dvname <- as.character(terms(lmerObj)[[2]])
  
  #populate the model-predicted estimates with 0s prior to building model matrix
  predData[[dvname]] <- 0
  
  #Develop a grid 
  predData <- do.call(expand.grid, list(predData))
  
  #handle yoked predictors after the expansion
  if (!is.null(yoked)) {
    for (i in 1:length(yoked)) {
      predData[[yoked[[i]]$dest]] <- sapply(predData[[yoked[[i]]$source]], yoked[[i]]$transform)
    }
  }
  
  browser()
  
  mm <- model.matrix(terms(lmerObj),predData)
  
  predData[[dvname]] <- mm %*% fixef(lmerObj)
  
  pvar1 <- diag(mm %*% tcrossprod(vcov(lmerObj),mm))
  tvar1 <- pvar1+VarCorr(lmerObj)[[1]][1] #assumes that the first element in VarCorr is subject
  
  #confidence and prediction intervals
  predData <- data.frame(
    predData, se=sqrt(pvar1),
    plo = predData[[dvname]]-2*sqrt(pvar1),
    phi = predData[[dvname]]+2*sqrt(pvar1),
    tlo = predData[[dvname]]-2*sqrt(tvar1),
    thi = predData[[dvname]]+2*sqrt(tvar1)
  )
  
  for (f in divide) {
    if (n.divide==3) { flevels <- c("-1 SD", "Mean", "+1 SD")
    } else if (n.divide==5) { flevels <- c("-2SD", "-1 SD", "Mean", "+1 SD", "+2 SD") }
    if (divide.prefix) flevels <- paste(Hmisc::capitalize(f), flevels, sep=": ")
    predData[[f]] <- factor(predData[[f]], levels=sort(unique(predData[[f]])), labels=flevels)
  }
  
  return(predData)
}

oneWayANOVA("eflank_noAS", "FlankerDisplay_RT_Trim", "emotion")



#Simple Effects Matrix Function
getSimpleEffectsMatrix <- function(coefNames, iv.slice, iv.test, levels.slice, levels.test) {
  require(gtools)
  
  numEffects <- choose(length(levels.test), 2)
  
  conMats <- list()
  for (l in levels.slice) {
    thisCon <- matrix(0, nrow=numEffects, ncol=length(coefNames))
    conNames <- rep(NA_character_, numEffects)
    colnames(thisCon) <- coefNames
    combs <- combinations(length(levels.test), 2)
    for (pair in 1:nrow(combs)) {
      p1 <- levels.test[combs[pair,1]]
      p2 <- levels.test[combs[pair,2]]
      me.p1 <- paste(iv.test, p1, sep="")
      me.p2 <- paste(iv.test, p2, sep="")
      
      #add main effects
      if(any(grepl(paste("^", me.p1, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair,me.p1] <- 1
      if(any(grepl(paste("^", me.p2, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair,me.p2] <- -1
      
      #add interactions
      thisSlice <- paste(iv.slice, l, sep="")
      #allow for either order (depending on how model was specified
      if (any(grepl(paste("^", me.p1, ":", thisSlice, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair, paste(me.p1, ":", thisSlice, sep="")] <- 1
      else if (any(grepl(paste("^", thisSlice, ":", me.p1, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair, paste(thisSlice, ":", me.p1, sep="")] <- 1
      if (any(grepl(paste("^", me.p2, ":", thisSlice, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair, paste(me.p2, ":", thisSlice, sep="")] <- -1
      else if (any(grepl(paste("^", thisSlice, ":", me.p2, "$", sep=""), coefNames, perl=TRUE))) thisCon[pair, paste(thisSlice, ":", me.p2, sep="")] <- -1
      conNames[pair] <- paste(p1, "-", p2, "AT", l)
    }
    rownames(thisCon) <- conNames
    conMats[[l]] <- thisCon
  }
  
  return(conMats)
  
}

#contrasts <- getSimpleEffectsMatrix(names(fixef(lmerAmp)), 
#    iv.slice="signal", iv.test="pipeline", levels.slice=levels(subsamp$signal), levels.test=levels(subsamp$pipeline))

#l_ply(contrasts, function(lev) {
#      cat("\n  GLHT, single-step p-values\n\n")
#      print(summary(glht(lmerAmp, linfct = lev)))
#      cat("\n  MCMC p-values, 1000 samples\n\n")
#      print(estimable(multcoherMLM_hetvar_ROI_f, lev, sim.mer=TRUE, n.sim=1000))      
#    })

#RT Congruent Vs. Emotion no ang Sad simple effects
lmeEff_congemoRT <- lmer(FlankerDisplay_RT_Trim~ congruent+emotion+EmotionPairs+(1|Subject), RTflank_noAS)
                       conmat <- do.call(rbind, 
                                         getSimpleEffectsMatrix(
                                           names(fixef(lmeEff_congemoRT)), 
                                           iv.slice="congruent",
                                           iv.test="emotion", 
                                           levels.slice=levels(RTflank_noAS$congruent),
                                           levels.test=levels(RTflank_noAS$emotion)
                                         )
                       )
                       
summary(lmeEff_congemoRT)
summary(glht(lmeEff_congemoRT, linfct=conmat))

#RT Emotion Vs. congruent no ang sad simple effects *Nothing significant
lmeEff_congemopairRT <- lmer(FlankerDisplay_RT_Trim~ congruent+emotion+EmotionPairs+(1|Subject), RTflank_noAS)
conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(lmeEff_congemopairRT)), 
                    iv.slice="emotion",
                    iv.test="congruent", 
                    levels.slice=levels(RTflank_noAS$emotion),
                    levels.test=levels(RTflank_noAS$congruent)
                  )
)

summary(lmeEff_congemopairRT)
summary(glht(lmeEff_congemopairRT, linfct=conmat), na.action=na.omit)
##^^^^This is the last function you left on RC, 5/5/14 3:41 pm


#RT Emotion Vs. Emotion Pair. No ang sad simple effects
#didnt work. nonfinite list is doubtful?
lmeEff_emoemopairRT <- lmer(FlankerDisplay_RT_Trim~ congruent+emotion+EmotionPairs+(1|Subject), eflank_noAS)
conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(lmeEff_emoemopairRT)), 
                    iv.slice="emotion",
                    iv.test="EmotionPairs", 
                    levels.slice=levels(eflank$emotion),
                    levels.test=levels(eflank$EmotionPairs)
                  )
)

summary(lmeEff_emoemopairRT)
summary(glht(lmeEff_emoemopairRT, linfct=conmat))

#RT Congruent Vs Emotion Simple Effects. No Ang-Sad
lmeEff_congemoACC <- lmer(FlankerDisplay_ACC~ congruent+emotion*EmotionPairs+ (1|Subject), eflank_noAS)
conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(lmeEff_congemoACC)), 
                    iv.slice="congruent",
                    iv.test="emotion", 
                    levels.slice=levels(eflank$congruent),
                    levels.test=levels(eflank$emotion)
                  )
)

summary(glht(lmeEff_congemoACC, linfct=conmat))

#ACC Congruent Vs Emotion Pair Simple EFfects. No Ang-Sad
lmeEff_congemoACC <- lmer(FlankerDisplay_ACC~ congruent+emotion*EmotionPairs+ (1|Subject), eflank_noAS)
conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(lmeEff_congemoACC)), 
                    iv.slice="congruent",
                    iv.test="emotion", 
                    levels.slice=levels(eflank$congruent),
                    levels.test=levels(eflank$emotion)
                  )
)

summary(glht(lmeEff_congemoACC, linfct=conmat))

#Run in lab. Multiple ivs?
################################## lmer cell means matrix new model###
fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank)
summary(fm1)
conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(fm1)), 
                    iv.slice="congruent",
                    iv.test="emotion", 
                    levels.slice=levels(eflank$congruent),
                    levels.test=levels(eflank$emotion)
                  )
)

summary(glht(fm1, linfct=conmat))

conmat <- do.call(rbind, 
                  getSimpleEffectsMatrix(
                    names(fixef(fm1)), 
                    iv.slice="congruent",
                    iv.test="emotion",
                    levels.slice=levels(eflank$congruent),
                    levels.test=levels(eflank$emotion),
                  )
)
summary(fm1)
summary(glht(fm1, linfct=conmat))

summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion*emotion/EmotionPairs*congruent+(1|Subject)+SubTrial+(1+SubTrial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+emotion/EmotionPairs*congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))
#summary(fm1<-lmer(FlankerDisplay_RT_Trim~emotion+EmotionPairs+congruent+(1|Subject)+Trial+(1+Trial|Subject), data=eflank))

cm <- lmerCellMeans(fm1)
plot(cm)

#new lmer cell means model works for both directions of nesting: emotion/EmoPair and EmoPair/emotion
#then subset out the nonsense predictions
#ggplot(cm)
#NExt Steop: plot the lmer cell means (predicted) and plot the observed data( already done):
# See how well the predicted plots match the observed plots: they should be simialar for good models
lmerCellMeans(fm1)
#Linear model for ACC on all trials
gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + SubTrial + congruent + (1+SubTrial | Subject), eflank, family="binomial")
summary(gm1 <- glmer(FlankerDisplay_ACC ~ emotion/EmotionPairs + Trial + congruent + (1+Trial | Subject), eflank, family="binomial"))
gm1diag<-glm.diag(gm1)
glm.diag.plots(gm1, gm1diag)


#Into binary linear model, fit the personality measures 
#1) as fixed variables of relevant task index
#2) as dependent var with task indeces as fixed vars 
#Also do the same into generalized mixed model


#Correlate with the ATQ, SNAP2, and MPQ

