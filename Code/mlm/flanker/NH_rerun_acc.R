
if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(car, nlme, lme4, readr, tidyverse, emmeans, cowplot, MplusAutomation, knitr, sas7bdat, ggcorrplot, psych, lavaan)


#flanker data processed in flanker_mlm.Rmd
#has the worst subjects removed, rt_inv transformation, etc.
load("../../../Data/preprocessed/flanker_for_traits.RData")
flanker <- flanker %>% mutate(rt_inv=if_else(correct==0, NA_real_, rt_inv),  #treat incorrect RTs as missing
                              prev_rt_z=as.vector(scale(prev_rt)))

#basedir <- "C:/Users/timot/Documents/GitHub/PD_Inhibition_DDM"
basedir <- "~/github_repos/PD_Inhibition_DDM"

load("/Users/natehall/Downloads/flanker_for_traits.RData")

flanker <- flanker %>% mutate(rt_inv=if_else(correct==0, NA_real_, rt_inv),  #treat incorrect RTs as missing
                              rt = if_else(correct==0, NA_real_, rt),
                              prev_rt_z = as.numeric(scale(prev_rt))) %>% arrange(id, trial)

mpq <- get(load(file.path(basedir, "Data/preprocessed/MPQ_all_scored_final.RData"))) %>% 
  filter(exclude_MPQ == 0) %>% select(subject, MPS_wbR:MPS_abR)
snap <- get(load(file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData"))) %>% filter(exclude_SNAP == 0) %>% 
  select(subject, NEGTEMP:HDWK, DISINHP)
k10 <- sas7bdat::read.sas7bdat(file.path(basedir, "Data/SAS Originals/k10.sas7bdat")) %>%
  dplyr::select(subject, k10Total) %>% dplyr::rename(K10=k10Total)
stai <- sas7bdat::read.sas7bdat(file.path(basedir, "Data/SAS Originals/stai.sas7bdat")) %>%
  dplyr::select(subject, staitotal) %>% dplyr::rename(STAI=staitotal)

self_reps <- inner_join(mpq, snap, by = "subject")
self_reps <- left_join(self_reps, k10, by = "subject")
self_reps <- left_join(self_reps, stai, by = "subject")
mpluscmd <- "/Users/natehall/Desktop/Mplus/mplus"
#mpluscmd <- "/Applications/Mplus/mplus"



# Initial model in Mplus (m15)

# This model mimics m15 from flanker ACC mlms for publication. Here, I start out with a two-level model to keep life simpler, rather than having blocks nested within subjects.


# ```
# m15 <- glmer(correct ~ cond + block + prev_rt_z + prev_error +
#                (1 + prev_rt_z + cond + prev_error | id) +
#                (1 | id:block),
#              na.action = na.exclude, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 30000)),
#              data = flanker, family=binomial)
# ```


m15_mplus <- mplusObject(
  TITLE = "Equivalent model to glmer m15",
  DEFINE = "
    cond = cond - 1; ! 0=congruent, 1=incongruent
    block = block - 1; ! 0=most_incon, 1=most_con
  ",
  VARIABLE = "
    WITHIN = cond block prev_rt_z prev_error;
    USEVARIABLES = id correct cond block prev_rt_z prev_error;
    CLUSTER = id;
    CATEGORICAL = correct;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(20000) 100000;
    CHAINS=2;
    PROCESSORS=4;
  ",
  MODEL = "
    %WITHIN%
    correct ON block;
    condslo | correct ON cond;        !random slope of cond
    prevslo | correct ON prev_rt_z;   !random slope of rt autocorrelation
    errslo | correct ON prev_error;   !random slope of previous error

    %BETWEEN%
    !means of random slopes
    [condslo];
    [prevslo];
    [errslo];

    !variances of random slopes
    condslo;
    prevslo;
    errslo;

    !slope correlations
    condslo prevslo errslo WITH
       condslo prevslo errslo;

    [correct$1]; !intercept: average correct
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = flanker
)

mout <- mplusModeler(m15_mplus, dataout="flanker_acc_mplus.dat",
                     modelout="acc_m15.inp", run=TRUE, hashfilename=FALSE,
                     Mplus_command = mpluscmd)

summary(mout)
# kable(mout$results$parameters$stdyx.standardized %>%
#         filter(!paramHeader %in% c("Intercepts", "Means", "Residual Variances")))



# ACC: Antagonism + Constraint MSEM


tofactor <- self_reps %>% select(subject, MPS_agR, AGG, MISTRUST, MPS_alR, MANIP,
                                 PROPER, MPS_clR, IMPUL, MPS_tdR, MPS_acR, HDWK, K10, STAI) %>% 
  # mutate(IMPUL =-1*IMPUL) %>% #score toward constraint to make loadings upright. Update 6/25: reverse coding to scale positively with high externalizing
  mutate(PROPER = -1*PROPER, MPS_tdR = -1*MPS_tdR, #P2
         MPS_clR = -1*MPS_clR, #P3
         HDWK = -1*HDWK, MPS_acR = -1*MPS_acR
  ) %>%
  mutate_at(vars(-subject), list(~as.vector(scale(.)))) %>%
  dplyr::rename(id=subject)

for_msem <- flanker %>% inner_join(tofactor, by="id")
anti_join(tofactor, flanker, by="id")

m15_mplus <- mplusObject(
  TITLE="Antagonism and Constraint BSEM CFA 3+3 Parcels, Flanker ACC",
  DEFINE="
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
    cond = cond - 1; ! 0=congruent, 1=incongruent
    block = block - 1; ! 0=most_incon, 1=most_con
  ",
  VARIABLE = "
    WITHIN = cond block prev_rt_z prev_error;
    USEVARIABLES = id correct cond block prev_rt_z prev_error
          MANIP p1 p2 p3 p4 p5;
    BETWEEN = MANIP p1 p2 p3 p4 p5;
    CLUSTER = id;
    CATEGORICAL = correct;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(20000) 100000;
    CHAINS=2;
    PROCESSORS=4;
  ",
  MODEL = "
    %WITHIN%
    correct ON block;
    condslo | correct ON cond;        !random slope of cond
    prevslo | correct ON prev_rt_z;   !random slope of rt autocorrelation
    errslo | correct ON prev_error;   !random slope of previous error

    %BETWEEN%
    !means of random slopes
    [condslo];
    [prevslo];
    [errslo];

    !variances of random slopes
    condslo;
    prevslo;
    errslo;

    !slope correlations
    condslo prevslo errslo WITH
       condslo prevslo errslo;

    [correct$1]; !intercept: average correct

    !TRAIT MODEL

    antag BY p5*1 ! MPS_agR AGG
      p1        ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    constraint BY
      p2*  !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p4; !MPS_acR HDWK;
    constraint@1;

    antag WITH constraint;

    !TRAIT MODERATION

    ! trait moderates flanker performance -- incongruency effect and stickiness of RT
    condslo ON antag constraint;
    prevslo ON antag constraint;
    errslo ON antag constraint;

    correct ON antag constraint;
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem
)

mout <- mplusModeler(m15_mplus, dataout="flanker_acc_mplus.dat",
                     modelout="acc_m15_AC.inp", run=TRUE, hashfilename = FALSE,
                     Mplus_command = mpluscmd)
beepr::beep()
summary(mout)


flextable(mout$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit() #%>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")





# AC MSEM with cross-loading

# Between-persons CFA suggests that P1 cross-loads on constraint. Including this cross-loading makes the fit good enough IMO, while without it, CFI is < .9. Do our effects here hold if we allow that?
  
  
m15_mplus <- mplusObject(
  TITLE="Antagonism and Constraint BSEM CFA 3+3 Parcels, Flanker ACC Mod",
  DEFINE="
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
    cond = cond - 1; ! 0=congruent, 1=incongruent
    block = block - 1; ! 0=most_incon, 1=most_con
  ",
  VARIABLE = "
    WITHIN = cond block prev_rt_z prev_error;
    USEVARIABLES = id correct cond block prev_rt_z prev_error
          MANIP p1 p2 p3 p4 p5;
    BETWEEN = MANIP p1 p2 p3 p4 p5;
    CLUSTER = id;
    CATEGORICAL = correct;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(10000);
    CHAINS=2;
    PROCESSORS=4;
  ",
  MODEL = "
    %WITHIN%
    correct ON block;
    condslo | correct ON cond;        !random slope of cond
    prevslo | correct ON prev_rt_z;   !random slope of rt autocorrelation
    errslo | correct ON prev_error;   !random slope of previous error

    %BETWEEN%
    !means of random slopes
    [condslo];
    [prevslo];
    [errslo];

    !variances of random slopes
    condslo;
    prevslo;
    errslo;

    !slope correlations
    condslo prevslo errslo WITH
       condslo prevslo errslo;

    [correct$1]; !intercept: average correct

    !TRAIT MODEL

    antag BY p5*1 ! MPS_agR AGG
      p1        ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    constraint BY
      p2*  !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p1        ! ADD IN MISTRUST MPS_alR
      p4; !MPS_acR HDWK;
    constraint@1;

    antag WITH constraint;

    constraint BY P1; !modification index cross-loading for fit

    !TRAIT MODERATION

    ! trait moderates flanker performance -- incongruency effect and stickiness of RT
    condslo ON antag constraint;
    prevslo ON antag constraint;
    errslo ON antag constraint;

    correct ON antag constraint;
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem
)

mout <- mplusModeler(m15_mplus, dataout="flanker_acc_mplus.dat",
                     modelout="acc_m15_AC_mod.inp", run=TRUE, hashfilename = FALSE,
                     Mplus_command = mpluscmd)

summary(mout)
# kable(mout$results$parameters$stdyx.standardized %>%
#         filter(!paramHeader %in% c("Intercepts", "Means", "Residual Variances")))
# 


flextable(mout$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit() #%>% save_as_docx(path = "~/Desktop/flanker_acc_traits.docx")


