
#basedir <- "C:/Users/timot/Documents/GitHub/PD_Inhibition_DDM"
# basedir <- "~/Data_Analysis/PD_Inhibition_DDM" 
basedir <- "~/github_repos/PD_Inhibition_DDM" 

# setwd("~/github_repos/PD_Inhibition_DDM/Code/mlm/flanker/")

if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(car, nlme, lme4, readr, tidyverse, emmeans, cowplot, MplusAutomation, knitr, sas7bdat, ggcorrplot, psych, flextable)
# knitr::opts_chunk$set(echo = TRUE) #print code by default
# #knitr::opts_chunk$set(cache.path = "../../Outputs/flanker_mlms/")
# options(digits=3)
# options(width=180)

#flanker data processed in flanker_mlm.Rmd
#has the worst subjects removed, rt_inv transformation, etc.
load(paste0(basedir,"/Data/preprocessed/flanker_for_traits.RData"))
# flanker <- 
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



# setup self-reports  -----------------------------------------------------


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



# run model ---------------------------------------------------------------
# AC MSEM with cross-loading

# Between-persons CFA suggests that P1 cross-loads on constraint. Including this cross-loading makes the fit good enough IMO, while without it, CFI is < .9. Do our effects here hold if we allow that?
#   
#   Indeed, the ANTAG effects are consistent, even slightly stronger. So, no problems with whether P1 is a secondary indicator of constraint. But, the stickiness of RT `PREVSLO ON constraint` becomes non-significant, p = .10.

tictoc::tic()
AC_msem_mod <- mplusObject(
  TITLE="Antagonism and Constraint BSEM CFA 3+3 Parcels, Flanker RT Mod",
  DEFINE="
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
    cond = cond - 1; ! 0=congruent, 1=incongruent
    block = block - 1; ! 0=most_incon, 1=most_con
    cond_block = cond*block; ! 1=most_con, incongruent
  ",
  
  VARIABLE = "
    WITHIN = cond block cond_block prev_rt_z prev_error trial_z run_trial_z;
    USEVARIABLES = id rt_inv cond block prev_rt_z 
      prev_error trial_z run_trial_z
      MANIP p1 p2 p3 p4 p5 cond_block;
    BETWEEN = MANIP p1 p2 p3 p4 p5;
    CLUSTER = id;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(15000);
    CHAINS=2;
    PROCESSORS=4;
  ",
  MODEL = "
  
  %WITHIN%
    rt_inv ON block (b_block)
      cond_block (b_cxb);
    condslo | rt_inv ON cond;        !random slope of cond
    prevslo | rt_inv ON prev_rt_z;   !random slope of rt autocorrelation
    errslo | rt_inv ON prev_error;   !random slope of previous error
    runslo | rt_inv ON run_trial_z;  !random slope of run_trial
    trialslo | rt_inv ON trial_z;    !random slope of trial

    
  %BETWEEN%
    !means of random slopes
    [condslo] (c_b0);
    [runslo];
    [trialslo];
    [prevslo] (p_b0);
    [errslo] (e_b0);
    
    !variances of random slopes
    condslo;
    runslo;
    trialslo;
    prevslo;
    errslo;
    
    !slope correlations
    condslo runslo trialslo prevslo errslo WITH
       condslo runslo trialslo prevslo errslo;
    
    [rt_inv] (b0); !mean average inverse RT
    
    !TRAIT MODEL
        
    antag BY p5* ! MPS_agR AGG
      p1        ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    constraint BY
      p2*  !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p4; !MPS_acR HDWK;
    constraint@1;

    antag WITH constraint;
    
    constraint BY P1; !modification index cross-loading for fit

    !TRAIT MODERATION
    
    ! trait moderates flanker performance
    condslo ON antag (c_bA)
      constraint (c_bC);
    
    prevslo ON antag (p_bA)
      constraint (p_bC);
    
    errslo ON antag (e_bA) 
      constraint (e_bC);
    
    !N.B. Leaving out the association of antag with rt_inv 
    !omits a hugely important relationship.
    !Thus, allow antag as a predictor of average (person) RT
    rt_inv ON antag (bA) 
      constraint (bC);
  ",
  MODELCONSTRAINT = "
    NEW (RT_I_loA RT_I_miA RT_I_hiA
         RT_I_loC RT_I_miC RT_I_hiC
         RT_C_loA RT_C_miA RT_C_hiA
         RT_C_loC RT_C_miC RT_C_hiC
         tLo tMi tHi
      );
      
    tLo = -1;
    tMi = 0;
    tHi = 1;
    
    !since trial_z and prev_rt_z are mean centered, they
    !drop out of the predicted values (marginalized out)
    !these are predicted values at prev_error = 0, which is more
    !representative anyhow (only 3% errors)
    !Use block = 0.5 to compute the average across blocks
    
    !pred RT for incongruent trials (cond=1), low Antag
    RT_I_loA = b0 + c_b0 + bA*tLo + c_bA*tLo + 0.5*b_block + 0.5*b_cxb;
    RT_I_miA = b0 + c_b0 + bA*tMi + c_bA*tMi + 0.5*b_block + 0.5*b_cxb;
    RT_I_hiA = b0 + c_b0 + bA*tHi + c_bA*tHi + 0.5*b_block + 0.5*b_cxb;
    
    !pred RT for congruent trials (cond=0), low Antag
    RT_C_loA = b0 + bA*tLo + 0.5*b_block;
    RT_C_miA = b0 + bA*tMi + 0.5*b_block;
    RT_C_hiA = b0 + bA*tHi + 0.5*b_block;
    
    !pred RT for incongruent trials (cond=1), low Constraint
    RT_I_loC = b0 + c_b0 + bC*tLo + c_bC*tLo + 0.5*b_block + 0.5*b_cxb;
    RT_I_miC = b0 + c_b0 + bC*tMi + c_bC*tMi + 0.5*b_block + 0.5*b_cxb;
    RT_I_hiC = b0 + c_b0 + bC*tHi + c_bC*tHi + 0.5*b_block + 0.5*b_cxb;
    
    !pred RT for congruent trials (cond=0), low Antag
    RT_C_loC = b0 + bC*tLo + 0.5*b_block;
    RT_C_miC = b0 + bC*tMi + 0.5*b_block;
    RT_C_hiC = b0 + bC*tHi + 0.5*b_block;
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT ="
    TECH1 TECH8 STANDARDIZED CINTERVAL;
  ",
  rdata=for_msem
)

mres <- mplusModeler(AC_msem_mod, modelout = file.path("rt_m17_AC_mod.inp"), 
                     run=TRUE, hashfilename = FALSE,
                     Mplus_command = mpluscmd)

tictoc::toc()

summary(mres)
kable(mres$results$parameters$stdyx.standardized %>% 
        filter(!paramHeader %in% c("Intercepts", "Means", "Residual Variances")))




flextable(mres$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit() %>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")





# run model ---------------------------------------------------------------


