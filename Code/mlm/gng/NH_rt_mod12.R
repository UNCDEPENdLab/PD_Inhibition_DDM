setwd("/Users/michael/Data_Analysis/PD_Inhibition_DDM/Code/mlm/gng")
if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(car, brms, nlme, lme4, loo, readr, tidyverse, emmeans, cowplot, glmmTMB, bbmle, broom, MplusAutomation)


# this is an admixture of Alison and Nate's code.

# Load data
# Data has already been cleaned and initial analyses completed by NH. 
# For RT analyses, use variable rt_log_trim_grp
# For outcome analyses, use response -> this indicates whether they correctly Go vs No Go'd (can figure out whether they spressed space or not based on the stim)
# Cond refers to number of prior go trials 
# Stim refers whether go or no go trial
# Use trial_z rather than trial (for scaling purposes)
# Will likely add in prev_rt in more complex analyses

#basedir <- "~/github_repos/PD_Inhibition_DDM"
basedir <- "~/Data_Analysis/PD_Inhibition_DDM"

#NH add from GNG RT analyses
gng <- read.csv(file.path(basedir, "Data/preprocessed/go_nogo_full_accCode.csv")) %>% mutate(bad_subj = if_else(subj_idx == 15, 1, 0)) #flagged 15 as bad as this was the one person who needed to be according to NH's preprocessing; will test how inclusion of 15 alters results
gng <- mutate(gng, cond= factor(cond, levels = c("OneGo","ThreeGo", "FiveGo", "SevenGo")))

gng <- gng %>% group_by(subj_idx) %>% mutate(
  prev_rt = ifelse(lag(rt_log_trim_grp) == 0, NA, lag(rt_log_trim_grp)),
  prev_rt_z = scale(prev_rt),
  prev_error = ifelse(lag(response) == 1, 0, 1),
  id = as.character(subj_idx))

# gng_rt <- dplyr::filter(gng, stim == "Go", response == 1) #only analyze GO RTs!
gng$cond_id <- with(gng, interaction(cond,id))
gng <- mutate(gng, block_trial_z = as.vector(scale(block_trial))) #%>% filter(response == 1, stim == "Go")

# gng <- mutate(gng, block_trial_z = as.vector(scale(block_trial))) 

#pull self reps

# fscores <- read.csv(paste0(basedir,"/Outputs/factor_structure/PDDDM_factor_scores_5-18-20.csv")) %>% rename(subj_idx = subject) 

fscores_update <- read.csv(file.path(basedir, "Code/factor_structure/AC_bsem_cfa_parcel_mod_savedata.csv")) %>% select(SUBJECT, DISINHIB.Median, ANTAG.Median) %>% rename(subj_idx = SUBJECT, DISINH = DISINHIB.Median, ANTAG = ANTAG.Median) #%>% mutate(DISINH = -1*DISINH) #reverse code Constraint to set in maladaptive direction. 
cor(fscores_update) #make sure scores positively correlate

gng <- gng %>% left_join(fscores_update, by = "subj_idx")

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
#mpluscmd <- "/Users/natehall/Desktop/Mplus/mplus"
mpluscmd <- "/Applications/Mplus/mplus"

# RT analyses
# block trial*trial and prev_rtXblock_trial effects on RT. homogonous covariance structure estimated per subject and variance estimated per cond. Random intercept and slope of trial and block_trial estmimated for each subject. Random slope of block trial estimated per condition within subject.


# Winning from NH updates: Drop nesting of condition for simplicity, as in flanker_acc_traits


# minfo[["m12"]] <- c(fixed="Trials from no-go, Trial number, Trials from no-go x Trial number", l2="Condition within subject: Trials from no-go\nSubject: Intercept, Trials from no-go, Trial")
# 
# m12 <- lmer(rt_log_trim_grp ~  block_trial_z*trial_z + (1 + block_trial_z + trial_z| id) + (1 + block_trial_z | id:cond),
#             na.action = na.exclude, 
#             data = gng_rt, control=lmerControl(optCtrl=list(maxfun=2e4)), REML = FALSE) #bump up iterations)
# summary(m12)

tofactor <- self_reps %>% select(subject, MPS_agR, AGG, MISTRUST, MPS_alR, MANIP,
                                 PROPER, MPS_clR, IMPUL, MPS_tdR, MPS_acR, HDWK, K10, STAI) %>% 
  # mutate(IMPUL =-1*IMPUL) %>% #score toward constraint to make loadings upright. Update 6/25: reverse coding to scale positively with high externalizing
  mutate(PROPER = -1*PROPER, MPS_tdR = -1*MPS_tdR, #P2
         MPS_clR = -1*MPS_clR, #P3
         HDWK = -1*HDWK, MPS_acR = -1*MPS_acR
  ) %>%
  mutate_at(vars(-subject), list(~as.vector(scale(.)))) %>% 
  dplyr::rename(id=subject)

ggcorrplot::ggcorrplot(cor(tofactor %>% select(-id), use="pairwise")) #looks right

for_msem <- gng %>% mutate(id = subj_idx) %>% inner_join(tofactor, by="id") %>% ungroup() %>%  
  select(id, response, cond, stim, trial_z, rt_log_trim_grp, block_trial_z,
         prev_error, MPS_agR, AGG, MISTRUST, MPS_alR, MANIP,
         PROPER, MPS_clR, IMPUL, MPS_tdR, MPS_acR, HDWK, K10, STAI) %>% 
  mutate(cond = ifelse(cond == "OneGo", 0, 
                       ifelse(cond == "ThreeGo", 1, 
                              ifelse(cond == "FiveGo", 2,3))),
         stim = ifelse(stim == "Go", 0,1))

head(for_msem)

for_msem_rt <- for_msem %>% dplyr::filter(stim == 0 & response == 1) #only analyze GO RTs!
lattice::histogram(~rt_log_trim_grp, for_msem_rt) #looks normal to me
xtabs(~id, for_msem_rt) #all similar

lattice::histogram(~block_trial_z, for_msem_rt) #okay
lattice::histogram(~trial_z, for_msem_rt) #relatively uniform given that trial unfolds linearly

# m12 <- lmer(rt_log_trim_grp ~ block_trial_z*trial_z + (1 + block_trial_z + trial_z| id) + (1 + block_trial_z | id:cond),

## MNH: Translate m12 to Mplus (no traits)
rt_with_mod_gng <- mplusObject(
  TITLE = "GnG M12 Mplus",
  DEFINE = "
    t_ixn = block_trial_z*trial_z; !ixn of trials from no-go and overall trial;
  ",
  VARIABLE = "
    WITHIN = trial_z block_trial_z t_ixn;
    USEVARIABLES = id rt_log_trim_grp block_trial_z trial_z t_ixn;
    CLUSTER = id;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(15000);
    CHAINS=4;
    PROCESSORS=4;
  ",
  MODEL = "
    %WITHIN%

    b_trial | rt_log_trim_grp ON trial_z;
    b_block_trial |  rt_log_trim_grp ON block_trial_z;
    rt_log_trim_grp ON t_ixn;

    %BETWEEN%
    [b_trial b_block_trial]; !slope means
    b_trial b_block_trial; !slope variances

    !slope correlations
    b_trial b_block_trial WITH 
      b_trial b_block_trial;

    [rt_log_trim_grp] (b0); !mean average log RT
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem_rt
)
mout <- mplusModeler(rt_with_mod_gng, modelout="gng_m12_mplus.inp", 
                     run=TRUE, Mplus_command = mpluscmd, hashfilename = FALSE)


# attempt 1, not looking too great. ---------------------------------------


# rt_with_mod_gng <- mplusObject(
#   TITLE = "Antagonism, Disinhibition moderates gng random slopes and fixed effects interaction",
#   DEFINE = "
#     t_ixn = block_trial_z*trial_z; !ixn of trials from no-go and overall trial;
#     p1=MEAN(MISTRUST MPS_alR);
#     p2=MEAN(PROPER MPS_tdR);
#     p3=MEAN(MPS_clR IMPUL);
#     p4=MEAN(MPS_acR HDWK);
#     p5=MEAN(MPS_agR AGG);
#   ",
#     
#   VARIABLE = "
#     WITHIN = trial_z block_trial_z;
#     USEVARIABLES = id rt_log_trim_grp block_trial_z trial_z 
#       MANIP p1 p2 p3 p4 p5
#       t_ixn;
#     BETWEEN = MANIP p1 p2 p3 p4 p5;
#     
#     CLUSTER = id;
#   ",
#   ANALYSIS = "
#     TYPE=TWOLEVEL RANDOM;
#     ESTIMATOR=BAYES;
#     BITERATIONS=(15000);
#     BCONVERGENCE=.02;
#     CHAINS=2;
#     PROCESSORS=2;
#   ",
#   
# 
#   MODEL = "
#     %WITHIN%
#     
#     b_trial | rt_log_trim_grp ON trial_z; 
#     b_block_trial |  rt_log_trim_grp ON block_trial_z; 
#     b_trial_ixn | rt_log_trim_grp ON t_ixn;
#       
#     
#     %BETWEEN%
#     [b_trial b_block_trial b_trial_ixn]; !slope means
#     b_trial b_block_trial b_trial_ixn; !slope variances
#     
#         !slope correlations
#     b_trial b_block_trial b_trial_ixn WITH
#        b_trial b_block_trial b_trial_ixn;
#     
#     [rt_log_trim_grp] (b0); !mean average log RT
#     ![t_ixn] (bixn); !mean ixn term at between subjects level
#     
#     
#     !TRAIT MODEL
#         
#     antag BY 
#       p5* ! MPS_agR AGG
#       p1        ! MISTRUST MPS_alR
#       MANIP;
#     antag@1;
# 
#     disinhib BY
#       p2*  !PROPER MPS_tdR
#       p3  !MPS_clR IMPUL
#       p4 !MPS_acR HDWK;
#       p1; 
#     disinhib@1;
# 
#     antag WITH disinhib;
#     
#     !disinhib BY p1; !modification index cross-loading for fit
# 
#     !TRAIT MODERATION
#     
#     ! trait moderates flanker performance
#     b_block_trial ON antag (b_bA)
#       disinhib (b_bD);
#     
#     b_trial ON antag (t_bA)
#       disinhib (t_bD);
#     
#     b_trial_ixn ON antag (ixn_bA) 
#       disinhib (ixn_bD);
#     
#     !N.B. Leaving out the association of antag with rt_inv 
#     !omits a hugely important relationship.
#     !Thus, allow antag as a predictor of average (person) RT
#     
#     rt_log_trim_grp ON antag (bA) 
#       disinhib (bD);
#     
#   ",
#   PLOT = "TYPE = PLOT2;",
#   OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
#   rdata = for_msem_rt
# )
# mout <- mplusModeler(rt_with_mod_gng,
#                      modelout="rt_with_mod_gng.inp", run=TRUE, Mplus_command = mpluscmd, hashfilename = FALSE)
# 
# 
# flextable(mout$results$parameters$stdyx.standardized %>% 
#             filter(!paramHeader %in% c("Variances")) %>%
#             select(-sig, -posterior_sd) %>%
#             mutate(pval=2*pval, #convert to 2-tailed
#                    pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
# ) %>% autofit()# %>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")
# 


# wrong to fit ixn as random ----------------------------------------------

rt_with_mod_gng2 <- mplusObject(
  TITLE = "Antagonism, Disinhibition moderates gng random slopes and fixed effects interaction",
  DEFINE = "
    t_ixn = block_trial_z*trial_z; !ixn of trials from no-go and overall trial;
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
    CENTER trial_z block_trial_z t_ixn (GRANDMEAN); !make intercepts easy to understand
  ",
  
  VARIABLE = "
    WITHIN = trial_z block_trial_z t_ixn;
    USEVARIABLES = id rt_log_trim_grp block_trial_z trial_z 
      MANIP p1 p2 p3 p4 p5 t_ixn;
    BETWEEN = MANIP p1 p2 p3 p4 p5;
    CLUSTER = id;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(30000);
    CHAINS=4;
    PROCESSORS=4;
  ",
  
  MODEL = "
    %WITHIN%

    b_trial | rt_log_trim_grp ON trial_z;
    b_block_trial |  rt_log_trim_grp ON block_trial_z;
    rt_log_trim_grp ON t_ixn;

    %BETWEEN%
    [b_trial b_block_trial]; !slope means
    b_trial b_block_trial; !slope variances

    !slope correlations
    b_trial b_block_trial WITH 
      b_trial b_block_trial;

    [rt_log_trim_grp] (b0); !mean average log RT
    
    !TRAIT MODEL
        
    antag BY 
      p5* ! MPS_agR AGG
      p1  ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    disinhib BY
      p2* !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p4  !MPS_acR HDWK
      p1; !cross-load
    disinhib@1;

    antag WITH disinhib;
    
    !TRAIT MODERATION
    
    ! trait moderates flanker performance
    b_block_trial ON antag (b_bA)
      disinhib (b_bD);
    
    b_trial ON antag (t_bA)
      disinhib (t_bD);
    
    !this would only make sense if we had meaningful between-person variation in
    !the interaction and were decomposing person emans of the interaction for analysis.
    !or, we could add a random slope of the interaction and model that here b_ixn ON ...
    !but as written, this doesn't make sense.
    
    !t_ixn ON antag (ixn_bA) 
    !  disinhib (ixn_bD);

    !average RT on traits    
    rt_log_trim_grp ON antag (bA) 
      disinhib (bD);
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem_rt
)
mout <- mplusModeler(rt_with_mod_gng2, modelout="rt_with_mod_gng2.inp", 
                     run=TRUE, Mplus_command = mpluscmd, hashfilename = FALSE)


flextable(mout$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit()# %>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")


# try dropping ixn? -------------------------------------------------------


rt_with_mod_gng3 <- mplusObject(
  TITLE = "Antagonism, Disinhibition moderates gng random slopes and fixed effects interaction",
  DEFINE = "
    !t_ixn = block_trial_z*trial_z; !ixn of trials from no-go and overall trial;
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
  ",
  
  VARIABLE = "
    WITHIN = trial_z block_trial_z;
    USEVARIABLES = id rt_log_trim_grp block_trial_z trial_z 
      MANIP p1 p2 p3 p4 p5;
      !t_ixn;
    BETWEEN = MANIP p1 p2 p3 p4 p5;
    
    CLUSTER = id;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(15000);
    BCONVERGENCE=.02;
    CHAINS=2;
    PROCESSORS=2;
  ",
  
  
  MODEL = "
    %WITHIN%
    
    b_trial | rt_log_trim_grp ON trial_z; 
    b_block_trial |  rt_log_trim_grp ON block_trial_z; 
      
    
    %BETWEEN%
    [b_trial b_block_trial]; !slope means
    b_trial b_block_trial; !slope variances
    
        !slope correlations
    b_trial b_block_trial  WITH
       b_trial b_block_trial;
    
    [rt_log_trim_grp] (b0); !mean average log RT
    
    
    !TRAIT MODEL
        
    antag BY 
      p5* ! MPS_agR AGG
      p1        ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    disinhib BY
      p2*  !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p4 !MPS_acR HDWK;
      p1; 
    disinhib@1;

    antag WITH disinhib;
    
    !disinhib BY p1; !modification index cross-loading for fit

    !TRAIT MODERATION
    
    ! trait moderates flanker performance
    b_block_trial ON antag (b_bA)
      disinhib (b_bD);
    
    b_trial ON antag (t_bA)
      disinhib (t_bD);
    
    !t_ixn ON antag (ixn_bA) 
    !  disinhib (ixn_bD);
    
    !N.B. Leaving out the association of antag with rt_inv 
    !omits a hugely important relationship.
    !Thus, allow antag as a predictor of average (person) RT
    
    rt_log_trim_grp ON antag (bA) 
      disinhib (bD);
    
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem
)
mout <- mplusModeler(rt_with_mod_gng3,
                     modelout="rt_with_mod_gng3.inp", run=TRUE, Mplus_command = mpluscmd, hashfilename = FALSE)


flextable(mout$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit()# %>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")




# ACCURACY ANALYSES -------------------------------------------------------

#winning model: m19

# fixed: stimulus, condition, trial, trialxstimulus
# random: stimulus


acc_with_mod_gng <- mplusObject(
  TITLE = "Antagonism, Disinhibition GNG ACC",
  DEFINE = "
    ts_ixn = stim*trial_z; !ixn of overall trial and stimulus;
    p1=MEAN(MISTRUST MPS_alR);
    p2=MEAN(PROPER MPS_tdR);
    p3=MEAN(MPS_clR IMPUL);
    p4=MEAN(MPS_acR HDWK);
    p5=MEAN(MPS_agR AGG);
  ",
  
  VARIABLE = "
    ! WITHIN = stim; !if my understanding is correct, this will only specify an effect at the within-level
    USEVARIABLES = id response stim cond trial_z 
      MANIP p1 p2 p3 p4 p5
      ts_ixn;
    BETWEEN = MANIP p1 p2 p3 p4 p5
      cond trial_z
      ts_ixn;
    
    CLUSTER = id;
  ",
  ANALYSIS = "
    TYPE=TWOLEVEL RANDOM;
    ESTIMATOR=BAYES;
    BITERATIONS=(15000);
    BCONVERGENCE=.02;
    CHAINS=2;
    PROCESSORS=2;
  ",
  
  
  MODEL = "
    %WITHIN%
    
    b_stim | response ON stim; 
    b_cond | response ON cond;
    b_trial | response ON trial_z
    b_ts_ixn | response ON ts_ixn
    
    %BETWEEN%
    [b_stim b_cond b_trial b_ts_ixn]; !slope means
    b_stim b_cond b_trial b_ts_ixn; !slope variances
    
        
    [response] (b0); !mean average ACC
    
    
    !TRAIT MODEL
        
    antag BY 
      p5* ! MPS_agR AGG
      p1        ! MISTRUST MPS_alR
      MANIP;
    antag@1;

    disinhib BY
      p2*  !PROPER MPS_tdR
      p3  !MPS_clR IMPUL
      p4 !MPS_acR HDWK;
      p1; 
    disinhib@1;

    antag WITH disinhib;
    
    !disinhib BY p1; !modification index cross-loading for fit. Included above.

    !TRAIT MODERATION
    
    ! trait moderates flanker performance
    b_stim ON antag (b_bA)
      disinhib (b_bD);
    
    
    ts_ixn ON antag (ixn_bA) 
      disinhib (ixn_bD);
    
    cond ON antag  
      disinhib; 
    
    !N.B. Leaving out the association of antag with rt_inv 
    !omits a hugely important relationship.
    !Thus, allow antag as a predictor of average (person) RT
    
    response ON antag (bA) 
      disinhib (bD);
    
  ",
  PLOT = "TYPE = PLOT2;",
  OUTPUT = "TECH1 TECH8 STANDARDIZED CINTERVAL;",
  rdata = for_msem
)
mout <- mplusModeler(acc_with_mod_gng,
                     modelout="acc_with_mod_gng.inp", run=TRUE, Mplus_command = mpluscmd, hashfilename = FALSE)


flextable(mout$results$parameters$stdyx.standardized %>% 
            filter(!paramHeader %in% c("Variances")) %>%
            select(-sig, -posterior_sd) %>%
            mutate(pval=2*pval, #convert to 2-tailed
                   pval=if_else(pval < .05, as.character(paste0(pval, "*")), as.character(pval))) 
) %>% autofit()# %>% save_as_docx(path = "~/Desktop/flanker_rt_traits.docx")




