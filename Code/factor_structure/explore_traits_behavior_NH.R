
# 5/13/20 NH combine flanker, rp model params, behavior, and traits -------

pacman::p_load(tidyverse, psych, ggcorrplot, lavaan, GPArotation)

# basedir <- "C:/Users/timot/Documents/GitHub/PD_Inhibition_DDM" # Tim
basedir <- "~/github_repos/PD_Inhibition_DDM" # Nate

setwd(basedir)

# alldat <- read_csv("/Data/alldat.csv")
mpq <- get(load(file.path(basedir, "Data/preprocessed/MPQ_all_scored_final.RData"))) %>% filter(exclude_MPQ == 0)
snap <- get(load(file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData"))) %>% filter(exclude_SNAP == 0)
alldat <- mpq %>% inner_join(snap, by = "subject") %>% tibble()

names(alldat)


# trim to Domains and Traits ----------------------------------------------



all_pers <- alldat %>% select(subject, ends_with("T",ignore.case = FALSE), starts_with("T_")) %>% select(-MISTRUST, -LOSLFEST, -T_VRIN, -T_TRIN, -T_DRIN, -T_RareVirtues, -T_Deviance, -T_InvalidityIndex, -T_BackDeviance)
names(all_pers)

pDomains <- all_pers %>% select(subject, MPS_PET, MPS_NET, MPS_COT, MPS_PGT, MPS_PCT, MPS_NGT, MPS_NLT, T_NegativeTemperament, T_PositiveTemperament, T_Disinhibition) %>% rename(subj_idx = subject)
tnames <- names(all_pers)[!names(all_pers) %in% names(pDomains)]
pTraits <- all_pers %>% select("subject", all_of(tnames)) %>% rename(subj_idx = subject)


# additional post-hoc trimming --------------------------------------------

# not interested (psychoticism) or messy scale
Tsr_fa <- pTraits %>% select(-subj_idx) %>% select(-T_EccentricPerceptions, -MPS_abT, -T_SelfHarm, -MPS_vrT, -MPS_uvT, -MPS_trT) #-T_SelfHarm)
Dsr_fa <- pDomains %>% select(-subj_idx, -MPS_PGT, -MPS_PCT, -MPS_NGT, -MPS_NLT) # Don't know what a number of the MPQ vars mean, but SNAP could be of interest. Update: just keep 3 domains from Tim's spreadsheet. These are what we care about.


# factor analyze just self-reports (traits) ----------------------------------------

map <- nfactors(Tsr_fa, n=10, rotate="oblimin", fm="ml")
map
par <- fa.parallel(sr_fa, fm = "ml", fa = "fa")

#store factor scores based on all tested solutions for visualization below
f_scores <- list()

# 9 factors ---------------------------------------------------------------

sr.fa <- fa(Tsr_fa, 9, rotate="oblimin", fm = "ml") #corrplot looks like maybe we could get 9 factors? 
print(sr.fa, cut=.3, sort = T) ##looks interpretable to me
#pa estimation
sr.fa <- fa(Tsr_fa, 9, rotate="oblimin", fm = "pa") 
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 

scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>% 
  rename(antag_machiavel = PA9,
         inverse_affiliativeE = PA8,
         achieve = PA5,
         extraversion_narcissism = PA3,
         impulsivity = PA2,
         rule_following = PA6,
         suspiciousness = PA1,
         anxious_sensitive = PA4,
         inverse_self_esteem = PA7
  )
f_scores[["factor_9"]] <- scores
# 8 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 8, rotate="oblimin", fm = "ml") #ml
# print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 
sr.fa <- fa(Tsr_fa, 8, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antag_machiavel = PA7,
         inverse_affiliativeE = PA8,
         extraversion_narcissism = PA3,
         impulsivity = PA2,
         achieve = PA5,
         rule_following = PA6,
         suspiciousness = PA1,
         neuroticism = PA4
  )
f_scores[["factor_8"]] <- scores
# 7 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 7, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 7, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T)
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antag_machiavel = PA7,
         inverse_affiliativeE = PA1,
         extraversion_narcissism = PA3,
         impulsivity = PA2,
         achieve = PA5,
         anxiousness = PA4,
         rule_following = PA6
  )
f_scores[["factor_7"]] <- scores
# 6 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 6, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 6, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antagonism = PA1,
         inverse_affiliativeE = PA3,
         impulsivity = PA2,
         achieve = PA5,
         neuroticism = PA4,
         conscientiousness = PA6
  )
f_scores[["factor_6"]] <- scores
# 5 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 5, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 5, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antag_machiavel = PA1,
         extraversion_narcissism = PA3,
         control = PA2,
         achieve = PA5,
         neuroticism = PA4
  )
f_scores[["factor_5"]] <- scores
# 4 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 4, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 4, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antag_externalizing = PA1,
         control = PA2,
         extraversion_narcissism = PA3,
         neuroticism_anx = PA4
  )
f_scores[["factor_4"]] <- scores
# 3 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 3, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 3, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(antag_externalizing = PA1,
         control = PA2,
         extraversion_narcissism = PA3
  )
f_scores[["factor_3"]] <- scores

# 2 factors ---------------------------------------------------------------
# sr.fa <- fa(Tsr_fa, 2, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(Tsr_fa, 2, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
scores <- data.frame(subj_idx = pTraits$subj_idx, factor.scores(Tsr_fa, sr.fa)$scores) %>%
  rename(inverse_closeness = PA1,
         impulsivity_externalizing = PA2
  )
f_scores[["factor_2"]] <- scores

# gather all relevant task behavior ---------------------------------------

R.utils::sourceDirectory("Code/Functions") # source all of Nate's functions for the project (won't use all of them)

behav <- summarise_task_behavior(basedir) %>% ungroup()


# pull traces -------------------------------------------------------------

source(paste0(basedir, "Code/Nate/temp_pull_posterior_summaries.R"))


#loop over task and plot corrs of task behavior with traits, dimensions, and factor scores
tasks <- c("flanker", "recent_probes", "go_nogo", "all")
correlations <- list() # to help guide towards large correlations

behav <- behav %>% left_join(ddm_params, by = "subj_idx")

for(t in tasks){
  pdf(file = paste0("Figures/", t, "_personality_corrs.pdf"), width = 11, height = 8)
  # t <- "flanker"
  if(t == "all"){
    tdf <- behav 
  } else{
    tdf <- behav %>% select(subj_idx, starts_with(t))  
  }
  
  
  # start with broad dimensions from SNAP and MPQ
  domdf <- pDomains %>% left_join(tdf, by = "subj_idx")
  f_cors <- cor(domdf, use = "pairwise.complete.obs")[-1,-1] # drop subj_idx
  p <- ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1, title = paste0("bivariate correlations between ", t, " behavior and broad personality domains"))
  plot(p)
  correlations[[t]][["Domains"]] <- f_cors
  
  # move down to traits
  traitdf <- pTraits %>% left_join(tdf, by = "subj_idx")
  f_cors <- cor(traitdf, use = "pairwise.complete.obs")[-1,-1] # drop subj_idx
  p <- ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1, title = paste0("bivariate correlations between ", t, " behavior and personality traits"))
  plot(p)
  correlations[[t]][["Traits"]] <- f_cors
  
  ## loop over factor solutions
  for(f in names(f_scores)){
    # f <- "factor_3"
    fdf <- f_scores[[f]]
    facdf <-  fdf %>% left_join(tdf, by = "subj_idx")
    f_cors <- cor(facdf, use = "pairwise.complete.obs")[-1,-1] # drop subj_idx
    p <- ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1, title = paste0("bivariate correlations between ", t, " behavior and factor scores based on a ", f, " solution"))
    plot(p)
    correlations[[t]][[f]] <- f_cors
  }
  dev.off()
}

save(correlations, file = "Outputs/behavior_traits_allcorrelations.RData")



# just look at behavioral indices -----------------------------------------
zero_sds <- c()
for(i in 1:ncol(behav)){
  nam <- names(behav)[i]
  if(sd(pull(behav, nam), na.rm = TRUE) == 0) {
    print(nam)
    zero_sds <- c(zero_sds, nam)
    }
}

drop_zsd <- select(behav, -zero_sds, -subj_idx)

fa_ddm <- select(behav, contains("MAP"))

map <- nfactors(fa_ddm, rotate = "oblimin", n=7,  fm="fiml")
map
par <- fa.parallel(fa_ddm, fm = "fiml", fa = "fa")

sr.fa <- fa(fa_ddm, 2, rotate="oblimin", fm = "fiml") #corrplot looks like maybe we could get 9 factors? 
print(sr.fa, cut=.3, sort = T) ##looks interpretable to me
#pa estimation
sr.fa <- fa(fa_ddm, 2, rotate="oblimin", fm = "pa") 
print(sr.fa, cut=.3, sort = T) 

sr.fa <- fa(fa_ddm, 3, rotate="oblimin", fm = "ml") 
print(sr.fa, cut=.3, sort = T) 

sr.fa <- fa(fa_ddm, 4, rotate="oblimin", fm = "ml") 
print(sr.fa, cut=.3, sort = T) 

sr.fa <- fa(fa_ddm, 5, rotate="oblimin", fm = "ml") 
print(sr.fa, cut=.3, sort = T) 

