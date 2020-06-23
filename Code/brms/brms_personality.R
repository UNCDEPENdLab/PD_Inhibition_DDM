
# Load Packages -----------------------------------------------------------


library(tidyverse)
library(brms)

# Load in data ------------------------------------------------------------


setwd("~/github_dirs/PD_Inhibition_DDM/code/brms/")
#load in and scale personality variables; makes it easier to provide weakly informative priors
pdata <- read.csv("~/github_dirs/PD_Inhibition_DDM/Code/factor_structure/AC_bsem_cfa_parcel_mod_savedata.csv") %>% dplyr::select(ANTAG.Median, ANTAG.Standard.Deviation, INHIB.Median, INHIB.Standard.Deviation, SUBJECT) %>% rename(id = SUBJECT)

#load in traces and wrangle so that only have subject specific parameters 
# may need to alter to be more extensible (in terms of wrangling so that getting everything you want
# potential room to have this be true across tasks (i.e., multiple fixed effects/random effects to allow for examining condition effects across tasks)
#n.b. this approach of combinging tasks may not be as preferable as a latent variable approach so that should be an explicit conversation
traces <- read.csv("posterior_summaries.csv") %>% dplyr::select(-X) %>% left_join(pdata) %>% mutate(a_available = if_else(is.na(a_MAP), FALSE, TRUE),
                                                                                                    v_available = if_else(is.na(v_MAP), FALSE, TRUE)) %>% mutate(a_MAP.c = if_else(is.na(a_MAP), 0, a_MAP),
                                                                                                                                                                 a_se.c = if_else(is.na(a_se), 0, a_se),
                                                                                                                                                                 v_MAP.c = if_else(is.na(v_MAP), 0, v_MAP),
                                                                                                                                                                 v_se.c= if_else(is.na(v_se), 0, v_se))


traces_base <- dplyr::filter(traces, beta_contrast == "intercept")
traces_rp <-dplyr::filter(traces, beta_contrast %in% c("response_conflict", "highly_familiar", "familiar"))
traces_uf <-dplyr::filter(traces, beta_contrast=="unfamiliar") #
# a condition effect doesn't mean anything > so only getting the intercept threshold effect
traces_ngn <-dplyr::filter(traces, beta_contrast == "nogo") #nothing
bf1 = bf(v_MAP | se(v_se) ~ ANTAG.Median  + INHIB.Median + (1|id/task))
bf2 = bf(a_MAP | se(a_se) ~ ANTAG.Median  + INHIB.Median + (1|id/task))

model <- brm(
  formula = bf1+bf2, 
  prior = c(set_prior("normal(0, 2)", class = "b")),
  iter = 1000, chains = 3,
  data = traces_base)
bf1 <- bf(v_MAP | se(v_se) ~  beta_contrast*me(ANTAG.Median, ANTAG.Standard.Deviation)  + beta_contrast*me(INHIB.Median, INHIB.Standard.Deviation) + (1|id/task))
bf2 <- bf(a_MAP | se(a_se) ~beta_contrast*me(ANTAG.Median, ANTAG.Standard.Deviation)   + beta_contrast*me(INHIB.Median, INHIB.Standard.Deviation) + (1|task/id))
bf1_out <- brm(bf1,
  iter = 50000, chains = 5,
  prior = c(set_prior("normal(0, 2)", class = "b")),
  data = dplyr::filter(traces, v_available = TRUE))
bf2_out <- brm(bf2,
               iter = iter_num, chains = chains_n,
               data = dplyr::filter(traces, a_available = TRUE))

# Loop over personality varaibles and append to list ----------------------

##set up iter number and number of chains
iter_num =1000
chains_n = 3

outputs <- list()




for(i in 1:length(p_vars)) {
  pvar_name <- p_vars[[i]]
  traces_tmp <- traces
  traces_tmp <- dplyr::select(traces_tmp, mean, se,key, subj)
  traces_tmp_pvar <- traces[[pvar_name]]
  traces_tmp$pvar <- traces_tmp_pvar
  v_out <- NULL
  a_out <- NULL
  t_out <- NULL
  v_out <- brm(
    mean | se(se) ~ pvar*key + (1|subj), 
    prior = c(set_prior("normal(0, 2)", coef = "pvar")),
    iter = iter_num, chains = chains_n,
    data = dplyr::filter(traces_tmp, str_detect(key, "^v"))
  )
  a_out <- brm(
    mean | se(se) ~ pvar + (1|subj), 
    prior = c(set_prior("normal(0, 2)", coef = "pvar")),
    iter = iter_num, chains = chains_n,
    data = dplyr::filter(traces_tmp, str_detect(key, "^a"))
  )
  
  t_out <- brm(
    mean | se(se) ~ pvar + (1|subj), 
    prior = c(set_prior("normal(0, 2)", coef = "pvar")),
    iter = iter_num, chains = chains_n,
    data = dplyr::filter(traces_tmp, str_detect(key, "^t"))
  )
  outputs[[pvar_name]] <- list("v_out" = v_out, "a_out" = a_out, "t_out" = t_out)
  
  
  if(!is.null(v_out)&& (slot(v_out$fit, "mode") != 2)){
    v_df <- data.frame(regressor = dimnames(fixef(v_out))[[1]],
                       mean = fixef(v_out)[,1],
                       ll = fixef(v_out)[,3],
                       ul = fixef(v_out)[,4],
                       rhat = summary(v_out)$fixed[,5],
                       ess = summary(v_out)$fixed[,6]) %>%
      mutate(sig = if_else( ll*ul > 0, 1,0), pvar = pvar_name, outcome = "v") 
    
  } else {
    v_df <- data.frame(regressor = NA_character_, mean = NA_real_, ll= NA_real_, ul = NA_real_, rhat = NA_real_, ess = NA_real_, sig = NA_real_, pvar = pvar_name, outcome = "v")
  }
  if(!is.null(a_out)&& (slot(a_out$fit, "mode") != 2)){
    a_df <- data.frame(regressor = dimnames(fixef(a_out))[[1]],
                       mean = fixef(a_out)[,1],
                       ll = fixef(a_out)[,3],
                       ul = fixef(a_out)[,4],
                       rhat = summary(a_out)$fixed[,5],
                       ess = summary(a_out)$fixed[,6]) %>%
      mutate(sig = if_else( ll*ul > 0, 1,0), pvar = pvar_name, outcome = "a") 
    
  } else {
    a_df <- data.frame(regressor = NA_character_, mean = NA_real_, ll= NA_real_, ul = NA_real_, rhat = NA_real_, ess = NA_real_, sig = NA_real_, pvar = pvar_name, outcome = "a")
  }
  
  if(!is.null(t_out) && (slot(t_out$fit, "mode") != 2)){
    t_df <- data.frame(regressor = dimnames(fixef(t_out))[[1]],
                       mean = fixef(t_out)[,1],
                       ll = fixef(t_out)[,3],
                       ul = fixef(t_out)[,4],
                       rhat = summary(t_out)$fixed[,5],
                       ess = summary(t_out)$fixed[,6]) %>%
      mutate(sig = if_else( ll*ul > 0, 1,0), pvar = pvar_name, outcome = "t") 
    
  } else {
    t_df <- data.frame(regressor = NA_character_, mean = NA_real_, ll= NA_real_, ul = NA_real_, rhat = NA_real_, ess = NA_real_, sig = NA_real_, pvar = pvar_name, outcome = "t")
  }
  tmp_df <- bind_rows(v_df, a_df, t_df)
  summary_df <- bind_rows(summary_df, tmp_df)
}
# this provides a summary of all the different models run 
summary_df <- dplyr::filter(summary_df, outcome != "999")


