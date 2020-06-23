
# Load Packages -----------------------------------------------------------


library(tidyverse)
library(coda)
library(brms)

# Load in data ------------------------------------------------------------


setwd("~/github_dirs/PD_Inhibition_DDM/code/brms/")
#load in and scale personality variables; makes it easier to provide weakly informative priors
pdata <- read.csv("~/github_dirs/PD_Inhibition_DDM/Data/summarydat.csv") %>% dplyr::select(-X, -contains("Acc"), -contains("RT")) %>% rename(subj = Subject) %>% mutate_at(vars(contains("MPQ", "T")), list(scale))

#load in traces and wrangle so that only have subject specific parameters 
# may need to alter to be more extensible (in terms of wrangling so that getting everything you want
# potential room to have this be true across tasks (i.e., multiple fixed effects/random effects to allow for examining condition effects across tasks)
#n.b. this approach of combinging tasks may not be as preferable as a latent variable approach so that should be an explicit conversation
traces <- read.csv("v_reg_traces.csv.gz") %>% dplyr::select(contains("subj")) %>% 
  mutate(trace = 1:nrow(.))  %>% gather(key = "key", value = "value", -trace) %>% 
  mutate(key = gsub("_C\\.stim\\.\\.Treatment\\.0\\.\\.\\.T\\.1\\.", "C", key)) %>% #this particularly may need to become more extensible
  mutate(key = gsub("_Intercept", "Int", key)) %>%
  separate(key, into = c("stat", "subj"), "_") %>% 
  mutate(subj = gsub("subj\\.", "", subj)) %>%
  group_by(subj, stat) %>% 
  summarise(mean = mean(value),se = summary(as.mcmc(value))$statistics[["Time-series SE"]])  %>% 
  ungroup() %>% arrange(subj) %>% gather(key = "key", value = "value", mean, se) %>% 
  mutate(stat = paste0(stat, "_", key)) %>% dplyr::select(-key) %>% spread(key = "stat", value = "value") %>%
  gather(key = "key", value = "value", -subj) %>% separate(key, into = c("key", "stat"), sep = "_") %>% 
  spread(key = "stat", value = "value") %>% mutate(subj = as.integer(subj)) %>% left_join(pdata)
p_vars <- names(pdata)[grepl("MPQ|T_", names(pdata))]



# Loop over personality varaibles and append to list ----------------------

##set up iter number and number of chains
iter_num =1000
chains_n = 3

outputs <- list()
summary_df <- data.frame(regressor = "999", mean = 999, ll = 999, ul = 999, rhat = 999, ess = 999, sig = 999, pvar = "999", outcome = "999")

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


