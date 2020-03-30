library(tidyverse)
library(coda)
library(brms)
setwd("~/github_dirs/PD_Inhibition_DDM/Outputs/practice_hddm/flanker/diagnostics/")
vreg_traces <- read.csv("v_reg_traces.csv") %>% dplyr::select(-X, -a, -t, -a_std, -t, -t_std, -v_Intercept, -v_Intercept_std) %>% 
  mutate(trace = 1:nrow(.))  %>% gather(key = "key", value = "value", -trace) %>% 
  mutate(key = gsub("v_C\\.stim\\.\\.Treatment\\.0\\.\\.\\.T\\.1\\.", "vC", key)) %>%
  mutate(key = gsub("v_Intercept", "vInt", key)) %>%
  dplyr::filter(!key %in% c("vC", "vC_std")) %>%
  separate(key, into = c("stat", "subj"), "_") %>% 
  mutate(subj = gsub("subj\\.", "", subj))

vreg_traces_summary <- vreg_traces %>% group_by(subj, stat) %>% 
  summarise(mean = mean(value),se = summary(as.mcmc(value))$statistics[["Time-series SE"]])  %>% 
  ungroup() %>% arrange(subj) %>% gather(key = "key", value = "value", mean, se) %>% 
  mutate(stat = paste0(stat, key)) %>% dplyr::select(-key) %>% spread(key = "stat", value = "value")


pdata <- read.csv("~/github_dirs/PD_Inhibition_DDM/Data/summarydat.csv") %>% dplyr::select(Subject, MPQAggression) %>% rename(subj = Subject)

vreg_traces_summary <- left_join(dplyr::mutate(vreg_traces_summary, subj = as.integer(subj)), pdata) 

vC_out <- brm(
  vCmean | se(vCse) ~ MPQAggression + (1 | subj), 
  prior = set_prior("uniform(0, 1000)", class = "sd"),
  iter = 4000, 
  data = vreg_traces_summary
)
