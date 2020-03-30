library(tidyverse)
library(coda)
library(brms)

pdddm_home <- "~/Data_Analysis/PD_Inhibition_DDM"

vreg_traces <- read.csv(file.path(pdddm_home, "Outputs/practice_hddm/flanker/diagnostics", "v_reg_traces.csv.gz")) %>% dplyr::select(-X, -a, -t, -a_std, -t, -t_std, -v_Intercept, -v_Intercept_std) %>% 
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

pdata <- read.csv(file.path(pdddm_home, "Data", "summarydat.csv")) %>% dplyr::select(Subject, MPQAggression) %>% rename(subj = Subject)

vreg_traces_summary <- left_join(dplyr::mutate(vreg_traces_summary, subj = as.integer(subj)), pdata) %>%
  mutate(MPQAggression_z = as.vector(scale(MPQAggression)))

vC_out <- brm(
  vCmean | se(vCse) ~ 1 + MPQAggression_z + (1 | subj), 
  prior = set_prior("uniform(0, 100)", class = "sd"),
  iter = 32000, thin = 2,
  data = vreg_traces_summary
)

vInt_out <- brm(
  vIntmean | se(vIntse) ~ 1 + MPQAggression_z + (1 | subj), 
  prior = set_prior("uniform(0, 100)", class = "sd"),
  iter = 32000, thin = 2,
  data = vreg_traces_summary
)
