
# quick and dirty pull posteriors and generate MAP estimates and stdErrs --------

pacman::p_load(janitor, tidyverse)

R.utils::sourceDirectory("~/ics/Nate/PD_Inhibition_DDM/Code/Functions")

# Flanker:
# DIC comparisons indicate the winning model in the flanker data as: v ~ stimulus + trial + prev_rt; st ~ 1

fl0 <- read.csv("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp10000/flanker/full_sample/diagnostics/v_trial_prev_rt_st_traces_0.csv")
fl1 <- read.csv("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp10000/flanker/full_sample/diagnostics/v_trial_prev_rt_st_traces_1.csv")
# fl2 <- read.csv("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp10000/flanker/full_sample/diagnostics/v_trial_prev_rt_st_traces.csv")

beepr::beep()
flank_posts <- list(fl0,fl1)

# Go/No-go:
# DIC comparisons indicate the winning model in the go/no-go data as: v ~ stimulus; st ~ 1; a ~ block_trial (number of trials since last no-go)

gng_posts <- list()

# for(i in 0:9){
for(i in 0:1){
  gng_posts[[i+1]] <- read.csv(paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp10000/go_nogo/full_sample/diagnostics/v_a_st_traces_",i,".csv"))
}

beepr::beep()

# Recent Probes:
# DIC comparisons indicate the winning model in the recent probes data as: v ~ stimulus; st ~ 1

rp_posts <- list()

# for(i in 0:9){
for(i in 0:1){
  rp_posts[[i+1]] <- read.csv(paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp10000/recent_probes/full_sample/diagnostics/v_st_traces_",i,".csv"))
}

beepr::beep()


# combine, diagnose, summarize, export -----------------------------------------------

uber <- list(flank_posts, gng_posts, rp_posts)



all_chains_descriptives <- list() #hold everything here
tictoc::tic()
for(i in 1:length(uber)){
  t_list <- uber[[i]]
  
  group_summaries <- list()
  subj_summaries <- list()
  group_combine_chains <- data.frame()
  subj_combine_chains <- data.frame()
  task_all <- list()
  for(c in 1:length(t_list)){
    chain <- t_list[[c]]
    chain <- chain %>% clean_names() %>% select(-x)
    # names(chain)
    group_nodes <- chain %>% select(!contains("subj")) 
    group_summaries[[c]] <- summarise_posteriors(posterior = group_nodes) %>% select(Parameter, MAP, CI_low, CI_high, pd, ESS, se)
    group_combine_chains <- rbind(group_combine_chains, group_nodes)
    
    subj_nodes <- chain %>% select(contains("subj"))
    subj_summaries[[c]] <- summarise_posteriors(posterior = subj_nodes) %>% select(Parameter, MAP, CI_low, CI_high, pd, ESS, se)
    subj_combine_chains <- rbind(subj_combine_chains, subj_nodes)
  }
  task_all[["group"]] <- list(group_combine_chains, group_summaries)
  task_all[["subj"]] <- list(subj_combine_chains, subj_summaries)
  all_chains_descriptives[[i]] <- task_all
}
tictoc::toc()
beepr::beep()


names(all_chains_descriptives) <- c("flanker", "go_nogo", "recent_probes")



subj_names <- sub("a_subj_", "", names(all_chains_descriptives[["flanker"]][["subj"]][[1]])[which(grepl("a_subj", names(all_chains_descriptives[["flanker"]][["subj"]][[1]])))])

for(s in subj_names){
  
  s <- "1" # for testing
  ess_tasks <- data.frame()
  for(t in names(all_chains_descriptives)){ # loop over tasks
    subj_chains <- list()
    
    for(i in 1:length(all_chains_descriptives[[t]][["subj"]][[2]])){ #loop over chains
      subj_posts_all <-  all_chains_descriptives[[t]][["subj"]][[2]][[i]] %>% filter(endsWith(Parameter, paste0("subj_", s)))  %>% filter(Parameter != paste0("t_subj_", s))
      
      if(t == "flanker"){subj_posts_all <- subj_posts_all %>% 
        mutate(Parameter = ifelse(Parameter == paste0("v_c_stim_treatment_congruent_t_incongruent_subj_",s), "v_incongruent",
                                  ifelse(Parameter == paste0("v_intercept_subj_", s), "v_congruent", Parameter)),
               Parameter = sub(paste0("_subj_", s), "", Parameter))}
      
      if(t == "go_nogo"){subj_posts_all <- subj_posts_all %>% 
        mutate(Parameter = ifelse(Parameter == paste0("v_c_stim_treatment_go_t_no_go_subj_",s), "v_nogo",
                                  ifelse(Parameter == paste0("v_intercept_subj_", s), "v_go",
                                         ifelse(Parameter == paste0("a_intercept_subj_", s), "a_intercept",
                                                ifelse(Parameter == paste0("a_block_trial_subj_", s), "a_trials_from_go",Parameter)))))}
      
      subj_chains[[i]] <- subj_posts_all 
      print(subj_posts_all)
    }
    
    ESS_total <- data.frame(Parameter = subj_chains[[1]]$Parameter, ESS = subj_chains %>%
                              map(~.[[6]]) %>% 
                              #or as suggested in the comments
                              #map(2) %>%
                              reduce(`+`)) %>% mutate(id = s,
                                                      task = t)
    ess_tasks <- rbind(ess_tasks, ESS_total)
  }
}


all_chains_descriptives[[t]][["subj"]][[2]][[i]] <- NULL
length(all_chains_descriptives[[t]][["subj"]][[2]])
nrow(all_chains_descriptives[[t]][["subj"]][[1]])
all_chains_descriptives[[t]][["subj"]][[1]] <- all_chains_descriptives[[t]][["subj"]][[1]][1:16000,]




for(t in names(all_chains_descriptives)){
  no_ter <- all_chains_descriptives[[t]][["subj"]][[1]] %>% select(!starts_with("t_subj"))  
  summ_all <- summarise_posteriors(no_ter)
  
  # unique(summ_all$Parameter)
  
  summ_all$id <- regmatches(summ_all$Parameter, regexpr("[[:digit:]]+", summ_all$Parameter))
  if(t == "flanker"){
    summ_flank <- summ_all %>%
      mutate(Parameter = ifelse(grepl("v_c_stim_treatment_congruent_t_incongruent_subj_",Parameter), "v_incongruent",
                                ifelse(grepl("v_intercept_subj_",Parameter), "v_intercept", 
                                       ifelse(grepl("a_subj_",Parameter), "a_intercept", 
                                              ifelse(grepl("v_prev_rt_subj_",Parameter), "v_prev_rt", 
                                                     ifelse(grepl("v_trial_z_subj_",Parameter), "v_trial", Parameter))))),
             task = t) %>% select(id, task, Parameter, MAP, se, CI_low, CI_high, pd) 
    # summ_flank %>% tibble() %>% print(n = 500)
    summ_flank[which(summ_flank$Parameter == "v_incongruent"),"MAP"] <- summ_flank[which(summ_flank$Parameter == "v_intercept"),"MAP"] + summ_flank[which(summ_flank$Parameter == "v_incongruent"),"MAP"]
    y <- abs(summ_gng[which(summ_flank$Parameter == "v_intercept"),"MAP"])
    
    x <- summ_gng[which(summ_flank$Parameter == "v_incongruent"),"MAP"]
    
    z <- data.frame(v_incongruent = x, v_congruent = y) %>% reshape2::melt() %>% rename(condition = variable, MAP = value)
    ggplot(data = z, aes(x = MAP, color = condition)) + geom_density()
    
  }
  
  if(t == "recent_probes"){
    summ_rp <- summ_all %>%
      mutate(Parameter = ifelse(grepl("v_c_stim_treatment_positive_t_negative_familiar_subj_",Parameter), "v_negative_familiar",
                                ifelse(grepl("v_intercept_subj_",Parameter), "v_intercept", 
                                       ifelse(grepl("a_subj_",Parameter), "a_intercept", 
                                              ifelse(grepl("v_c_stim_treatment_positive_t_negative_unfamiliar_subj_",Parameter), "v_negative_unfamiliar", 
                                                     ifelse(grepl("v_c_stim_treatment_positive_t_negative_rc_subj_",Parameter), "v_negative_response_conflict", 
                                                            ifelse(grepl("v_c_stim_treatment_positive_t_negative_highly_familiar_subj_",Parameter), "v_negative_highly_familiar", Parameter)))))),
             task = t) %>% select(id, task, Parameter, MAP, se, CI_low, CI_high, pd) 
    
    summ_rp[which(summ_flank$Parameter == "v_negative_familiar"),"MAP"] <- summ_flank[which(summ_flank$Parameter == "v_intercept"),"MAP"] + summ_flank[which(summ_flank$Parameter == "v_negative_familiar"),"MAP"]
    summ_rp[which(summ_flank$Parameter == "v_negative_unfamiliar"),"MAP"] <- summ_flank[which(summ_flank$Parameter == "v_intercept"),"MAP"] + summ_flank[which(summ_flank$Parameter == "v_negative_unfamiliar"),"MAP"]
    summ_rp[which(summ_flank$Parameter == "v_negative_highly_familiar"),"MAP"] <- summ_flank[which(summ_flank$Parameter == "v_intercept"),"MAP"] + summ_flank[which(summ_flank$Parameter == "v_negative_highly_familiar"),"MAP"]
    summ_rp[which(summ_flank$Parameter == "v_negative_response_conflict"),"MAP"] <- summ_flank[which(summ_flank$Parameter == "v_intercept"),"MAP"] + summ_flank[which(summ_flank$Parameter == "v_negative_response_conflict"),"MAP"]
    
    y <- summ_rp[which(summ_rp$Parameter == "v_intercept"),"MAP"]
    x <- summ_rp[which(summ_rp$Parameter == "v_negative_familiar"),"MAP"]
    a <- summ_rp[which(summ_rp$Parameter == "v_negative_unfamiliar"),"MAP"]
    b <- summ_rp[which(summ_rp$Parameter == "v_negative_highly_familiar"),"MAP"]
    c <- summ_rp[which(summ_rp$Parameter == "v_negative_response_conflict"),"MAP"]

    z <- data.frame(v_negative_familiar = x, v_positive = y, v_negative_unfamiliar = a, v_negative_highly_familiar = b, v_negative_response_conflict = c) %>% reshape2::melt() %>% rename(condition = variable, MAP = value)
    ggplot(data = z, aes(x = MAP, color = condition)) + geom_density()
    
  }
  
  if(t == "go_nogo"){
    summ_gng <- summ_all %>% 
      mutate(Parameter = ifelse(grepl("v_c_stim_treatment_go_t_no_go_subj_",Parameter), "v_nogo",
                                ifelse(grepl("v_intercept_subj_", Parameter), "v_intercept",
                                       ifelse(grepl("a_intercept_subj_", Parameter), "a_intercept",
                                              ifelse(grepl("a_block_trial_subj_", Parameter), "a_trials_from_go",Parameter)))),
             task = t) %>% select(id, task, Parameter, MAP, se, CI_low, CI_high, pd) 
    
    summ_gng[which(summ_gng$Parameter == "v_nogo"),"MAP"] <- summ_gng[which(summ_gng$Parameter == "v_intercept"),"MAP"] + summ_gng[which(summ_gng$Parameter == "v_nogo"),"MAP"]
    
    
    y <- abs(summ_gng[which(summ_gng$Parameter == "v_nogo"),"MAP"])
    # y <- summ_gng[which(summ_gng$Parameter == "v_nogo"),"MAP"]
    x <- summ_gng[which(summ_gng$Parameter == "v_intercept"),"MAP"]
    
    z <- data.frame(v_go = x, v_nogo = y) %>% reshape2::melt() %>% rename(condition = variable, MAP = value)
    ggplot(data = z, aes(x = MAP, color = condition)) + geom_density()
    
    # x <- summ_gng %>% group_by(id) %>% mutate(MAP = ifelse(Parameter == "v_nogo", select(filter(select(summ_gng, Parameter, MAP), Parameter == "v_go"), MAP) - select(filter(select(summ_gng, Parameter, MAP), Parameter == "v_nogo"), MAP), MAP))
    # x %>% print(n = 200)
    #  head(summ_gng)
    
  }
  
  
}


summs_all <- rbind(summ_flank, summ_gng, summ_rp)

summs_all <- summs_all %>% mutate(parameter = ifelse(startsWith(Parameter, "v_"), "v", "a"),
                     Parameter = ifelse(startsWith(Parameter, "v_"), sub("v_", "",Parameter), sub("a_", "",Parameter))) %>% rename(beta_contrast = Parameter)


summs_all <- summs_all %>% dplyr::select(-CI_low, -CI_high, -pd) %>% gather(key ="key", value = "value",MAP, se) %>% mutate(key = paste0(parameter,"_",key)) %>% dplyr::select(-parameter) %>% spread(key = "key", value =  "value")


write.csv(summs_all, file = "~/github_repos/PD_Inhibition_DDM/Code/brms/posterior_summaries.csv")
