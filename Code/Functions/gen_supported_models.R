
# create all possible HDDM parameterizations for PD DDM project -----------
gen_supported_models <- function(...){
supported_models <- list(flanker = c('v',
                                     'v_stimblock',
                                     'v_block',
                                     'v_trial',
                                     'v_runtrial',
                                     'v_prev_rt',
                                     'v_stimblock_trial',
                                     'v_stimblock_runtrial',
                                     'v_stimblock_prev_rt',
                                     'v_block_trial',
                                     'v_block_runtrial',
                                     'v_block_prev_rt',
                                     'v_trial_runtrial',
                                     'v_trial_prev_rt',
                                     'v_runtrial_prev_rt',
                                     'v_stimblock_trial_runtrial',
                                     'v_stimblock_trial_prev_rt',
                                     'v_stimblock_runtrial_prev_rt',
                                     'v_block_trial_runtrial',
                                     'v_block_trial_prev_rt',
                                     'v_block_runtrial_prev_rt',
                                     'v_stimblock_trial_runtrial',
                                     'v_stimblock_trial_prev_rt',
                                     'v_stimblock_runtrial_prev_rt',
                                     'v_block_trial_runtrial_prev_rt',
                                     'v_stimblock_trial_runtrial_prev_rt',
                                     'v_st',
                                     'v_stimblock_st',
                                     'v_block_st',
                                     'v_trial_st',
                                     'v_runtrial_st',
                                     'v_prev_rt_st',
                                     'v_stimblock_trial_st',
                                     'v_stimblock_runtrial_st',
                                     'v_stimblock_prev_rt_st',
                                     'v_block_trial_st',
                                     'v_block_runtrial_st',
                                     'v_block_prev_rt_st',
                                     'v_trial_runtrial_st',
                                     'v_trial_prev_rt_st',
                                     'v_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_st',
                                     'v_stimblock_trial_prev_rt_st',
                                     'v_stimblock_runtrial_prev_rt_st',
                                     'v_block_trial_runtrial_st',
                                     'v_block_trial_prev_rt_st',
                                     'v_block_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_st',
                                     'v_stimblock_trial_prev_rt_st',
                                     'v_stimblock_runtrial_prev_rt_st',
                                     'v_block_trial_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_prev_rt_st',
                                     'v_st',
                                     'v_stimblock_st',
                                     'v_block_st',
                                     'v_trial_st',
                                     'v_runtrial_st',
                                     'v_prev_rt_st',
                                     'v_stimblock_trial_st',
                                     'v_stimblock_runtrial_st',
                                     'v_stimblock_prev_rt_st',
                                     'v_block_trial_st',
                                     'v_block_runtrial_st',
                                     'v_block_prev_rt_st',
                                     'v_trial_runtrial_st',
                                     'v_trial_prev_rt_st',
                                     'v_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_st',
                                     'v_stimblock_trial_prev_rt_st',
                                     'v_stimblock_runtrial_prev_rt_st',
                                     'v_block_trial_runtrial_st',
                                     'v_block_trial_prev_rt_st',
                                     'v_block_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_st',
                                     'v_stimblock_trial_prev_rt_st',
                                     'v_stimblock_runtrial_prev_rt_st',
                                     'v_block_trial_runtrial_prev_rt_st',
                                     'v_stimblock_trial_runtrial_prev_rt_st',
                                     'v_sv_st',
                                     'v_stimblock_sv_st',
                                     'v_block_sv_st',
                                     'v_trial_sv_st',
                                     'v_runtrial_sv_st',
                                     'v_prev_rt_st',
                                     'v_stimblock_trial_sv_st',
                                     'v_stimblock_runtrial_sv_st',
                                     'v_stimblock_prev_rt_sv_st',
                                     'v_block_trial_sv_st',
                                     'v_block_runtrial_sv_st',
                                     'v_block_prev_rt_sv_st',
                                     'v_trial_runtrial_sv_st',
                                     'v_trial_prev_rt_sv_st',
                                     'v_runtrial_prev_rt_sv_st',
                                     'v_stimblock_trial_runtrial_sv_st',
                                     'v_stimblock_trial_prev_rt_sv_st',
                                     'v_stimblock_runtrial_prev_rt_sv_st',
                                     'v_block_trial_runtrial_sv_st',
                                     'v_block_trial_prev_rt_sv_st',
                                     'v_block_runtrial_prev_rt_sv_st',
                                     'v_stimblock_trial_runtrial_sv_st',
                                     'v_stimblock_trial_prev_rt_sv_st',
                                     'v_stimblock_runtrial_prev_rt_sv_st',
                                     'v_block_trial_runtrial_prev_rt_sv_st',
                                     'v_stimblock_trial_runtrial_prev_rt_sv_st'),
                         recent_probes = c("v", 
                                           "v_st", 
                                           "v_sv",
                                           "v_sv_st",
                                           "v_cond", 
                                           "v_cond_st", 
                                           "v_cond_sv",
                                           "v_cond_sv_st"
                                           ),
                         go_nogo = c("v",
                                     "v_st",
                                     "v_sv",
                                     "v_sv_st",
                                     "v_cond",
                                     "v_cond_st",
                                     "v_cond_sv",
                                     "v_cond_sv_st",
                                     "v_stim_cond",
                                     "v_stim_cond_st",
                                     "v_stim_cond_sv",
                                     "v_stim_cond_sv_st",
                                     "a",
                                     "a_st",
                                     "a_sv",
                                     "a_sv_st",
                                     "v_a",
                                     "v_a_st",
                                     "v_a_sv",
                                     "v_a_sv_st",
                                     "v_a_cond",
                                     "v_a_cond_st",
                                     "v_a_cond_sv",
                                     "v_a_cond_sv_st",
                                     "v_a_stim_cond",
                                     "v_a_stim_cond_st",
                                     "v_a_stim_cond_sv",
                                     "v_a_stim_cond_sv_st"))
return(supported_models)
}


# if some models time out pull these  -------------------------------------

gen_missing_mods <- function(task,
                             mod_log = "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv",
                             full_sample = "full_sample",
                             nsamples = "samp2000"){
  
  message("Only configured to handle one value of nsamples!!")
  # browser()
  df <- read.csv(mod_log) %>% select(-X) %>% filter(outputs_located != "all" & TASK %in% task & NSAMP %in% nsamples & SAMPLE %in% full_sample) 
  
  out <- list()
  for(s in full_sample){
    for (t in task) {
      out[[s]][[t]] <- as.character(df[ which(df$SAMPLE %in% s & df$TASK %in% t),"MODEL"])
    }
  }
  
  
  return(out)
  
  
  
}
