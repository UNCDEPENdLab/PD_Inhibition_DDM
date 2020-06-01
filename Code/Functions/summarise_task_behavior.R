# load and calc summaries of behavioral variables -------------------------
summarise_task_behavior <- function(basedir = "~/github_dirs/PD_Inhibition_DDM/"){
  # flanker -----------------------------------------------------------------
  
  
  fraw <- read.csv(paste0(basedir, "/Data/preprocessed/flanker_full_sample_accCode.csv")) %>% tibble()
  flank_rt_stim <- fraw %>% group_by(subj_idx, stim) %>% summarise(flank_rt = mean(rt_inv_trim_grp, na.rm = TRUE)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = flank_rt) %>% 
    rename(flanker_rt_cong = congruent, flanker_rt_incong = incongruent)
  flank_acc_stim <-fraw %>% group_by(subj_idx, stim) %>% summarise(flank_acc = sum(response, na.rm = TRUE)/length(response)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = flank_acc) %>%
    rename(flanker_acc_cong = congruent, flanker_acc_incong = incongruent)
  flank_rt_stimblock <- fraw %>% group_by(subj_idx, stim,block) %>% summarise(flank_rt = mean(rt_inv_trim_grp, na.rm = TRUE)) %>%
    pivot_wider(id_cols = subj_idx, names_from = c(stim,block), values_from = flank_rt) %>% 
    rename(flanker_rt_Scong_Bcong = congruent_most_con,
           flanker_rt_Scong_Bincong = congruent_most_incon,
           flanker_rt_Sincong_Bcong = incongruent_most_con,
           flanker_rt_Sincong_Bincong = incongruent_most_incon)
  flank_acc_stimblock <- fraw %>% group_by(subj_idx, stim,block) %>% summarise(flank_acc = sum(response, na.rm = TRUE)/length(response)) %>%
    pivot_wider(id_cols = subj_idx, names_from = c(stim,block), values_from = flank_acc) %>% 
    rename(flanker_acc_Scong_Bcong = congruent_most_con,
           flanker_acc_Scong_Bincong = congruent_most_incon,
           flanker_acc_Sincong_Bcong = incongruent_most_con,
           flanker_acc_Sincong_Bincong = incongruent_most_incon)
  
  flanker_summaries <- flank_rt_stim %>% left_join(flank_rt_stimblock, by = "subj_idx") %>% left_join(flank_acc_stim, by = "subj_idx") %>% left_join(flank_acc_stimblock, by = "subj_idx")
  f_cors <- cor(flanker_summaries, use = "pairwise.complete.obs")
  ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
  ggsave(paste0(basedir, "/Figures/flanker_behav_corrs.pdf"), width = 11, height = 8)
  
  
  
  # recent probes -----------------------------------------------------------
  
  
  rraw <- read.csv(paste0(basedir, "/Data/preprocessed/recent_probes_full_sample_accCode.csv")) %>% tibble()
  rp_rt_stim <- rraw %>% group_by(subj_idx, stim) %>% summarise(rp_rt = mean(rt_log_trim_grp, na.rm = TRUE)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = rp_rt) 
  names(rp_rt_stim)[2:6] <- paste0("recent_probes_rt_", names(rp_rt_stim)[2:6])
  
  rp_acc_stim <-rraw %>% group_by(subj_idx, stim) %>% summarise(rp_acc = sum(response, na.rm = TRUE)/length(response)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = rp_acc) 
  names(rp_acc_stim)[2:6] <- paste0("recent_probes_acc_", names(rp_acc_stim)[2:6])
  
  recent_probes_summaries <- rp_rt_stim %>% left_join(rp_acc_stim, by = "subj_idx")
  f_cors <- cor(recent_probes_summaries, use = "pairwise.complete.obs")
  ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
  ggsave(paste0(basedir, "/Figures/recent_probes_behav_corrs.pdf"), width = 11, height = 8)
  
  
  # recent probes -----------------------------------------------------------
  
  
  graw <- read.csv(paste0(basedir, "/Data/preprocessed/go_nogo_full_sample_accCode.csv")) %>% tibble()
  #drop no-gos for RT calculations
  graw <- graw %>% filter(rt_log_trim_grp != 0)
  gng_rt_stim <- graw %>% group_by(subj_idx, stim) %>% summarise(gng_rt = mean(rt_log_trim_grp, na.rm = TRUE)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = gng_rt) 
  names(gng_rt_stim)[2:3] <- paste0("go_nogo_rt_", names(gng_rt_stim)[2:3])
  
  gng_acc_stim <-graw %>% group_by(subj_idx, stim) %>% summarise(gng_acc = sum(response, na.rm = TRUE)/length(response)) %>% 
    pivot_wider(id_cols = subj_idx, names_from = stim, values_from = gng_acc) 
  names(gng_acc_stim)[2:3] <- paste0("go_nogo_acc_", names(gng_acc_stim)[2:3])
  
  gng_rt_cond <- graw %>% group_by(subj_idx, stim,cond) %>% summarise(gng_rt = mean(rt_log_trim_grp, na.rm = TRUE)) %>%
    pivot_wider(id_cols = subj_idx, names_from = c(stim,cond), values_from = gng_rt) 
  names(gng_rt_cond)[2:9] <- paste0("go_nogo_rt_", names(gng_rt_cond)[2:9])
  
  gng_acc_cond <- graw %>% group_by(subj_idx, stim,cond) %>% summarise(gng_acc = sum(response, na.rm = TRUE)/length(response)) %>%
    pivot_wider(id_cols = subj_idx, names_from = c(stim,cond), values_from = gng_acc) 
  names(gng_acc_cond)[2:9] <- paste0("go_nogo_acc_", names(gng_acc_cond)[2:9])
  
  go_nogo_summaries <- gng_rt_stim %>% left_join(gng_rt_cond, by = "subj_idx") %>% left_join(gng_acc_stim, by = "subj_idx") %>% left_join(gng_acc_cond, by = "subj_idx")
  f_cors <- cor(go_nogo_summaries, use = "pairwise.complete.obs")
  ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
  ggsave(paste0(basedir, "/Figures/go_nogo_behav_corrs.pdf"), width = 11, height = 8)
  
  
  
  
  # output df with everything  ----------------------------------------------
  
  behav_summs_all <- flanker_summaries %>% full_join(recent_probes_summaries, by = "subj_idx") %>% full_join(go_nogo_summaries, by = "subj_idx")
  f_cors <- cor(behav_summs_all, use = "pairwise.complete.obs")
  ggcorrplot(f_cors, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
  ggsave(paste0(basedir, "/Figures/behav_corrs.pdf"), width = 11, height = 8)
  
  #v pretty
  return(behav_summs_all)
}












