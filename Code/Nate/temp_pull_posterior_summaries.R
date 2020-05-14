

# interim messy solution for pulling model params -------------------------



rp_traces <- read.csv(paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp2000/recent_probes/full_sample/diagnostics/v_st_traces.csv"))
rp_traces <- rp_traces %>% select(contains("subj"))
rp_ddm_summary <- summarise_posteriors(rp_traces)#, split_subjects = TRUE)



rp_ddm_summary$Parameter <- ifelse(grepl("a_subj", rp_ddm_summary$Parameter), sub("a_subj.", "a_subj_", rp_ddm_summary$Parameter),
                                   ifelse( grepl("t_subj", rp_ddm_summary$Parameter), sub("t_subj.", "t_subj_",rp_ddm_summary$Parameter),
                                           ifelse( grepl("negative_familiar", rp_ddm_summary$Parameter), sub("v_C.stim..Treatment..positive....T.negative_familiar._subj.", "v.negative.familiar_subj_",rp_ddm_summary$Parameter),
                                                   ifelse( grepl("negative_rc", rp_ddm_summary$Parameter), sub("v_C.stim..Treatment..positive....T.negative_rc._subj.", "v.negative.rc_subj_",rp_ddm_summary$Parameter),
                                                           ifelse( grepl("negative_highly_familiar", rp_ddm_summary$Parameter), sub("v_C.stim..Treatment..positive....T.negative_highly_familiar._subj.", "v.negative.highly.familiar_subj_",rp_ddm_summary$Parameter),
                                                                   ifelse( grepl("negative_unfamiliar", rp_ddm_summary$Parameter), sub("v_C.stim..Treatment..positive....T.negative_unfamiliar._subj.", "v.negative.unfamiliar_subj_",rp_ddm_summary$Parameter),
                                                                           rp_ddm_summary$Parameter))))))



rp_ddm_summary$Parameter <- ifelse(grepl("Intercept", rp_ddm_summary$Parameter), gsub("v_Intercept_subj.", "v.positive_subj_",rp_ddm_summary$Parameter), rp_ddm_summary$Parameter)


rp_ddm_summary$subj_idx <-  do.call(rbind,strsplit(rp_ddm_summary$Parameter, split = "_"))[,3]
rp_ddm_summary$Parameter <- do.call(rbind,strsplit(rp_ddm_summary$Parameter, split = "_"))[,1]

rp_ddm_summary <- rp_ddm_summary %>% pivot_wider(id_cols = subj_idx, names_from = Parameter, values_from = c(MAP, SD))

names(rp_ddm_summary)[2:15] <- paste0("recent_probes_", names(rp_ddm_summary)[2:15])



# flanker -----------------------------------------------------------------

flank_dic <- read.csv(paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp2000/flanker/full_sample/diagnostics/dics_all.csv"))



dics  <- flank_dic %>% select(-X)

dics <- dics %>% mutate(DIC_diff = DIC -dics[1,2]) %>% arrange(DIC_diff) %>% distinct()


dic_diff <- ggplot(dics, aes(x = model, y = DIC_diff)) + geom_bar(stat = "identity", fill = "steelblue") + theme_bw() +
  geom_text(aes(label = round(DIC_diff,2)), vjust = -.5, color = "black")  +
  labs(y = expression(Delta*DIC), x = "Model")
print(dic_diff)

v_stimblock_trial_prev_rt_st
flank_traces <- read.csv(paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/samp2000/flanker/full_sample/diagnostics/v_stimblock_trial_prev_rt_st_traces.csv"))


flank_traces <- flank_traces %>% select(contains("subj"))
flank_ddm_summary <- summarise_posteriors(flank_traces)#, split_subjects = TRUE)



flank_ddm_summary$Parameter <- ifelse(grepl("a_subj", flank_ddm_summary$Parameter), sub("a_subj.", "a_subj_", flank_ddm_summary$Parameter), flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("v_Intercept", flank_ddm_summary$Parameter), sub("v_Intercept_subj.", "v.congruent_subj_",flank_ddm_summary$Parameter), flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("t_subj", flank_ddm_summary$Parameter), sub("t_subj.", "t_subj_",flank_ddm_summary$Parameter),flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("v_trial_z", flank_ddm_summary$Parameter), sub("v_trial_z_subj.", "v.trial_subj_",flank_ddm_summary$Parameter),flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("v_prev_rt", flank_ddm_summary$Parameter), sub("v_prev_rt_subj.", "v.prev.rt_subj_",flank_ddm_summary$Parameter), flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("most_incon", flank_ddm_summary$Parameter), sub("v_C.stimblock..Treatment..congruent....T.incongruent_most_incon._subj.", "v.Sincongruent.Bincon_subj_",flank_ddm_summary$Parameter), flank_ddm_summary$Parameter)
flank_ddm_summary$Parameter <- ifelse( grepl("most_con", flank_ddm_summary$Parameter), sub("v_C.stimblock..Treatment..congruent....T.incongruent_most_con._subj.", "v.Sincongruent.Bcon_subj_",flank_ddm_summary$Parameter), flank_ddm_summary$Parameter)

flank_ddm_summary$subj_idx <-  do.call(rbind,strsplit(flank_ddm_summary$Parameter, split = "_"))[,3]
flank_ddm_summary$Parameter <- do.call(rbind,strsplit(flank_ddm_summary$Parameter, split = "_"))[,1]

flank_ddm_summary <- flank_ddm_summary %>% pivot_wider(id_cols = subj_idx, names_from = Parameter, values_from = c(MAP, SD))

names(flank_ddm_summary)[2:15] <- paste0("flanker", names(flank_ddm_summary)[2:15])
flank_ddm_summary$subj_idx <- as.numeric(flank_ddm_summary$subj_idx)

# combine -----------------------------------------------------------------


ddm_params <- rp_ddm_summary %>% mutate(subj_idx = as.numeric(subj_idx)) %>% left_join(flank_ddm_summary, by = "subj_idx")
