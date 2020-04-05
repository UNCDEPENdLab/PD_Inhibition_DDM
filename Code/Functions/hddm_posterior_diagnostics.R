## examine HDDM diagnostics and posteriors

task <- "flanker"
# task <- "recent_probes"
full_models <- FALSE; full <- ifelse(full_models, "full", "practice")
create_plots <- TRUE
basedir <- "~/ics/Nate/PD_Inhibition_DDM/"; setwd(basedir)
diagnosdir_full <- paste0(basedir, "Outputs/full_hddm/", task, "/diagnostics")
diagnosdir_prac <- paste0(basedir, "Outputs/practice_hddm/", task, "/diagnostics")
# diagnosdir_full <- paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/full_hddm/", task, "/diagnostics")
# diagnosdir_prac <- paste0("~/ics/Nate/HDDM_outputs_PD_Inhibition/practice_hddm/", task, "/diagnostics")

list.files(diagnosdir_prac)

if(task == "recent_probes"){
  models <- c("v_reg", "vst_reg", "vsv_reg", "vsvst_reg")    
} else if(task == "flanker"){
  models <- c("v_reg","vst_reg", "vsv_reg", "vsvst_reg", "v_block_reg", "v_blockst_reg", "v_blocksv_reg", "v_blockst_reg", "v_blocksvst_reg")
}

pacman::p_load(tidyverse, tictoc, bayestestR)

# compare DICS and select winning model ----------------------------------
if(full_models){
  dics  <- read.csv(paste0(diagnosdir_full, "/dics_all.csv"))   
} else{
  dics <- read.csv(paste0(diagnosdir_prac, "/dics_all.csv"))  
 }

dics <- dics %>% mutate(DIC_diff = DIC -dics[1,3])

if(create_plots){
  pdf(paste0(basedir,"Figures/DIC_diff_plot_", task, "_", full, ".pdf"), width = 11, height = 8) #this is so dumb.
  # pdf(paste0("/Users/natehall/ics/Nate/PD_Inhibition_DDM/Figures/DIC_diff_plot_", task, "_", full, ".pdf"), width = 11, height = 8)
  
  dic_diff <- ggplot(dics, aes(x = model, y = DIC_diff)) + geom_bar(stat = "identity", fill = "steelblue") + theme_bw() + geom_text(aes(label = round(DIC_diff,2)), vjust = -.5, color = "black")  
  print(dic_diff)
  dev.off()
}

win_mod <- as.character(dics$model[which(dics$DIC == min(dics$DIC))])
# Gelman-rubin Rhat statistics --------------------------------------------

if(full_models){
  grs <- read.csv(paste0(diagnosdir_full, "/gr_",win_mod,".csv")) %>% select(-X)
} else {
  grs <- read.csv(paste0(diagnosdir_prac, "/gr_",win_mod,".csv")) %>% select(-X)
}

if(create_plots){
  pdf(paste0(basedir,"Figures/gelman_rubin_hist_", task, "_", full, ".pdf"), width = 11, height = 8)
  gr_dist <- ggplot(grs, aes(x = rhat)) + geom_histogram()
  print(gr_dist)
  dev.off()
}


# grs %>% arrange(-rhat) %>% head(10)


# plot out traces ---------------------------------------------------------
if(full_models){
  tic(); traces <- read.csv(file = paste0(diagnosdir_full, "/",win_mod, "_traces.csv")); toc(); beepr::beep()  
} else{
  tic(); traces <- read.csv(file = paste0(diagnosdir_prac, "/",win_mod, "_traces.csv")); toc(); beepr::beep()  
}

# str(traces)
traces <- traces %>% select(-X)

subj_traces <- list()

if(task == "recent_probes"){
  # split up traces for ease of use
  subj_traces[["a"]] <- traces %>% select(starts_with("a_subj"))
  subj_traces[["t"]] <- traces %>% select(starts_with("t_subj"))
  subj_traces[["v_positive"]] <- traces %>% select(starts_with("v_Intercept_subj"))
  subj_traces[["v_neg_fam"]] <- traces %>% select(starts_with("v_C.Condition..Treatment..positive....T.n_fam._subj")); subj_traces[["v_neg_fam"]] <- subj_traces[["v_neg_fam"]] + subj_traces[["v_positive"]]
  subj_traces[["v_neg_highfam"]] <- traces %>% select(starts_with("v_C.Condition..Treatment..positive....T.n_hfam._subj")); subj_traces[["v_neg_highfam"]] <- subj_traces[["v_neg_highfam"]] + subj_traces[["v_positive"]]
  subj_traces[["v_neg_rc"]] <- traces %>% select(starts_with("v_C.Condition..Treatment..positive....T.n_rc._subj")); subj_traces[["v_neg_rc"]] <- subj_traces[["v_neg_rc"]] + subj_traces[["v_positive"]]
  subj_traces[["v_neg_unfam"]] <- traces %>% select(starts_with("v_C.Condition..Treatment..positive....T.n_unfam._subj")); subj_traces[["v_neg_unfam"]] <- subj_traces[["v_neg_unfam"]] + subj_traces[["v_positive"]]
  
  
  group_traces <- traces %>% select(-starts_with("a_subj"), -starts_with("t_subj"), -starts_with("v_Intercept_subj"),
                                    -starts_with("v_C.Condition..Treatment..positive....T.n_fam._subj"),
                                    -starts_with("v_C.Condition..Treatment..positive....T.n_hfam._subj"),
                                    -starts_with("v_C.Condition..Treatment..positive....T.n_rc._subj"),
                                    -starts_with("v_C.Condition..Treatment..positive....T.n_unfam._subj")) %>% mutate(nsample = 1:nrow(.))
  
  
}


pdf(file = paste0(basedir,"Figures/posteriors_",task,"/group_traces_",task,"_",win_mod, "_",full, ".pdf"), width =11, height =8)
for(i in names(group_traces)[which(names(group_traces) != "nsample")]){
  map <- as.numeric(map_estimate(group_traces[,i]))
  posterior_distribution <- ggplot(group_traces, aes_string(x = i)) + geom_histogram() + geom_vline(xintercept = map) + labs(title = paste0("MAP estimate:", round(map,3)))
  posterior_trace <- ggplot(group_traces, aes_string(x = "nsample", y = i)) + geom_line()
  comb_plot <- cowplot::plot_grid(posterior_trace, posterior_distribution)
  print(comb_plot)
}
dev.off()

# pdf(file = paste0(basedir,"Figures/posteriors_",task,"/subject_traces_",task,".pdf"), width =11, height =8)

subj_suff_stats <- list()
for(p in names(subj_traces)){
  if(create_plots){pdf(file = paste0(basedir,"Figures/posteriors_",task,"/subject_traces_",task,"_", p,".pdf"), width =11, height =8)}
  
  print(p)
  subj_traces[[p]] <- subj_traces[[p]] %>% mutate(nsample = 1:nrow(.))
  
  param_df <- data.frame()
  for(sub in colnames(subj_traces[[p]])[which(colnames(subj_traces[[p]]) != "nsample")]){
    # describe_posterior(subj_traces[[p]][,sub], centrality = "map")
    map <- as.numeric(map_estimate(subj_traces[[p]][,sub]))
    if(create_plots){
      posterior_distribution <- ggplot(subj_traces[[p]], aes_string(x = sub)) + geom_histogram() + geom_vline(xintercept = map) + labs(title = paste0("MAP estimate:", round(map,3)))
      posterior_trace <- ggplot(subj_traces[[p]], aes_string(x = "nsample", y = sub)) + geom_line()
      comb_plot <- cowplot::plot_grid(posterior_trace, posterior_distribution)
      print(comb_plot)
    }
        
    sd <- sd(subj_traces[[p]][,sub])
    number <- as.numeric(gsub("\\D+", "",sub))
    sub.df <- data.frame(id = number, map = map, sd = sd)
    
    colnames(sub.df)[c(2,3)] <- c(paste0(p,"_map"), paste0(p,"_sd"))
    param_df <- rbind(param_df, sub.df)
  }
  subj_suff_stats[[p]] <- param_df
  if(create_plots){dev.off()}
}


save(subj_suff_stats,file = paste0(basedir, "/Data/cache/", task, "_subj_sufficient_stats.RData"))  
subj_suff_stats  
  
  
 
  
  
  
  
  
  
