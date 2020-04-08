
# wrapper script to run model diagnostics on all HDDM analysis pipelines --------

########## still a working prototype!!
####################
## pipeline parameters to test
ics <- 0

nsamples <- paste0("samp",c(1000, 2000, 5000, 10000, 20000, 40000, 80000))
tasks <- c("flanker", "recent_probes", "go_nogo")
full_sample <- c("full_sample", "clean_sample")
data_type <- c("model_objects", "diagnostics")
create_plots <- TRUE

pacman::p_load(tidyverse)

# initialize large parameterization list ----------------------------------

##generated a 3-deep list that takes the structure [[task]][[parameterization]][[group/subjs]] that is utilized internally in the function when extracting and relabelling subject and group level posteriors.

## modify if alternate parameterizations are fit
parameterizations <- list()
for(task in tasks){
  if(task == "recent_probes"){
    
    p_labels <- c("v_reg", "vst_reg", "vsv_reg", "vsvst_reg")   
    message(cat("Recent probes models to examine: ", p_labels))
    
    for(param in p_labels){
      ## group-params. Be careful when specifying these!!
      gps_orig <- gps <- data.frame(param = c("a", "a_std", "t", "t_std", 
                                              "v_positive", "v_positive_std", 
                                              "v_neg_fam", "v_neg_fam_std", 
                                              "v_neg_highfam", "v_neg_highfam_std", 
                                              "v_neg_rc", "v_neg_rc_std", 
                                              "v_neg_unfam", "v_neg_unfam_std"),
                                    label = c("a", "a_std", "t", "t_std", 
                                              "v_Intercept", "v_Intercept_std", 
                                              "v_C.Condition..Treatment..positive....T.n_fam.", "v_C.Condition..Treatment..positive....T.n_fam._std", 
                                              "v_C.Condition..Treatment..positive....T.n_hfam.", "v_C.Condition..Treatment..positive....T.n_hfam._std", 
                                              "v_C.Condition..Treatment..positive....T.n_rc.", "v_C.Condition..Treatment..positive....T.n_rc._std",
                                              "v_C.Condition..Treatment..positive....T.n_unfam.", "v_C.Condition..Treatment..positive....T.n_unfam._std"))
      
      #adapt for parameterizations that go beyond simple stimulus type.
      if(grepl("sv",param)){
        gps <- gps %>% rbind(data.frame(param = "sv", label = "sv_Intercept"))  
      }
      if(grepl("st",param)){
        gps <- gps %>% rbind(data.frame(param = "st", label = "st_Intercept"))  
      }
      
      parameterizations[[task]][[param]][["group"]] <- gps
      
      ##subject-params
      parameterizations[[task]][[param]][["subjects"]] <- data.frame(param = paste0(gps_orig$param, "_subj"), label = paste0(gps_orig$label, "_subj.")) %>% filter(!grepl("_std", label))
    }
  } else if(task == "flanker"){
    p_labels <- c("v_reg","vst_reg", "vsv_reg", "vsvst_reg", "v_block_reg", "v_blockst_reg", "v_blocksv_reg", "v_blockst_reg", "v_blocksvst_reg")
    message(cat("Flanker models to examine: ", p_labels))
    
    for(param in p_labels){
      ## group-params. Be careful when specifying these!!
      if(!grepl("block", param)){
        gps_orig <- gps <- data.frame(param = c("a", "a_std", "t", "t_std", 
                                                "v_congruent", "v_congruent_std", 
                                                "v_incongruent", "v_incongruent_std"),
                                      label = c("a", "a_std", "t", "t_std", 
                                                "v_Intercept", "v_Intercept_std", 
                                                "v_C.stim..Treatment.0...T.1.", "v_C.stim..Treatment.0...T.1._std"))
      } else{
        gps_orig <- gps <- data.frame(param = c("a", "a_std", "t", "t_std", 
                                                "v_congruent_blockC", "v_congruent_blockC_std", 
                                                "v_incongruent_blockC", "v_incongruent_blockC_std",
                                                "v_congruent_blockI", "v_congruent_blockI_std", 
                                                "v_incongruent_blockI", "v_incongruent_blockI_std"
        ),
        label = c("a", "a_std", "t", "t_std", 
                  "v_Intercept", "v_Intercept_std", 
                  "v_C.stim..Treatment.0...T.1.", "v_C.stim..Treatment.0...T.1._std",
                  "v_C.CongruentBlock..Treatment.0...T.1.", "v_C.CongruentBlock..Treatment.0...T.1._std",
                  "v_C.stim..Treatment.0...T.1..C.CongruentBlock..Treatment.0...T.1.", "v_C.stim..Treatment.0...T.1..C.CongruentBlock..Treatment.0...T.1._std"
        ))
      }
      
      #adapt for parameterizations that go beyond simple stimulus type.
      if(grepl("sv",param)){
        gps <- gps %>% rbind(data.frame(param = "sv", label = "sv_Intercept"))  
      }
      if(grepl("st",param)){
        gps <- gps %>% rbind(data.frame(param = "st", label = "st_Intercept"))  
      }
      
      parameterizations[[task]][[param]][["group"]] <- gps
      
      ##subject-params
      parameterizations[[task]][[param]][["subjects"]] <- data.frame(param = paste0(gps_orig$param, "_subj"), label = paste0(gps_orig$label, "_subj.")) %>% filter(!grepl("_std", label))
    }
    
    
  } else if(task == "go_nogo"){
    message("Go/No-go models not fit yet.")
  }
}

########## 
####################
##

##setup directory strings and source HDDM functions

basedir <- ifelse(ics == 0,"~/github_repos/PD_Inhibition_DDM", "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM"); setwd(basedir)
hddm_outputdir <- ifelse(ics == 0,"~/ics/Nate/HDDM_outputs_PD_Inhibition", "/gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition")

R.utils::sourceDirectory(file.path(basedir, "Code/Functions/"))


# uber loop that performs posterior checks --------------------------------

suff_stats_all_pipelines <- list()

for(nsamp in nsamples){
  for(task in tasks){
    for(subs in full_sample){
      
      diagnosdir <- file.path(hddm_outputdir,nsamp, task, subs, "diagnostics")
      figuredir <- file.path(basedir,"Figures",nsamp, task, subs)
      outdir <- file.path(basedir, "Outputs", "posterior_summaries", nsamp, task, subs)
      
      suff_stats_all_pipelines[[nsamp]][[task]][[subs]] <-  hddm_posterior_diagnostics(diagnosdir, parameterizations[[task]], outdir, allowCache = TRUE)
    }
  }
}
