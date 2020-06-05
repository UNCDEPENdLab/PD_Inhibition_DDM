
# wrapper script to run model diagnostics on all HDDM analysis pipelines --------

########## still a working prototype!!
####################
## pipeline parameters to test
ics <- 0

# nsamples <- paste0("samp",c(1000, 2000, 5000, 10000, 20000, 40000, 80000))
nsamples <- "samp2000"
# tasks <- c("flanker", "recent_probes", "go_nogo")
tasks <- c("flanker", "recent_probes")
full_sample <- c("full_sample")#, "clean_sample")
data_type <- c("model_objects", "diagnostics")
create_plots <- TRUE

pacman::p_load(tidyverse)

# initialize large parameterization list ----------------------------------


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
      DIC_path <- file.path(diagnosdir, "dics_all.csv")
      
      for(mod in names(parameterizations[[task]])){
        parameterizations[[task]][["traces"]] <- paste0(diagnosdir, "/",mod, "_traces.csv")
        parameterizations[[task]][["gelman-rubin"]] <- paste0(diagnosdir, "/gr_",mod,".csv")
        parameterizations[[task]][["outdir"]] <- outdir
        parameterizations[[task]][["figuredir"]] <- figuredir
      }
      
      
      suff_stats_all_pipelines[[nsamp]][[task]][[subs]] <-  hddm_posterior_diagnostics(parameterizations[[task]], DICs = DIC_path, v_contrasts = TRUE,m_digest = "win")
    }
  }
}
