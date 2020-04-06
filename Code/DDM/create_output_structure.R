
# create hierarchical output directory structure  -------------------------

nsamples <- paste0("samp",c(1000, 2000, 5000, 10000, 20000, 40000, 80000))
task <- c("flanker", "recent_probes", "go_nogo")
full_sample <- c("full_sample", "clean_sample")
data_type <- c("model_objects", "diagnostics")

for(i in nsamples){
  dir_string <- file.path("~/ics/Nate/HDDM_outputs_PD_Inhibition",i)
  if(!dir.exists(dir_string)){dir.create(dir_string)}
  for(j in task){
    dir_string <- file.path("~/ics/Nate/HDDM_outputs_PD_Inhibition",i,j)
    if(!dir.exists(dir_string)){dir.create(dir_string)}
    for(k in full_sample){
      dir_string <- file.path("~/ics/Nate/HDDM_outputs_PD_Inhibition",i,j,k)
      if(!dir.exists(dir_string)){dir.create(dir_string)}
      for (l in data_type) {
        dir_string <- file.path("~/ics/Nate/HDDM_outputs_PD_Inhibition",i,j,k,l)
        if(!dir.exists(dir_string)){dir.create(dir_string)}
      }
    }
  }
}


# figures and outputs -----------------------------------------------------------------

nsamples <- paste0("samp",c(1000, 2000, 5000, 10000, 20000, 40000, 80000))
task <- c("flanker", "recent_probes", "go_nogo")
full_sample <- c("full_sample", "clean_sample")
# data_type <- c("diagnostics")

for(i in nsamples){
  # dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Figures/",i)
  dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Outputs/posterior_summaries",i)
  if(!dir.exists(dir_string)){dir.create(dir_string)}
  for(j in task){
    # dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Figures/",i,j)
    dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Outputs/posterior_summaries",i,j)
    if(!dir.exists(dir_string)){dir.create(dir_string)}
    for(k in full_sample){
      # dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Figures/",i,j,k)
      dir_string <- file.path("~/github_repos/PD_Inhibition_DDM/Outputs/posterior_summaries",i,j,k)
      if(!dir.exists(dir_string)){dir.create(dir_string)}
      # for (l in data_type) {
      #   dir_string <- file.path("~/ics/Nate/HDDM_outputs_PD_Inhibition",i,j,k,l)
      #   if(!dir.exists(dir_string)){dir.create(dir_string)}
      # }
    }
  }
}
