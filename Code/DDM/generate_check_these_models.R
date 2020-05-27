
# examine which 10,000 sample chains we need to look through --------------

basedir <- "~/ics/Nate/HDDM_outputs_PD_Inhibition"; setwd(basedir)
# tasks <- list.dirs(basedir)
tasks <- c("flanker", "go_nogo", "recent_probes")
samples <- c("samp10000", "samp2000")

df.all <- data.frame(task = character(0), samples = character(0), model = character(0))


for(t in tasks){
  for(s in samples){
    mods <- list.files(paste0(basedir, "/", s, "/",   t, "/full_sample/model_objects"))
    
    #check to make sure all chains have finished running for all models
    mds_all <- list()
    for(i in 0:9){
      nam <- paste0("chain", i)
      mds_all[[nam]] <- sub(paste0("_chain",i,"_accCode.model"), "", mods[grepl(paste0("_chain",i,"_accCode.model"), mods)])
      print(mds_all[[nam]]) 
      cat("\n")
      
    }
    mds_check <- unique(do.call(c, mds_all))
    
    df.check <- data.frame(task = rep(t, length(mds_check)), samples = rep(s, length(mds_check)), model = mds_check )
    df.all <- rbind(df.all, df.check)
    
  }
    
}




#trim to specific models from the past that we want to look at


df.all$post_processed <- "no"

write.csv(df.all, file = "~/ics/Nate/PD_Inhibition_DDM/Data/cache/check_these_models.csv", row.names = FALSE)


