
# digest log and check output directories ---------------------------------

library(pacman); p_load(tidyverse)

ics <- 0 
if(ics){
  log <- log_backup <- read.csv("/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log.csv") %>% select(-X) %>% filter(QSUB_STRING != "initialize")
  setwd("/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/")
} else{
  log <- log_backup <-  read.csv("/Users/natehall/ics/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log.csv") %>% select(-X) %>% filter(QSUB_STRING != "initialize")  
  setwd("/Users/natehall/ics/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/")
}
init_str <- "\\/gpfs\\/group\\/mnh5174\\/default\\/Nate\\/HDDM_outputs_PD_Inhibition\\/samp"

#demo on single string
# onstr <- as.character(log[2,"QSUB_STRING"])
# regmatches(onstr, gregexpr(pattern = paste0(init_str,".*model_objects"),text = onstr))[[1]]

# now try on whole column
log$OUTDIR <-  do.call(rbind,regmatches(as.character(log[,"QSUB_STRING"]), gregexpr(pattern = paste0(init_str,".*model_objects"),text = as.character(log[,"QSUB_STRING"]))))

if(!ics){
  log <- log %>% mutate(OUTDIR = sub("/gpfs/group/mnh5174/default","/Users/natehall/ics", OUTDIR, fixed = TRUE))
}


# loop over models tested -------------------------------------------------

verbose = TRUE
log$outputs_located <- NA

for (mod in 1:nrow(log)) {
  
  modR <- log[mod,] %>% select(-QSUB_STRING)
  
  path <- log$OUTDIR[mod]
  
  
  # list.dirs(path)
  files <- list.files(path)
  
  models <- paste0(modR$MODEL, "_chain", seq(0,modR$NCHAINS-1,1), "_", modR$CODE, "Code.model")
  dbs <- paste0(modR$MODEL, "_chain", seq(0,modR$NCHAINS-1,1), "_", modR$CODE, "Code.db")
  
  expected_files <- c(models, dbs)
  
  if(all(expected_files %in% files)){
    cat(paste0("mod ", mod, ": ", modR$MODEL, " --      All good!\n"))
    log$outputs_located[mod] <- "all"
  } else if(any(expected_files %in% files)){
    cat(paste0("mod ", mod, ": ", modR$MODEL, " --      Some missing\n"))
    log$outputs_located[mod] <- "some"
  } else{
    cat(paste0("mod ", mod, ": ", modR$MODEL, " --      All missing\n"))
    log$outputs_located[mod] <- "none"
  }
  
  
}

log <- log %>% select(-QSUB_STRING, -OUTDIR)

##get file info including most recent modification, this will equate roughly to run time for reference.


# this will only work on ics due to permissions not transferring over sshfs. If on local, logon to ACI and run manually and read below
if(ics) system("RScript gen_file_info.R")


finfo <- read.csv("pbs_file_info.csv") %>% dplyr::rename(job_id = X, end_time = ctime) %>% filter(grepl("DDM_job",job_id)) %>%
  mutate(job_id = sub("\\.out", "",sub("DDM_job_", "", job_id))) %>%
  select(job_id, end_time)
head(finfo,40)

log$job_id <- as.character(log$job_id)

log <- tibble(log)
log <- log %>% left_join(finfo, by = "job_id") %>% 
  mutate(end_time = as.POSIXct(end_time),
         TIME_SUB = as.POSIXct(TIME_SUB),
         run_time = as.numeric(end_time - TIME_SUB))


log %>% filter(outputs_located != "all") %>% print(n = 150)
log %>% filter(TASK == "recent_probes") %>% print(n = 150)

if(ics){
  write.csv(log, "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv")  
} else {
  write.csv(log, "/Users/natehall/ics/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv")  
}
log_backup