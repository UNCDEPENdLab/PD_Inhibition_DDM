
# try running an Rscript that submits an array of qsubs using qsub -v --------

# initialize job log?
initialize <- FALSE

# for local debugging
ics <- 1

# run models or test job submission loop?
RUN = TRUE

#################
### setup params
pacman::p_load(tidyverse)

data_base <- ifelse(ics == 1,"/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed", "~/github_repos/PD_Inhibition_DDM/Data/preprocessed")
output_base <- "/gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition"
log_file <- ifelse(ics == 1, "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log.csv", "~/ics/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log.csv")
log_file_comp <- ifelse(ics == 1, "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv", "~/ics/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log_completed.csv")

nchains <- 2
nsamples <- paste0("samp",c(#1000, 
  #2000))#,#, 
#5000, 
10000))#, 20000, 40000, 80000))
tasks <- c(#"flanker")#,
  "recent_probes")#, "go_nogo")

full_sample <- c(#"clean_sample",
  "full_sample")
wt_scaling_factor <- .01
nburn_percentile <- .2
coding <- "acc"

#################
### done setting up params


## load in string of all possible models and select the ones that you wish to run
source(file.path(data_base,"../../Code/Functions/gen_supported_models.R"))
all_models <- list()

for(f in full_sample){all_models[[f]] <- gen_supported_models()}

models <- all_models

## load completed log and see who never finished running. More customized

# models <- gen_missing_mods(tasks,
#                            mod_log = log_file_comp, #this simply contains info on if jobs completed that is run after the fact.
#                            nsamples = nsamples,
#                            full_sample = full_sample)









# handle any necessary job-logging details --------------------------------


if (initialize) {
  job_info_log <- data.frame(job_id = 0, TASK = "initialize", SAMPLE = "initialize",MODEL = "initialize",CODE = "initialize",NCHAINS = 0, NBURN = 0, NSAMP = 0, DAY_SUB = Sys.Date(), TIME_SUB = Sys.time(), WT_REQUEST = 0, QSUB_STRING = "initialize")
  write.csv(job_info_log, file = log_file)
}

if(!RUN & initialize){
  job_info_log <- data.frame(job_id = 0, TASK = "initialize", SAMPLE = "initialize",MODEL = "initialize",CODE = "initialize",NCHAINS = 0, NBURN = 0, NSAMP = 0, DAY_SUB = Sys.Date(), TIME_SUB = Sys.time(),  WT_REQUEST = 0, QSUB_STRING = "initialize")
  job_id = 0} #if testing, just initialize a job counter and job log


# main worker function: ---------------------------------------------------


for(nsamp in nsamples){
  
  nsamp_numeric <- as.numeric(sub("samp", "", nsamp))
  nburn <- nsamp_numeric*nburn_percentile
  wt = wt_scaling_factor*nsamp_numeric
  
  #check if directory exists and if not make it exist.
  dir_string <- file.path(output_base,nsamp)
  if(!dir.exists(dir_string) & RUN){dir.create(dir_string)}
  
  for (task in tasks) {
    #check if directory exists and if not make it exist.
    dir_string <- file.path(output_base,nsamp, task)
    if(!dir.exists(dir_string) & RUN){dir.create(dir_string)}
    
    for(c in coding){
      for (subject_sample in full_sample) {
        #check if directory exists and if not make it exist.
        dir_string <- file.path(output_base,nsamp, task, subject_sample, "model_objects")
        if(!dir.exists(dir_string) & RUN){dir.create(dir_string)}
        
        #load relevant data. 
        rawdf <- file.path(data_base, paste0(task,"_",subject_sample, "_nafilt_",c,"Code.csv"))
        
        out <- dir_string
        for (m in models[[subject_sample]][[task]]) { # 5/6/20: expand to enable running different models depending on we want the full or the clean sample
          
          if(RUN){ #only mess with the log if you actually run the model, otherwise use a dummy counter.
            
            # the log keeps track of all parameters of the models being run for monitoring purposes
            job_info_log <- read.csv(log_file) %>% select(-X)
            
            job_id <- paste0("DDM_job_",job_info_log$job_id[nrow(job_info_log)] + 1)  
          } else{
            if(!initialize){
              job_info_log <- read.csv(log_file) %>% select(-X)
              
              job_id <- paste0("DDM_job_",job_info_log$job_id[nrow(job_info_log)] + 1)  
            }else{job_id <- job_id + 1}
          }
          
          
          cat("\nSubmitting ", job_id, "\n")
          
          # 4/9/20: per MH suggestion, run atomic python script that simply parallelizes over nchains, this ensures that the scheduler does not need to allocate a single job across nodes.
          cmd <- paste0("qsub -l nodes=1:ppn=",nchains," -l walltime=",wt,":00:00 -N ",job_id, " -A mnh5174_c_g_sc_default -v ",
                        "DF=",rawdf,
                        ",OUTDIR=",out,
                        ",TASK=",task,
                        ",MODELS='", m, "'",
                        ",CODE=", c,
                        ",NCHAINS=", as.character(nchains),
                        ",NBURN=", as.character(nburn),
                        ",NSAMP=", as.character(nsamp_numeric),
                        " qsub_inhibition_ddm_args.bash")
          
          
          job_info <- data.frame(job_id = job_info_log$job_id[nrow(job_info_log)] + 1, TASK = task,SAMPLE = subject_sample, MODEL = m, CODE = c, NCHAINS = nchains, NBURN = nburn, NSAMP = nsamp, DAY_SUB = as.character(Sys.Date()), TIME_SUB = as.character(Sys.time()), WT_REQUEST = wt, QSUB_STRING = cmd)
          
          job_info_log <- rbind(job_info_log, job_info)
          
          
          cat(cmd, "\n")
          if (RUN) {
            write.csv(job_info_log, file = log_file)
            system(cmd) # rock and roll
          }
        }
      }
    }
  }
}


