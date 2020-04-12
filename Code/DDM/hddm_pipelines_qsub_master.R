
# try running an Rscript that submits an array of qsubs using qsub -v --------

#initialize job log?
initialize <- TRUE



#################
### setup params
pacman::p_load(tidyverse)

RUN = TRUE

data_base <- "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed"
output_base <- "/gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition"
log_file <- "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/pbs_outputs/PD_Inhibition_DDM_job_info_log.csv"

nchains <- 3
nsamples <- paste0("samp",c(1000))#, 2000, 5000, 10000, 20000, 40000, 80000))
tasks <- c("flanker", "recent_probes")#, "go_nogo")
full_sample <- c("clean_sample")#, "full_sample")
wt_scaling_factor <- .01
nburn_percentile <- .2
coding <- "acc"

models <- list(flanker = c("v"))#, 'vsv', 'v_block', 'v_blocksv', 'vst', 'vsvst', 'v_blockst', "v_blocksvst"), 
               # recent_probes = c("v", "vst", "vsv", "vsvst"))

if (initialize) {
  job_info_log <- data.frame(job_id = 0, TASK = "initialize", MODEL = "initialize",CODE = "initialize",NCHAINS = 0, NBURN = 0, NSAMP = 0, DAY_SUB = Sys.Date(), TIME_SUB = Sys.time(), qsub_string = "initialize")
  write.csv(job_info_log, file = log_file)
}
#################
### done setting up params



# main worker function: ---------------------------------------------------


for(nsamp in nsamples){

  nsamp_numeric <- as.numeric(sub("samp", "", nsamp))
  nburn <- nsamp_numeric*nburn_percentile
  wt = wt_scaling_factor*nsamp_numeric

  #check if directory exists and if not make it exist.
  dir_string <- file.path(output_base,nsamp)
  if(!dir.exists(dir_string)){dir.create(dir_string)}

  for (task in tasks) {
    #check if directory exists and if not make it exist.
    dir_string <- file.path(output_base,nsamp, task)
    if(!dir.exists(dir_string)){dir.create(dir_string)}

    for(c in coding){

      rawdf <- file.path(data_base, paste0(task,"_",c,"Code.csv"))

      for (subject_sample in full_sample) {
        #check if directory exists and if not make it exist.
        dir_string <- file.path(output_base,nsamp, task, subject_sample, "model_objects")
        if(!dir.exists(dir_string)){dir.create(dir_string)}

        out <- dir_string
        for (m in models[[task]]) {

          # the log keeps track of all parameters of the models being run for monitoring purposes
          job_info_log <- read.csv(log_file) %>% select(-X)

          job_id <- paste0("DDM_job_",job_info_log$job_id[nrow(job_info_log)] + 1)

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

          job_info <- data.frame(job_id = job_info_log$job_id[nrow(job_info_log)] + 1, TASK = task, MODEL = m, CODE = c, NCHAINS = nchains, NBURN = nburn, NSAMP = nsamp, DAY_SUB = as.character(Sys.Date()), TIME_SUB = as.character(Sys.time()), qsub_string = cmd)

          job_info_log <- rbind(job_info_log, job_info)
          write.csv(job_info_log, file = log_file)

          cat(cmd, "\n")
          if (RUN) {system(cmd)}
        }
      }
    }
  }
}


