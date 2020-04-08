
# try running an Rscript that submits an array of qsubs using qsub -v --------


nodes <- 1
ppn <- 1
wt <- 4
rawdf <- "/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv"
out <- "/gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition/samp1000/flanker/clean_sample"
task <- "flanker"
job_id <- "test"
nchains <- 10
nburn <- 100
nsamp <- 1000
models <- paste("v" , 'vsv', 'v_block', 'v_blocksv', 'vst', 'vsvst', 'v_blockst', "v_blocksvst")


# system(paste0("qsub -v NODES=\"",nodes,
#        "\",PPN"=,ppn,
#        ",WT=", wt,
#        ",JOB_ID_NH=",job_id,
#        ",DF=",rawdf,
#        ",OUTDIR=",out,
#        ",TASK=",task,
#        ",MODELS=", models,
#        ",NCHAINS=", as.character(nchains),
#        ",NBURN=", as.character(nburn),
#        ",NSAMP=", as.character(nsamp),
#        " qsub_inhibition_ddm_args.bash"))



command <- shQuote(paste0("qsub -v NODES=\"",nodes,"\",PPN=\"",ppn, "\" test_script.bash"))
system(command)
# paste0('qsub -v NODES=\"',nodes),'\",PPN'=,ppn)
# 
# paste0('qsub -v NODES=\"',nodes,'\",PPN=\"',ppn,"\" test_script.bash")
