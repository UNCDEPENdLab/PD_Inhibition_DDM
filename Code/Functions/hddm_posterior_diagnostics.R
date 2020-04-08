## examine HDDM diagnostics and posteriors
## generate plots for group and subject posteriors and extract summary statistics for the requested posteriors

hddm_posterior_diagnostics <- function(diagnosdir, #where should the function look for model diagnostics? This should also include the model traces.
                                       models, # list with structure [[model]][[group/subject]] where the lowest level is a data.frame with 2 columns: parameter name (should be more descriptive) and label (output from HDDM) and a row for every (fixed and random, respectively) parameter in the model
                                       outdir = NULL, # where to write file with summary stats. If null, will simply return without writing.
                                       output_specifiers = c("",""), # string to append to outputs if requested. Takes two concatenated strings, the first of which is used to append to plots and second of which is appended to the summary statistics output.
                                       allowCache = FALSE, # allows for caching and reloading summary stats files. This will check outdir and return the relevant file. Importantly, if this is set to true, diagnostic figures will not be created, as these all require the full traces, which tend to be quite large and have decently long read times.
                                       m_digest = "win", # either "win" (only digest the winning model based on DIC comparisons), "all" (digest every model), or string or vector of strings with specific models to extract statistics from.
                                       fixed_random = "all", # group-level nodes or subject-specific nodes? (takes fixed, random, or all)
                                       create_plots = FALSE, # T/F print pdfs along the way?
                                       nCores = 1, # for processing models in parallel if desired
                                       gr = FALSE, # create histograms of gelman-rubin statistics for parameters in models? Doesn't calculate, just plots.
                                       v_contrasts = FALSE, # plot contrasts for drift rate. Finds v_intercept and appends other contrasts. All contrasts must start with "v_"
                                       ...
){
  
  # load required packages  ----------------------------
  
  
  pacman::p_load(tidyverse, bayestestR, coda, beepr, cowplot, psych, R.utils, qdap, foreach, doParallel, data.table, tictoc) 
  
  
  # compare DICS and select winning model ----------------------------------
  dics  <- read.csv(paste0(diagnosdir, "/dics_all.csv"))   
  
  dics <- dics %>% mutate(DIC_diff = DIC -dics[1,3])
  
  if(create_plots){
    pdf(file.path(figuredir,paste0("DIC_diff_plot",output_specifiers[1],".pdf")), width = 11, height = 8) 
    dic_diff <- ggplot(dics, aes(x = model, y = DIC_diff)) + geom_bar(stat = "identity", fill = "steelblue") + theme_bw() +
      geom_text(aes(label = round(DIC_diff,2)), vjust = -.5, color = "black")  +
      labs(y = expression(Delta*DIC), x = "Model")
    print(dic_diff)
    dev.off()
  }
  
  win_mod <- as.character(dics$model[which(dics$DIC == min(dics$DIC))])
  
  
  # select which models to loop over ----------------------------------------
  
  
  if(m_digest == "win"){
    m_digest <- win_mod
  } else if(p_digest == "all"){
    m_digest <- names(parameterizations)
  }
  
  # outer loop that iterates over all models requested -----------------------
  
  summary_stats_allmods <- list()   
  
  if(nCores != 1){
    message("Multi-thread processing not yet implemented")
  }
  # registerDoParallel(nCores)
  # 
  # summary_stats <- foreach(mod = m_digest, .packages = c("tidyverse", "bayestestR", "coda", "beepr", "cowplot", "psych", "R.utils", "qdap", "foreach", "doParallel", "data.table")) %dopar% {
  
  for(mod in m_digest){
    # Gelman-rubin Rhat statistics --------------------------------------------
    if(gr){
      
      grs <- read.csv(paste0(diagnosdir, "/gr_",win_mod,".csv")) %>% select(-X)
      
      if(create_plots){
        pdf(file.path(figuredir,paste0("gelman_rubin_hist",output_specifiers[1],".pdf")), width = 11, height = 8)
        gr_dist <- ggplot(grs, aes(x = rhat)) + geom_histogram() + theme_bw() + labs(y = "Frequency", x = "R")
        print(gr_dist)
        dev.off()
      }
      
    }
    # check for cached summary_stats files ------------------------------------
    
    sfile_path <- file.path(outdir, paste0("summary_stats_", mod, ifelse(output_specifiers[2] == "", "", paste0("_",output_specifiers[2])), ".RData"))
    if(allowCache && file.exists(sfile_path)){
      load(sfile_path)
    } else{
      # read and relabel group traces based on model input specification --------
      
      
      
      
      
      # tic(); 
      traces <- read.csv(file = paste0(diagnosdir, "/",mod, "_traces.csv")) %>% select(-X) 
      # toc(); beepr::beep()
      
      if(fixed_random %in% c("all", "fixed")){
        group_traces <- traces %>% select(as.character(models[[mod]][["group"]]$label))
        names(group_traces) <-  as.character(models[[mod]][["group"]]$param)
        
        summary_stats[["group"]] <- summarise_posteriors(group_traces)
        # group diagnostic plots --------------------------------------------------
        
        if(create_plots){
          pdf(file = paste0(figuredir,"/group_diagnostics_",mod, "_", output_specifiers[1],".pdf"), width =11, height =8)
          for(i in names(group_traces)){
            p <- diagnostic_plot(group_traces[,i], i)
            print(p)
          }
          dev.off()
          
          # compare drift rate between conditions: group ----------------------------
          
          if(v_contrasts){
            pdf(file = paste0(figuredir,"/v_contrasts_",mod, "_", output_specifiers[1],".pdf"), width =11, height =8)
            v_traces <- group_traces %>% select(starts_with("v_")) %>% select(-ends_with("_std"))
            intercept <- as.character(models[[mod]][["group"]]$param[which(models[[mod]][["group"]]$label == "v_Intercept")])
            contrasts <- names(select(v_traces, -intercept))
            
            contrast_plot <- plot_v_contrasts(v_traces, intercept, contrasts)
            print(contrast_plot)
            dev.off()
          }
          
          
        }
        
      }
      
      # grab subject-level traces and process  ---------------------------------------------------------
      
      if(fixed_random %in% c("all", "random")){
        subj_traces <- list()
        sub_summary_stats <- list()
        for(s_param in as.character(models[[mod]][["subjects"]]$label)){
          
          # relabel subject traces based on model input specification ---------------
          
          
          s_param_name <- as.character(models[[mod]][["subjects"]][which(models[[mod]][["subjects"]][,2] == s_param),1])
          print(paste0("processing posteriors for: ", s_param_name))
          
          subj_traces[[s_param_name]] <- traces %>% select(starts_with(s_param))
          
          # convert colnames to desired parameter names 
          names(subj_traces[[s_param_name]]) <- mgsub(as.character(models[[mod]][["subjects"]]$label), paste0(as.character(models[[mod]][["subjects"]]$param), "_"),names(subj_traces[[s_param_name]]))
          subj_traces[[s_param_name]] <- subj_traces[[s_param_name]] %>% mutate(nsample = 1:nrow(.))
          
          # extract summary of posterior distribution: subject-level  ------------------------------
          
          
          sub_summary_stats[[s_param]] <- summarise_posteriors(select(subj_traces[[s_param_name]], -nsample), split_subjects = TRUE) 
          
          if(create_plots){
            
            # subject-level diagnostic plots ------------------------------------------
            
            pdf(file = file.path(figuredir,paste0("subject_diagnostics_",mod, "_", sub("_subj", "", s_param_name), ifelse(output_specifiers[1] == "","","_"),output_specifiers[1], ".pdf")), width =11, height =8)
            for(p in names(subj_traces[[s_param_name]])[which(names(subj_traces[[s_param_name]]) != "nsample")]){
              dplot <- diagnostic_plot(subj_traces[[s_param_name]][,p], p)
              print(dplot)
            }
            dev.off()
          }
        }  
        
        # combine subject-level summaries
        summary_stats[["subjects"]] <- do.call(rbind,sub_summary_stats) %>% arrange(as.numeric(subject))
        
        
        # compare drift rate between conditions: subjects ----------------------------
        
        if(v_contrasts){
          pdf(file = paste0(figuredir,"/subjects_v_contrasts_",mod, "_", output_specifiers[1],".pdf"), width =11, height =8)
          
          v_names <- names(subj_traces)[grepl("v_", names(subj_traces))]
          
          
          intercept <- as.character(models[[mod]][["subjects"]]$param[which(models[[mod]][["subjects"]]$label == "v_Intercept_subj.")])
          contrasts <- v_names[which(v_names != intercept)]
          
          for(sub in 1:ncol(select(subj_traces[[v_names[1]]], -nsample))){
            # sub <- 1
            v_traces <- data.frame()
            
            v_traces <- do.call(cbind, lapply(subj_traces[v_names], function(cond){
              column_name <- colnames(cond)[sub]
              select(cond, column_name)
            }))
            
            snum <- sub(intercept, "", colnames(v_traces)[1])
            
            
            contrast_plot <- plot_v_contrasts(v_traces, paste0(intercept, snum), paste0(contrasts,snum))
            print(contrast_plot)
            
          }
          dev.off()
        }
      }
      
    }
    
    if(!is.null(outdir)){save(summary_stats,file = sfile_path)}
    summary_stats_allmods[[mod]] <- summary_stats
  }
  
  
  
  
  return(summary_stats_allmods)
}
