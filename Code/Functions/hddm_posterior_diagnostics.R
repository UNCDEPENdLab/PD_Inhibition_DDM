## examine HDDM diagnostics and posteriors
## generate plots for group and subject posteriors and extract summary statistics for the requested posteriors

hddm_posterior_diagnostics <- function(models, # list with structure [[model]][[group/subject/traces/gelman-rubin/outdir/figuredir]] where the lowest level is a data.frame with 2 columns: parameter name (should be more descriptive) and label (output from HDDM) and a row for every (fixed and random, respectively) parameter in the model. Must include an element [[model]][["path"]] to the traces object(s) output from HDDM and calculated gelman-rubin statistics if requested (optional). Can also include where to write the output (summary statistics) and the figures generated buy the function. If both are not included will set these both to the working directory.
                                       output_specifiers = c("",""), # string to append to outputs if requested. Takes two concatenated strings, the first of which is used to append to plots and second of which is appended to the summary statistics output.
                                       allowCache = FALSE, # allows for caching and reloading summary stats files. This will check outdir and return the relevant file. Importantly, if this is set to true, diagnostic figures will not be created, as these all require the full traces, which tend to be quite large and have decently long read times.
                                       DICs = NULL, # path to .csv containing information about deviance information criterion (DIC) for model comparison if requested, otherwise bypass.
                                       m_digest = "all", # either "win" (only digest the winning model based on DIC comparisons), "all" (digest every model), or string or vector of strings with specific models to extract statistics from.
                                       fixed_random = "all", # group-level nodes or subject-specific nodes? (takes fixed, random, or all)
                                       create_plots = TRUE, # T/F print pdfs along the way?
                                       nCores = 1, # for processing models in parallel if desired
                                       v_contrasts = FALSE, # plot contrasts for drift rate. Finds v_intercept and appends other contrasts. All contrasts must start with "v_". Still need to figure out what to do with interactions, as these are less straightforward.
                                       ...
){
  
  # load required packages  ----------------------------
  
  
  pacman::p_load(tidyverse, bayestestR, coda, beepr, cowplot, psych, R.utils, qdap)#, foreach, doParallel, data.table, tictoc) 
  
  # compare DICS and select winning model ----------------------------------
  if(!is.null(DICs)){
    dics  <- read.csv(DICs)   
    
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
  }
  
  
  # select which models to loop over ----------------------------------------
  
  
  if(m_digest == "win"){
    m_digest <- win_mod
  } else if(m_digest == "all"){
    m_digest <- names(models)
  }
  
  
  # outer loop that iterates over all models requested
  
  summary_stats_allmods <- list()   
  
  if(nCores != 1){
    message("Multi-thread processing not yet implemented")
  } else{
    # registerDoParallel(nCores)
    # 
    # summary_stats <- foreach(mod = m_digest, .packages = c("tidyverse", "bayestestR", "coda", "beepr", "cowplot", "psych", "R.utils", "qdap", "foreach", "doParallel", "data.table")) %dopar% {
    
  }
  
  for(mod in m_digest){
    
    
    # check for elements with paths that lead to desired output/figure directories--------
    summary_stats_outstring <- paste0("summary_stats_", mod, ifelse(output_specifiers[2] == "", "", paste0("_",output_specifiers[2])), ".RData")
    
    if(!"outdir" %in% names(models[[mod]])){
      models[[mod]][["outdir"]] <- getwd()
    }
    outdir <- models[[mod]][["outdir"]]
    message("summary statistics will write to: ", outdir, "/", summary_stats_outstring)
    
    if(create_plots){
      if(!"figuredir" %in% names(models[[mod]])){
        models[[mod]][["figuredir"]] <- getwd()
      }
      figuredir <- models[[mod]][["figuredir"]]
      message("figure .pdfs will be written to: ", figuredir)
    }
    
    # Gelman-rubin Rhat statistics --------------------------------------------
    if("gelman-rubin" %in% names(models[[mod]])){
      
      grs <- read.csv(models[[mod]][["gelman-rubin"]]) %>% select(-X)
      
      if(create_plots){
        message("Creating histogram of Gelman-Rubin statistics")
        pdf(file.path(figuredir,paste0("gelman_rubin_hist",output_specifiers[1],".pdf")), width = 11, height = 8)
        gr_dist <- ggplot(grs, aes(x = rhat)) + geom_histogram() + theme_bw() + labs(y = "Frequency", x = "R")
        print(gr_dist)
        dev.off()
      }
      
    }
    
    sfile_path <- file.path(outdir, summary_stats_outstring)
    if(allowCache && file.exists(sfile_path)){ #check for cached summary_stats files
      load(sfile_path)
    } else{
      # grab group-level traces and process --------
      
      tictoc::tic()
      traces <- read.csv(file = models[[mod]][["traces"]]) %>% select(-X) 
      tictoc::toc(); #beepr::beep()
      
      
      
      summary_stats <- list()
      # browser()
      if(fixed_random %in% c("all", "fixed")){
        group_traces <- traces %>%   select(as.character(models[[mod]][["group"]]$label))
        names(group_traces) <-  as.character(models[[mod]][["group"]]$param)
        
        message("Digesting group posteriors")
        
        summary_stats[["group"]] <- summarise_posteriors(group_traces)
        # group diagnostic plots --------------------------------------------------
        
        if(create_plots){
          group_pdf_str <- ifelse(output_specifiers[1] == "",
                                  paste0(figuredir,"/group_diagnostics_",mod,".pdf"),
                                  paste0(figuredir,"/group_diagnostics_",mod, "_", output_specifiers[1],".pdf"))
          
          pdf(file =group_pdf_str, width =11, height =8)
          for(i in names(group_traces)){
            p <- diagnostic_plot(group_traces[,i], i)
            print(p)
          }
          dev.off()
          
          # compare drift rate between conditions: group ----------------------------
          
          if(v_contrasts){
            
            vcont_pdf_str <- ifelse(output_specifiers[1] == "",
                                    paste0(figuredir,"/v_contrasts_",mod,".pdf"),
                                    paste0(figuredir,"/v_contrasts_",mod, "_", output_specifiers[1],".pdf"))
            
            pdf(file = vcont_pdf_str, width =11, height =8)
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
          message(paste0("processing subject posteriors for: ", s_param_name))
          
          subj_traces[[s_param_name]] <- traces %>% select(starts_with(s_param))
          
          # convert colnames to desired parameter names 
          names(subj_traces[[s_param_name]]) <- mgsub(as.character(models[[mod]][["subjects"]]$label), paste0(as.character(models[[mod]][["subjects"]]$param), "_"),names(subj_traces[[s_param_name]]))
          subj_traces[[s_param_name]] <- subj_traces[[s_param_name]] %>% mutate(nsample = 1:nrow(.))
          
          # extract summary of posterior distribution: subject-level  ------------------------------
          
          
          sub_summary_stats[[s_param]] <- summarise_posteriors(select(subj_traces[[s_param_name]], -nsample), split_subjects = TRUE) 
          
          if(create_plots){
            
            # subject-level diagnostic plots ------------------------------------------
            
            sub_pdf_str <- ifelse(output_specifiers[1] == "",
                                    paste0(figuredir,"/subject_diagnostics_",s_param_name,"_",mod,".pdf"),
                                    paste0(figuredir,"/subject_diagnostics_",s_param_name,"_", mod, "_", output_specifiers[1],".pdf"))
            
            pdf(file = file.path(sub_pdf_str), width =11, height =8)
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
          sub_vcont_str <- ifelse(output_specifiers[1] == "",
                                paste0(figuredir,"/subjects_v_constrasts_",mod,".pdf"),
                                paste0(figuredir,"/subject_v_contrasts_",mod, "_", output_specifiers[1],".pdf"))
          
          pdf(file = sub_vcont_str, width =11, height =8)
          
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
