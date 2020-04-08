summarise_posteriors <- function(posterior, split_subjects = FALSE){
  postMCMC <- coda::as.mcmc(posterior)
  p_ESS <- data.frame(coda::effectiveSize(postMCMC)) %>% rownames_to_column()
  names(p_ESS) <- c("Parameter", "ESS")
  
  summary_post <- describe_posterior(posterior,
                                     centrality = "all",
                                     dispersion = TRUE,
                                     ci = 0.95,
                                     diagnostic = "all") 
  
  summary_post <- summary_post %>% left_join(p_ESS, by = "Parameter")
  
  if(split_subjects){
    summary_post$subject <- do.call(rbind,strsplit(summary_post$Parameter, split = "_subj_"))[,2]
    summary_post$Parameter <- do.call(rbind,strsplit(summary_post$Parameter, split = "_subj_"))[,1]  
  }
  
  return(summary_post)
}
