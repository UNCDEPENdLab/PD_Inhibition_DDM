diagnostic_plot <- function(post, label){
  postMCMC <- coda::as.mcmc(post)
  p_ESS <- as.numeric(coda::effectiveSize(postMCMC))
  
  summary_post <- describe_posterior(post,
                                     centrality = "map",
                                     # test = c("p_direction", "p_significance"),
                                     ci = 0.95,
                                     diagnostic = "all")
  
  
  
  
  auto <- data.frame(autocorrelation = as.numeric(autocorr(postMCMC, 0:100)), lag = 0:100)
  
  autoplot <- ggplot(auto, aes(x = lag, y = autocorrelation)) + geom_line() + theme_bw() + coord_cartesian(ylim = c(0,1))
  
  
  trace_df <- data.frame(p = post, nsample = 1:length(post)) 
  names(trace_df)[1] <- label
  
  posterior_distribution <- ggplot(trace_df, aes_string(x = label)) + geom_histogram(bins = 30) + geom_vline(xintercept = summary_post$MAP) + geom_vline(xintercept = summary_post$CI_low) + geom_vline(xintercept = summary_post$CI_high) + 
    labs(title = paste0("MAP estimate: ", round(summary_post$MAP,3)), subtitle = paste("ESS: ", round(p_ESS,1))) + theme_bw()
  posterior_trace <- ggplot(trace_df, aes_string(x = "nsample", y = label)) + geom_line() + theme_bw()
  
  left <- cowplot::plot_grid(posterior_trace, autoplot, ncol = 1)
  comb_plot <- cowplot::plot_grid(left, posterior_distribution, ncol = 2)
  # print(comb_plot)
  return(comb_plot)
}