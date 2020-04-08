#plots posterior distributions of drift rates that can vary by some contrast of interest (group, condition)

plot_v_contrasts <- function(v_traces,intercept,contrasts){
  for(i in contrasts){v_traces[,i] <- v_traces[,i] + v_traces[,intercept]} 
  v_melt <- reshape2::melt(v_traces) %>% rename(condition = variable, v_posterior = value)
  plot_obj <- ggplot(v_melt, aes(x = v_posterior, color = condition)) + geom_density()
  return(plot_obj)
}
