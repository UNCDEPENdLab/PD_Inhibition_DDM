#####synthesize MCMC chains
setwd("~/Box/DEPENd/Projects/PD_Inhibition_DDM/")
library(dplyr)
library(psych)

##check DICs
DIC.dat <- read.csv("DICs_flanker_SNE.csv")
DIC.dat <- DIC.dat %>% select(-X)

data <- read.csv("flanker_model_outputs_COMBINED_INTERACTIONS_SNE.csv")
data <- data %>% dplyr::select(-X)
SNAP_dim <- unique(data$dimension)
dim_list_traces <- list()
dim_list_traces_betas <- list()

for (d in SNAP_dim){
  d_df <- data[which(data$dimension==d),]
  d_df <- d_df %>% dplyr::select(-dimension)
  dim_list_traces[[d]] <- d_df
  B0A <- mean(d_df$comb_stim_model_a_Int)
  B1A <- mean(d_df$comb_stim_model_a_SNAP)
  B2A <- mean(d_df$comb_stim_model_a_C)
  B3A <- mean(d_df$comb_stim_model_a_ixn)
  B0V <- mean(d_df$comb_stim_model_v_Int)
  B1V <- mean(d_df$comb_stim_model_v_SNAP)
  B2V <- mean(d_df$comb_stim_model_v_C)
  B3V <- mean(d_df$comb_stim_model_v_ixn)
  dim_list_traces_betas[[paste0(d, "_betas")]] <- data.frame(B0A, B1A, B2A, B3A, B0V, B1V, B2V, B3V)
}

#descriptives_df <- do.call(rbind, dim_list_traces_descriptives)
#write.csv(descriptives_df, file = "flanker_models_descriptives.csv")#, row.names = TRUE)



SNAP_info <- read.csv("Flanker_reduced_nooutliers_mnh.csv")
SNAP_info <- SNAP_info %>% group_by(subj_index) %>% filter(row_number() == 1) %>% dplyr::select(EXHIB, AGG, MISTRUST, MANIP, DISINH) %>% ungroup()
SNAP_info <- data.frame(SNAP_info)

preds <- c("EXHIB", "AGG", "MISTRUST", "MANIP", "DISINH")
df <- c()
for (v in preds) {
  qs <- quantile(as.matrix(SNAP_info[v]), c(.1, .5, .8))
  sumstat <- dim_list_traces_betas[[paste0(v, "_betas")]]
  vhatcon <- sumstat$B0V + sumstat$B1V*qs
  vhatincon <- sumstat$B0V + sumstat$B1V*qs + sumstat$B2V + sumstat$B3V*qs
  ahatcon <- sumstat$B0A + sumstat$B1A*qs
  ahatincon <- sumstat$B0A + sumstat$B1A*qs + sumstat$B2A + sumstat$B3A*qs
  stim <- c(rep("incon",3), rep("con",3))
  df_plot <- rbind(df, data.frame(pred=v, pred.quant = qs, stimulus = stim, vhat = c(vhatincon, vhatcon), ahat = c(ahatincon, ahatcon)))
}
df_plot

library(ggplot2)

pdf("SNE_interaction_plots.pdf", width = 10, height = 7)
for(p in preds){
  thisp <- df_plot[which(df_plot$pred == p) ,]
  v <- ggplot(thisp, aes(x=pred.quant, y = vhat, colour = stimulus)) +geom_line() + ylab("Drift Rate") +xlab(p)
  print(v)
  a <- ggplot(thisp, aes(x=pred.quant, y = ahat, colour = stimulus)) +geom_line() + ylab("Threshold") + xlab(p)
  print(a)
}
dev.off()


