#####synthesize MCMC chains
#setwd("~/Box Sync/DEPENd/Projects/PD_Inhibition_DDM/")
setwd("~/Box Sync/DEPENd/Projects/PD_Inhibition_DDM/Outputs/")
library(dplyr)
library(psych)
library(tidyr)
library(tidyverse)

##check DICs
# DIC.dat <- read.csv("DICs_baseline_flanker_092617.csv")
# DIC.dat[order(DIC.dat$DIC),]
# #DIC.dat <- DIC.dat %>% select(-X)

#setwd("/mnt/ics/PD_Inhibition_DDM_bashfiles_outputs/full_models/")
setwd("~/Desktop/")
data <- read.csv("flanker_posteriors_10052017.csv")
data <- data %>% dplyr::select(-X)
SNAP_dim <- as.character(unique(data$dimension))
SNAP_dim <- c("AGG", "DISINH", "ECCPERC", "ENTITL", "EXHIB", "IMPUL", "LOSLFEST", "POSTEMP", "SELFHARM")
dim_list_traces <- list()
dim_list_descriptives <- list()
dim_list_traces_betas <- list()
dim_list_traces_betas_SEs <- list()

#StE <- function(x){sd(x)/sqrt(length(x))}

#pull in snap data
SNAP_info <- read.csv("~/Box Sync/DEPENd/Projects/PD_Inhibition_DDM/Flanker_reduced_stim_coded_nooutliers.csv")

# ggplot_v_list <- list()
# ggplot_a_list <- list()

#MH: prototype one var
df_all_v <- data.frame()
for (dim in SNAP_dim){
d_df <- dplyr::filter(data, dimension== dim) %>% dplyr::select(-dimension)

#plots of individual distributions for manual checks
# hist(d_df$vincon)
# hist(d_df$vcon)
# hist(d_df$v_SNAP)
# hist(d_df$v_ixn)
# hist(d_df$v_SNAP + d_df$v_ixn)
# hist(d_df$v_SNAP)

snapvar <- SNAP_info %>% group_by(subj_index) %>% 
  filter(row_number() == 1) %>% dplyr::select(subj_index, dim) %>% arrange(subj_index) %>% ungroup()
qs <- data.frame(qs=c(.2, .5, .8), vals=quantile(snapvar[[dim]], c(.2, .5, .8)))

#cartesian project: replicate posteriors for each of the quantile targets
d_df <- merge(d_df, qs)
d_df_v <- d_df %>% mutate(
  Congruent=v_Intercept + v_SNAP*vals,
  Incongruent=v_Intercept + v_SNAP*vals + v_C + v_ixn*vals
) %>% gather(key="condition", value="vhat", Congruent, Incongruent)

d_df_v$dim <- rep(dim, length(d_df_v[,1]))
#d_df_a <- d_df %>% mutate(ahat = a_Intercept + a_SNAP*vals)
df_all_v <- rbind(df_all_v, d_df_v)

}

#data %>% group_by(dimension) %>% summarize(avSNAPcorr=cor(a_SNAP, v_SNAP), intcor=cor(a_Intercept, v_Intercept))

####try some bullshit 

compute_stats <- function(vec) {
  require(coda)
  post <- coda::mcmc(vec)
  hp <- coda::HPDinterval(post, prob = .9)
  c(m=mean(vec), lo=hp[1], hi=hp[2])
}
raw.SNAP <- read.csv("/Users/natehall/Box Sync/DEPENd/Projects/PD_Inhibition_DDM/Data_from_MNH/alldat.csv")
raw.SNAP <- raw.SNAP %>% select(SNAP_dim)

# dim <- "DISINH"
##insert for loop here
toplot <- data.frame()
for (dim in SNAP_dim){
SNAP_vec <- seq(from = min(SNAP_info[,dim]), to = max(SNAP_info[,dim]), by = .05)
SNAP_vec <- data.frame(SNAP_vec)
big_df <- data[which(data$dimension == dim),] 
big_df <- merge(big_df, SNAP_vec)
big_df_t <- big_df %>% transmute(
  ahat = a_Intercept + a_SNAP*SNAP_vec,
  vhat_inc = v_Intercept + v_C + v_SNAP + v_ixn*SNAP_vec,
  vhat_con = v_Intercept + v_SNAP*SNAP_vec,
  dim = SNAP_vec
) %>% ungroup()

bigsplit <- split(big_df_t, big_df$SNAP_vec)

# hist(big_df$a_SNAP)
slope_rework_allstats <- do.call(rbind, lapply(bigsplit, function(subdf) {
  #browser()
  vars <- grep("ahat", names(subdf), value=TRUE)
  df <- sapply(vars, function(v) {
    compute_stats(subdf[[v]])
    #data.frame(age=subdf$age[1], stats)
  })
  df <- as.data.frame(df)
  df$statistic <- rownames(df)
  #make the data frame sane
  df$SNAP <- subdf$dim[1]
  #df <- t(df)
  
  df_reshape <- df %>% gather(key="slope", value="value", -statistic, -SNAP) %>%
    spread(key="statistic", value="value")
  
  return(df_reshape)
}))


slope_rework_allstats$SNAP_nc <- slope_rework_allstats$SNAP + mean(raw.SNAP[,dim], na.rm = TRUE)
slope_rework_allstats$dimension <- dim

toplot <- rbind(toplot, slope_rework_allstats)  
}

pdf("ahat_facet_SNE_errors_10062017.pdf", width = 13, height =7)
a <- ggplot(toplot) + aes(x = SNAP_nc, y=m, ymin = lo, ymax = hi) + geom_line(size = 1.5) + 
  geom_ribbon(aes(color=NULL), alpha=0.15)  + ylab("Predicted threshold (a)") + xlab("SNAP dimension") + theme_bw(base_size=26) + 
  scale_fill_brewer("Group", palette = "Set1") + scale_color_brewer("Group", palette = "Set1") + facet_wrap(~dimension, scales = "free") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
print(a)
dev.off()

pdf("wtfv_vhat_facet_SNE_errors_10062017.pdf", width = 13, height =7)
a <- ggplot(toplot) + aes(x = SNAP_nc, y=m, ymin = lo, ymax = hi, colour = as.factor(slope), fill = as.factor(slope)) + geom_line(size = 1.5) + 
  geom_ribbon(aes(color=NULL), alpha=0.15)  + ylab("Predicted drift rate (v)") + xlab("SNAP dimension") + theme_bw(base_size=26) + 
  scale_fill_brewer("Condition", palette = "Set1", labels = c("Congruent", "Incongruent"))  +scale_color_brewer("Condition", palette = "Set1", labels = c("Congruent", "Incongruent")) + facet_wrap(~dimension, scales = "free") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
print(a)
dev.off()

# pdf("vhat_ixns_facet_SNE.pdf", width = 13, height =7)
# plot <- ggplot(df_all_v, aes(x=qs, y=vhat, color=condition)) + stat_smooth(method="lm") + labs(title ="Drift Rate predicted by SNAP x condition", color = "Condition") +facet_wrap(~dim) +
#   ylab("Predicted drift rate") + xlab("Predcted SNAP quantile") + scale_colour_brewer(palette = "Set1") + theme_bw()
# print(plot)
# dev.off()

df_all_a <- data.frame()
for (dim in SNAP_dim){
  d_df <- dplyr::filter(data, dimension== dim) %>% dplyr::select(-dimension)
  
  snapvar <- SNAP_info %>% group_by(subj_index) %>% 
    filter(row_number() == 1) %>% dplyr::select(subj_index, dim) %>% arrange(subj_index) %>% ungroup()
  qs <- data.frame(qs=c(.2, .5, .8), vals=quantile(snapvar[[dim]], c(.2, .5, .8)))
  
  #cartesian project: replicate posteriors for each of the quantile targets
  d_df <- merge(d_df, qs)
  # d_df_v <- d_df %>% mutate(
  #   Congruent=v_Intercept + v_SNAP*vals,
  #   Incongruent=v_Intercept + v_SNAP*vals + v_C + v_C*vals
  # ) %>% gather(key="condition", value="vhat", Congruent, Incongruent)
  #d_df_v$dim <- rep(dim, length(d_df_v[,1]))
  #df_all_v <- rbind(df_all_v, d_df_v)
  
  d_df_a <- d_df %>% mutate(ahat = a_Intercept + a_SNAP*vals)
  d_df_a$dim <- rep(dim, length(d_df_a[,1]))
  df_all_a <- rbind(df_all_a, d_df_a)
}

pdf("ahat_facet_SNE.pdf", width = 13, height =7)
plot <- ggplot(df_all_a, aes(x=qs, y=ahat)) + stat_smooth(method="lm") + labs(title ="Decision Threshold predicted by SNAP") +facet_wrap(~dim) +
  ylab("Predicted decision threshold") + xlab("Predcted SNAP quantile") + scale_colour_brewer(palette = "Set1") + theme_bw()
print(plot)
dev.off()

pdf("vhat_ixns_SNE.pdf")

  #d_df <- dplyr::filter(data, dimension== dim) %>% dplyr::select(-dimension)
  data_temp <- data
  snapvar <- SNAP_info %>% group_by(subj_index) %>% dplyr::select(subj_index, SNAP_dim)   %>%arrange(subj_index) %>% ungroup()
  qs <- data.frame()
  for (dim in SNAP_dim){
    qs_merge <- data.frame(qs = c(.2,.5, .8), vals = quantile(snapvar[[dim]], c(.2,.5,.8)), dim = rep(dim, 3))
    qs <- rbind(qs, qs_merge)
  }
  
  #cartesian project: replicate posteriors for each of the quantile targets
  data_temp <- merge(data_temp, qs)
  d_df_v <- d_df %>% mutate(
    vhatcon=v_Intercept + v_SNAP*vals,
    vhatincon=v_Intercept + v_SNAP*vals + v_C + v_C*vals
  ) %>% gather(key="condition", value="vhat", vhatcon, vhatincon)
  
  d_df_a <- d_df %>% mutate(ahat = a_Intercept + a_SNAP*vals)
  
  plot <- ggplot(d_df, aes(x=qs, y=vhat, color=condition)) + stat_smooth(method="lm") + labs(title =dim)
  print(plot)

dev.off()

pdf("ahat_SNE.pdf")
for (dim in SNAP_dim){
  d_df <- dplyr::filter(data, dimension== dim) %>% dplyr::select(-dimension)
  snapvar <- SNAP_info %>% group_by(subj_index) %>% 
    filter(row_number() == 1) %>% dplyr::select(subj_index, dim) %>% arrange(subj_index) %>% ungroup()
  qs <- data.frame(qs=c(.2, .5, .8), vals=quantile(snapvar[[dim]], c(.2, .5, .8)))
  
  #cartesian project: replicate posteriors for each of the quantile targets
  d_df <- merge(d_df, qs)
  # d_df_v <- d_df %>% mutate(
  #   vhatcon=v_Intercept + v_SNAP*vals,
  #   vhatincon=v_Intercept + v_SNAP*vals + v_C + v_C*vals
  # ) %>% gather(key="condition", value="vhat", vhatcon, vhatincon)
  
  d_df_a <- d_df %>% mutate(ahat = a_Intercept + a_SNAP*vals)
  
  plot <- ggplot(d_df_a, aes(x=qs, y=ahat)) + stat_smooth(method="lm") + labs(title =dim)
  print(plot)
}
dev.off()




# 
# ###################
# for (d in SNAP_dim){
#   d_df <- dplyr::filter(data, dimension==d) #dplyr equivalent
#   #d_df <- data[which(data$dimension==d),]
#   d_df <- d_df %>% dplyr::select(-dimension)
#   dim_list_traces[[d]] <- d_df
#   B0A <- mean(d_df$a_Intercept)
#   B0A.SE <- StE(d_df$a_Intercept)
#   B1A <- mean(d_df$a_SNAP)
#   B1A.SE <- StE(d_df$a_SNAP)
#   B0V <- mean(d_df$v_Intercept)
#   B0V.SE <- StE(d_df$v_Intercept)
#   B1V <- mean(d_df$v_SNAP)
#   B1V.SE <- StE(d_df$v_SNAP)
#   B2V <- mean(d_df$v_C)
#   B2V.SE <- StE(d_df$v_C)
#   B3V <- mean(d_df$v_ixn)
#   B3V.SE <- StE(d_df$v_ixn)
#   dim_list_descriptives[[paste0(d, "_descriptives")]] <- describe(d_df)
#   dim_list_traces_betas[[paste0(d, "_betas")]] <- data.frame(B0A, B1A, B0V, B1V, B2V, B3V)
#   dim_list_traces_betas_SEs[[paste0(d, "_betas_SEs")]] <- data.frame(B0A.SE, B1A.SE, B0V.SE, B1V.SE, B2V.SE, B3V.SE)
# }
# dim_traces_betaSEs <- do.call(rbind,dim_list_traces_betas_SEs)
# descriptives_df <- do.call(rbind, dim_list_descriptives)
# 
# write.csv(descriptives_df, file = "flanker_models_descriptives_10012017.csv")#, row.names = TRUE)
# 
# 
# 
# #range(SNAP_info$EXHIB)
# 
# #preds <- c("EXHIB", "AGG", "MISTRUST", "MANIP", "DISINH")
# df <- c()
# df_all <- c()
# for (v in SNAP_dim) {
#   SNAP_info.pred <- SNAP_info %>% group_by(subj_index) %>% filter(row_number() == 1) %>% dplyr::select(v) %>% ungroup()
#   SNAP_info.pred <- data.frame(SNAP_info.pred)
#   qs <- quantile(SNAP_info.pred[[v]], c(.3, .5, .8))
#   sumstat <- dim_list_traces_betas[[paste0(v, "_betas")]]
#  
#   vhatcon <- sumstat$B0V + sumstat$B1V*qs
#   vhatincon <- sumstat$B0V + sumstat$B1V*qs + sumstat$B2V + sumstat$B3V*qs
#   
#   vhatall <- sumstat$B0V + sumstat$B1V*qs + sumstat$B2V
#   
#   ahatall <- sumstat$B0A + sumstat$B1A
#   ##only when a varies by condition
#   # ahatcon <- sumstat$B0A + sumstat$B1A*qs
#   # ahatincon <- sumstat$B0A + sumstat$B1A*qs + sumstat$B2A + sumstat$B3A
#   stim <- c(rep("incon",3), rep("con",3))
#  
#  
#   df <- rbind(df, data.frame(pred=v, pred.quant = qs, stimulus = stim, vhat = c(vhatincon, vhatcon), ahat = ahatall))
#   df_all <- rbind(df_all, data.frame(pred=v, pred.quant = qs, vhat = vhatall, ahat = ahatall))                           
# }
# 
# 
# library(ggplot2)
# ggplot(df, aes(x=pred.quant, y = vhat, colour = stimulus)) +geom_line() + geom_point() + facet_wrap(~pred)
# ggplot(df_all, aes(x = pred.quant, y = vhat)) + geom_line()  + facet_wrap(~pred)
# 
# 
