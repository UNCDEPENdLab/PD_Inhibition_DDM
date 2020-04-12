
# HDDM-personality simple analyses -----------------------------------------------

basedir <- "~/ics/Nate/PD_Inhibition_DDM/"; setwd(basedir)
task <- "recent_probes"
full_models <- TRUE; full <- ifelse(full_models, "full", "practice")
create_plots <- TRUE
pacman::p_load(tidyverse, sas7bdat, skimr, dependlab, brmstools, brms)


SNAP <- read.csv(paste0(basedir, "Data/alldat.csv")) %>% select(Subject, starts_with("Z_")) %>% select(-Z_VRIN, -Z_TRIN, -Z_DRIN, -Z_RareVirtues, -Z_Deviance, -Z_InvalidityIndex, -Z_BackDeviance)
load(paste0(basedir, "/Data/cache/", task, "_subj_sufficient_stats.RData"))

SNAP_dims <- names(SNAP)[-1]


# names(subj_suff_stats)
all_models <- data.frame()
for(p in names(subj_suff_stats)){
  df <- subj_suff_stats[[p]]
  colnames(df)[1] <- "Subject"
  SNAP <- SNAP %>% left_join(df, by = "Subject")
  some_models <- data.frame()
  for(dim in SNAP_dims){
    form <- formula(paste0(dim, "~", p,"_map" ))
    output <- lm(data = SNAP, formula = form)
    summ <- broom::tidy(output) %>% mutate(SNAP_dim  = dim, parameter = p)
    some_models <- rbind(some_models, summ)
  }
  all_models <- rbind(all_models, some_models)
}

all_models %>% arrange(p.value) %>% filter(term != "(Intercept)")
cor_heatmap(SNAP)


v_maps <- SNAP %>% select(starts_with("v_")) %>% select(ends_with("_map"))
cor(v_maps, use = "complete.obs")

hddm_maps <- SNAP %>% select(ends_with("_map"))

some_models <- data.frame()
for(dim in SNAP_dims){
  form <- formula(paste0(dim, "~", paste(names(hddm_maps), collapse = " + ")))
  output <- lm(data = SNAP, formula = form)
  summ <- broom::tidy(output) %>% mutate(SNAP_dim  = dim)#, parameter = p)
  some_models <- rbind(some_models, summ)
}

some_models %>% filter(term != "(Intercept)") %>% arrange(p.value)


# devtools::install_github("mvuorre/brmstools")
# 
# # brm()

pdf(pdf(file = paste0(basedir,"Figures/subj_point_ests_sd_recent_probes.pdf"), width =11, height =8))
for(i in sub("_map", "",names(hddm_maps))){
  x <- ggplot(SNAP, aes(x=paste0(i,"_map"), y=Subject)) +
    geom_segment(aes(x = paste0(i,"_map")-paste0(i,"_sd")*2, xend = paste0(i,"_map")+paste0(i,"_sd")*2, y=Subject, yend=Subject)) +
    geom_point()
  print(x)
}
dev.off()