## scratchy code for looking at factoring parameters

# devtools::install_github("PennStateDEPENdLab/dependlab")
# 
# devtools::install_github("RobinHankin/gsl")
# 
# install.packages("glue")
pacman::p_load(tidyverse,psych, psychTools, dependlab)




# x <- get(load("~/github_repos/PD_Inhibition_DDM/Outputs/Data/cache/flanker_full_subset_subj_sufficient_stats.RData"))
x <- get(load("/Users/natehall/github_repos/PD_Inhibition_DDM/Data/cache/flanker_full_subset_subj_sufficient_stats.RData"))
do.call(cbind,x) 
flanker <- data.frame(id = x[["a"]][,1])
for(n in names(x)){
  flanker <- cbind(flanker, x[[n]][,c(2,3)])
}

flanker <- flanker %>% select(id, ends_with("_map"))
names(flanker)[2:7] <- paste0(names(flanker)[2:7], "_flank")
str(flanker)

load("~/github_repos/PD_Inhibition_DDM/Outputs/posterior_summaries/samp20000/recent_probes/clean_sample/summary_stats_vst_reg.RData")
summary_stats

df <- summary_stats[["subjects"]] %>% select(subject, Parameter, MAP) 

df <- df %>% spread(key = Parameter, value = MAP) %>% arrange(as.numeric(subject)) %>% rename(id = subject) #%>% select(-subject)
df$id <- as.numeric(df$id)
df <- df %>% full_join(flanker, by = "id")


df_scaled <- data.frame(scale(df))


pairs.panels(df_scaled)


cor_heatmap(df_scaled)



data("sat.act")
  