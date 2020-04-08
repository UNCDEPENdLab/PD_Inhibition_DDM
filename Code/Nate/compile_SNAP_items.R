
# compile item-level SNAP data for Tim ------------------------------------


pacman::p_load(tidyverse, sas7bdat)
basedir <- "~/github_repos/PD_Inhibition_DDM"; setwd(basedir)

if(file.exists(file.path(basedir,"Data/R_Originals/SNAP_all.csv"))){
  SNAP_all <- read.csv(file.path(basedir,"Data/R_Originals/SNAP_all.csv"))
} else {
  #to rerun, swap dir calls
  
  SNAP_all <- read.sas7bdat(file.path(basedir, "MH_Diss_Analyses/Data/snap2.sas7bdat")) %>% select(subject, starts_with("SNAP_"))
  
  # subject numbers in 1-4 do not line up with snap2.sas7bat. Similar experience with MPQ.
  
  # snap_1 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_1.sas7bdat")
  # snap_2 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_2.sas7bdat")
  # snap_3 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_3.sas7bdat")
  # snap_4 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_4.sas7bdat")
  # 
  # SNAP_all <-  snap_1 %>% select(-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1)
  # SNAP_all <- left_join(SNAP_all, select(snap_2,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")
  # SNAP_all <- left_join(SNAP_all, select(snap_3,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")
  # SNAP_all <- left_join(SNAP_all, select(snap_4,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")
  
  write.csv(SNAP_all, file = "~/ics/Nate/PD_Inhibition_DDM/Data/R_Originals/SNAP_all.csv",row.names = FALSE)
  
}

# names(SNAP_all)

output <- percentage_checks(SNAP_all, "subject") %>% select(subj, perc_sums) %>% rename(subject = subj) %>% mutate(perc_sum_rank = 1:nrow(.))

head(output)

# assess validity ---------------------------------------------------------


SNAP_all_scored <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/snap2scores.sas7bdat"))
invalid_plot <- ggplot(SNAP_all_scored, aes(x = Z_InvalidityIndex)) + geom_histogram()
print(invalid_plot)

validity_indices <- c("Z_VRIN", "Z_TRIN","Z_InvalidityIndex")#, "Z_DRIN", "Z_RareVirtues", "Z_Deviance",  "Z_BackDeviance")

vdf <- SNAP_all_scored %>% select(subject, validity_indices)


vmelt <- reshape2::melt(vdf, id.vars = "subject" )
v_plot <- ggplot(vmelt, aes(x = value)) + geom_histogram() + facet_wrap(~variable, scales = "free")
print(v_plot)



vdf <- vdf %>% left_join(output, by = "subject")
cor(vdf)

vdf %>% arrange(-abs(Z_InvalidityIndex)) %>% head(10)
vdf %>% arrange(-abs(Z_TRIN)) %>% head(10)
vdf %>% arrange(-abs(Z_VRIN)) %>% head(10)

# From Patrick et al (2002) MPQ-BF paper: a) VRIN above 3, b) TRIN above 3.21, c) VRIN above 2 and TRIN above 2.28
# I'm going to vote for dropping the two subjects that satisfy condition c) (e.g 51 and 100), and the one subject that had an invalidity index above 3.16 bc they seem a little whack (e.g.23).

SNAP_all_scored$exclude_SNAP <- ifelse(SNAP_all_scored$subject %in% c(23,51,100), 1, 0); save(SNAP_all_scored, file = file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData"))
SNAP_all$exclude_SNAP <- ifelse(SNAP_all$subject %in% c(23,51,100), 1, 0); save(SNAP_all_scored, file = file.path(basedir, "Data/preprocessed/SNAP_all_item_level_final.RData"))
