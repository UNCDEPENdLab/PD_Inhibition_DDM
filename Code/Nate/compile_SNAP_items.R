
# compile item-level SNAP data for Tim ------------------------------------


pacman::p_load(tidyverse, sas7bdat, assertr, dependlab)
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

R.utils::sourceDirectory(file.path(basedir, "Code/Functions"))

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
vdf %>% arrange(-Z_VRIN) %>% head(10)

# From Patrick et al (2002) MPQ-BF paper: a) VRIN above 3, b) TRIN above 3.21, c) VRIN above 2 and TRIN above 2.28
# I'm going to vote for dropping the two subjects that satisfy condition c) (e.g 51 and 100),
# and the one subject that had an invalidity index above 3.16 bc their scores on this scale, which is derived via all other invalidity scales is markedly different than the rest of the distribution (e.g.23).


SNAP_all_scored$exclude_SNAP <- ifelse(SNAP_all_scored$subject %in% c(23,51,100), 1, 0); 

save(SNAP_all_scored, file = file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData"))
SNAP_all$exclude_SNAP <- ifelse(SNAP_all$subject %in% c(23,51,100), 1, 0); 
save(SNAP_all_scored, file = file.path(basedir, "Data/preprocessed/SNAP_all_item_level_final.RData"))



# MPQ validity checks ------------------------------------------------------------


# compile item-level MPQ data------------------------------------


mpq_1 <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq_1.sas7bdat")) %>% select(id, starts_with("MPQ_")) %>% dplyr::rename(subject = id); nrow(mpq_1)
mpq_2 <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq_2.sas7bdat"))%>% select(id, starts_with("MPQ_")) %>% dplyr::rename(subject = id); nrow(mpq_2)
mpq <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq.sas7bdat")); nrow(mpq)
mpq_scores <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpqscores.sas7bdat")) %>% dplyr::rename(subject = SUBID); nrow(mpq_scores)

#no subject names a little strange for MPQ and that mpq1 has 12 respondents

# lets see who's missing in each data set
missing_mpq <- data.frame(subject = which(!1:114 %in% mpq$subject), mpq = 1)
missing_mpq_scores <- data.frame(subject = which(!1:114 %in% mpq_scores$subject), mpq_scores = 1)
missing_mpq_1 <- data.frame(subject = which(!1:114 %in% mpq_1$subject), mpq_1 = 1)
missing_mpq_2 <- data.frame(subject = which(!1:114 %in% mpq_2$subject), mpq_2 = 1)

all_missing_mpq <- missing_mpq %>% full_join(missing_mpq_scores, by = "subject") %>% full_join(missing_mpq_1, by = "subject") %>% full_join(missing_mpq_2, by = "subject")

# so... not entirely sure how item-level or scored MPQ data was salvaged given the missingness of some subjects in the mpq1 and 2 data sets, but for distance checks lets use just the MPQ item-level data that are contained in the mpq df.


mpq <- mpq %>% select(subject,starts_with("MPQ_"))

mpq_percs <- percentage_checks(mpq,inds = "subject", plot_hist = TRUE) %>% rename(subject = subj) %>% select(subject, perc_sums)

mpq_percs %>% arrange(perc_sums)


# datadir <- paste0(basedir, "/Data/")
# alldat <- read.csv(paste0(datadir, "alldat.csv"))
# names(alldat)



# calculate VRIN and TRIN -------------------------------------------------
mpq_valids <- calc_MPQSF_VRIN_TRIN(mpq) %>% select(subject, starts_with("Z_"))
mpq_valids



mpq_valids %>% arrange(-abs(Z_TRIN)) %>% head(10)
mpq_valids %>% arrange(-Z_VRIN) %>% head(10)

v <- ggplot(mpq_valids, aes(x = Z_VRIN)) + geom_histogram(bins = 15)
t <- ggplot(mpq_valids, aes(x = Z_TRIN)) + geom_histogram(bins = 15)

cowplot::plot_grid(v,t)

# From Patrick et al (2002) MPQ-BF paper: a) VRIN above 3, b) TRIN above 3.21, c) VRIN above 2 and TRIN above 2.28
# I'm going to vote for dropping the two subjects that satisfy conditions a) and b) (e.g 48, 82, 85),


### to Tim's point earlier: if SNAP and MPQ were administered at the same time, we should hope that invalidity indices should be correlated:
all_valids <- SNAP_all_scored %>% select(subject, Z_InvalidityIndex, Z_TRIN, Z_VRIN) %>% dplyr::rename(Z_II_SNAP = Z_InvalidityIndex,
                                                                                       Z_TRIN_SNAP = Z_TRIN,
                                                                                       Z_VRIN_SNAP = Z_VRIN) %>% left_join(mpq_valids, by = "subject")
cor_heatmap(all_valids) # haha, not quite.. still I don’t think that really calls the trimming procedures into much doubt, at least I hope not…


#to keep naming consistent
MPQ_all_scored <- mpq_scores %>% left_join(mpq_valids, by = "subject") %>% mutate(exclude_MPQ = ifelse(subject %in% c(48,82,85), 1,0))
save(MPQ_all_scored, file = file.path(basedir, "Data/preprocessed/MPQ_all_scored_final.RData"))
MPQ_all <-  mpq  %>% left_join(mpq_valids, by = "subject") %>% mutate(exclude_MPQ = ifelse(subject %in% c(23,51,100), 1, 0)) 
save(MPQ_all_scored, file = file.path(basedir, "Data/preprocessed/MPQ_all_item_level_final.RData"))

