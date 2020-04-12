
# compile item-level MPQ data------------------------------------


pacman::p_load(tidyverse, sas7bdat, assertr)

basedir <- "~/github_repos/PD_Inhibition_DDM"; setwd(basedir)

mpq_1 <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq_1.sas7bdat")) %>% select(id, starts_with("MPQ_")) %>% dplyr::rename(subject = id); nrow(mpq_1)
mpq_2 <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq_2.sas7bdat"))%>% select(id, starts_with("MPQ_")) %>% dplyr::rename(subject = id); nrow(mpq_2)
mpq <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpq.sas7bdat")); nrow(mpq)
mpq_scores <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpqscores.sas7bdat")) %>% dplyr::rename(subject = SUBID); nrow(mpq_scores)

#no subject names a little strange for MPQ and that mpq1 has 12 respondants

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
