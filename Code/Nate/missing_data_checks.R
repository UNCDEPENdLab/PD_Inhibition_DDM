
pacman::p_load(sas7bdat, tidyverse)

SNAP <- read.sas7bdat("~/Downloads/snap2scores.sas7bdat")
# nrow(SNAP)
# tail(SNAP)
basedir <- "~/github_repos/PD_Inhibition_DDM/"; setwd(basedir)




# missing data checks -----------------------------------------------------

### quick confirmation that these subjects do/don't have task data
missing_snap <- data.frame(subject = which(!1:114 %in% SNAP$subject), SNAP = 1) 

mpq <- read.sas7bdat(file.path(basedir,"MH_Diss_Analyses/Data/mpqscores.sas7bdat")) %>% dplyr::rename(subject = SUBID)
missing_mpq <- data.frame(subject = which(!1:114 %in% mpq$subject), mpq = 1)

flank <- read.sas7bdat(file.path(basedir, "Data/SAS Originals/flanker.sas7bdat"))
missing_flank <- data.frame(subject = which(!1:114 %in% flank$Subject), flank = 1)

rp <- read.sas7bdat(file.path(basedir, "Data/SAS Originals/recentprobes.sas7bdat"))
missing_rp <- data.frame(subject = which(!1:114 %in% rp$Subject), rp = 1)

gng <- read.sas7bdat(file.path(basedir, "Data/SAS Originals/gonogo.sas7bdat"))
missing_gng <- data.frame(subject = which(!1:114 %in% gng$Subject), gng = 1)

all_missing <- missing_snap %>% full_join(missing_flank, by = "subject") %>% full_join(missing_rp, by = "subject") %>% full_join(missing_gng, by = "subject") %>% full_join(missing_mpq, by = "subject")

print(all_missing)
