
# Screening on SNAP, MPQ, and behavioral indices --------------------------


basedir <- "~/ics/Nate/PD_Inhibition_DDM/"; setwd(basedir)
datadir <- paste0(basedir, "Data/")
figuredir <- paste0(basedir,"Figures/")

alldat <- read.csv(paste0(datadir, "alldat.csv"))


self_reps <- alldat %>% select(Subject, starts_with("Z_"), starts_with("ATQ"), starts_with("MPS_"))


length(unique(alldat$Subject))
## looks like with full data we only have 108 to start with. This means we'll have to rerun all the descriptive checks (gender, age, race, etc)


self_reps %>% na.omit() %>% nrow()
# View(self_reps)

# Subject 7 missing MPQ
# Subject 38 missing SNAP
# Full data for ATQ
############
# 106 with complete personality data
############

hist(self_reps$Z_InvalidityIndex)
quantile(self_reps$Z_InvalidityIndex, na.rm = TRUE)
# looks pretty normal, nothing goes above a 2 though... I wonder if these are due to the folks that Michael dropped for his diss.
# according to his records 2, 23, 51, and 109 were dropped bc of high SNAP invalidity index 
# 48, 62, 82, 85 were dropped due to MPQ invalidity
unique(alldat$Subject)
# okay, yup. Those are the folks he cut. Can go back to the raw data and see how egregious they were and maybe we can keep them in. 
# all the bad MPQ folks still remain, this is likely bc they did not get included in the diss manuscript. 
# Perhaps we touch base with Tim to see if he thinks we are warranted in dropping MPQ from analyses or if they add info that will be useful for our purposed. 
# My prior going in to such a conversation will be to drop MPQ, so as not to overpopulate the personality trait space. 
# We could try creating our own factors, but that may not be worth it for the scope of this project. 
hist(alldat$MPS_inv)



# who are the folks that got dropped in my behavioral cleaning scripts? --------------

bad_SNAPs <- c(2, 23, 38, 51, 109) #include 38 for having NA responses
(flanker_rts <- get(load(paste0(datadir, "cache/rt_drop_flanker.RData"))))
(flanker_accs <- get(load(paste0(datadir, "cache/acc_drop_flanker.RData"))))
(rp_accs <- get(load(paste0(datadir, "cache/acc_drop_recent_probes.RData"))))
(rp_nas <- get(load(paste0(datadir, "cache/na_drop_recent_probes.RData"))))

bad_behaviorals <- unique(c(flanker_rts, flanker_accs, rp_accs, rp_nas))
length(bad_behaviorals)

(messy_subs <- unique(c(bad_SNAPs, bad_behaviorals)) %>% sort())
112- length(messy_subs)

#3/24/20: this leaves 93 subjects with purportedly valid SNAP profiles and solid behavioral data on flanker and recent probes
