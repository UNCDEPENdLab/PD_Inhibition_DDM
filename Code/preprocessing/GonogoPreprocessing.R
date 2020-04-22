
# cleaning of go-nogo data ------------------------------------------------

## clean RTs for gng task for PD_Inhibition DDM project.
#### from MH dissertation MS:  For this task, participants viewed single letter trials and were required to press a key (i.e., “go”) or to withhold a key press (i.e., “no go”). 
#When a letter other than ‘X’ was presented, participants were instructed to press the space bar as quickly as possible. 
#When an ‘X’ appeared, participants were instructed not to press the space bar (i.e., “no go”). 
#Letters were presented for 500ms followed by a 1500 interstimulus interval (Casey et al., 1997) and 
#the frequency of no go trials was fixed at 20% within a block in order to promote a tendency to perform a key press. 
#Extending the work of Durston and colleagues (2002), the number of go trials preceding a no go trial was parametrically manipulated within each block to create varying 
#levels of difficulty. No go trials were preceded by one, three, five, or seven go trials, and difficulty withholding key presses on no go trials was 
#thought to increase as a function of the number of preceding go trials. ACC, VLPFC, and superior parietal cortex exhibit greater activation as the number of preceding 
#go trials increases (Durston et al., 2002). Trials of varying difficulty were randomized within the trial block. 
#Participants completed three blocks of 64 trials with two 30-second breaks. The frequency of key presses on no go trials was the primary dependent variable.


# setup key variables and load in required packages -----------------------

basedir <- "~/github_repos/PD_Inhibition_DDM/"; setwd(basedir)
datadir <- paste0(basedir, "Data/SAS\ Originals")
# list.files(datadir)

pacman::p_load(sas7bdat)

data <- read.sas7bdat(paste0(basedir, "Data/SAS\ Originals/gonogo.sas7bdat"))


# pull one subject and figure out data structure ---------------------------



ng <- data %>% group_by(Subject, Block) %>% filter(row_number()==1) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("NoGoDisplay")) %>% ungroup()
names(ng) <- sub("NoGoDisplay_", "", names(ng))
ng$TrialType <- "NoGo"
ng$ACC <- ifelse(ng$RESP == "{SPACE}", 0, 1)


go <- data %>% group_by(Subject, Block) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("GoDisplay")) %>% ungroup()
names(go) <- sub("GoDisplay_", "", names(go))
go$TrialType <- "Go"

gng_all <- data.frame()



# interleave go and no-go trials to replicate presentation order ----------


for(sub in unique(data$Subject)){
  df_sub <- go %>% filter(Subject == sub)
  sub_stacked <- data.frame()
  for(bl in unique(df_sub$Block)){
    this.bl <- df_sub %>% filter(Block == bl)
    ng.bl <- ng %>% filter(Subject == sub, Block == bl)
    bl_tot <- rbind(this.bl, ng.bl)
    sub_stacked <- rbind(sub_stacked, bl_tot)
  }
  
  sub_stacked <- sub_stacked %>% select(-Procedure_Block_, -DurationError, -OnsetDelay, -OnsetTime, - RTTime)
  
  
  gng_all <- rbind(gng_all, sub_stacked)
}

gng <- gng_all %>% dplyr::rename(id=Subject, correct=ACC, rt=RT, 
                                  stim=TrialType, cond = Procedure_Trial_, cresp = CRESP, resp = RESP) %>%
  mutate(
    # trial = trial -8,  #1:10 are practice trials
    stim = factor(stim, levels=c("Go", "NoGo")),
    cond = factor(cond, levels = c("GoProc", "ThreeGoProc", "FiveGoProc", "SevenGoProc"), labels = c("OneGo", "ThreeGo", "FiveGo", "SevenGo")),
    stim_cond = as.factor(paste0(stim, "_", cond)),
    trial = rep(1:240, length(unique(gng_all$Subject))), #risky but looks fine when implemented.
    trial_z = as.vector(scale(trial)),
    cresp = if_else(cresp == "{SPACE}", 1,0),
    resp = if_else(resp == "{SPACE}", 1,0),
    rt_inv=-1000/rt #puts Rts on the -4 -- -1 range (makes parameter estimates a bit easier to see)
  ) %>% group_by(id) %>% arrange(id,trial) %>% mutate(
  prev_rt=dplyr::lag(rt, 1)/1000, #convert to seconds to put it on a smaller scale, roughly comparable to the -1000/rt rt_inv
  prev_rt_inv=dplyr::lag(rt_inv, 1)
) %>% arrange(id, trial) %>% dplyr::group_by(cond) %>%
  dplyr::mutate(block_trial = dplyr::row_number()) %>% #in case we want to allow for the potential for the number of trials since the last nogo to influence RT rather than simply blocking by "goproc" or adding overall trial number (blind to sequencing of go and nogo trials) as covariate
  ungroup()


#sanity checks on structure
xtabs(~cond +id, gng) # 12 of each condition.
xtabs(~stim +cond , gng)


# sanity checks -----------------------------------------------------------


# A few super-short RTs, and then a long RT tail

ggplot(gng, aes(x=stim_cond, y=rt)) + geom_boxplot() + theme_cowplot() + theme(axis.text.x = element_text(angle = 45))



accuracies <- gng %>% group_by(id, stim_cond) %>%
  dplyr::summarise(nacc=sum(correct==1, na.rm=T), nerr=sum(correct==0, na.rm=T)) %>%
  pivot_wider(names_from="stim_cond", values_from=c("nerr", "nacc")) %>%
  mutate(
    p_acc=(nacc_Go_OneGo + nacc_NoGo_OneGo + nacc_Go_ThreeGo + nacc_NoGo_ThreeGo + nacc_Go_FiveGo + nacc_NoGo_FiveGo + nacc_Go_SevenGo + nacc_NoGo_SevenGo)/(nacc_Go_OneGo + nacc_NoGo_OneGo + nacc_Go_ThreeGo + nacc_NoGo_ThreeGo + nacc_Go_FiveGo + nacc_NoGo_FiveGo + nacc_Go_SevenGo + nacc_NoGo_SevenGo + nerr_Go_OneGo + nerr_NoGo_OneGo + nerr_Go_ThreeGo + nerr_NoGo_ThreeGo + nerr_Go_FiveGo + nerr_NoGo_FiveGo + nerr_Go_SevenGo + nerr_NoGo_SevenGo),
    p_acc_Go_OneGo=(nacc_Go_OneGo)/(nacc_Go_OneGo + nerr_Go_OneGo),
    p_acc_NoGo_OneGo=(nacc_NoGo_OneGo)/(nacc_NoGo_OneGo + nerr_NoGo_OneGo),
    p_acc_Go_ThreeGo=(nacc_Go_ThreeGo)/(nacc_Go_ThreeGo + nerr_Go_ThreeGo),
    p_acc_NoGo_ThreeGo=(nacc_NoGo_ThreeGo)/(nacc_NoGo_ThreeGo + nerr_NoGo_ThreeGo),
    p_acc_Go_FiveGo=(nacc_Go_FiveGo)/(nacc_Go_FiveGo + nerr_Go_FiveGo),
    p_acc_NoGo_FiveGo=(nacc_NoGo_FiveGo)/(nacc_NoGo_FiveGo + nerr_NoGo_FiveGo),
    p_acc_Go_SevenGo=(nacc_Go_SevenGo)/(nacc_Go_SevenGo + nerr_Go_SevenGo),
    p_acc_NoGo_SevenGo=(nacc_NoGo_SevenGo)/(nacc_NoGo_SevenGo + nerr_NoGo_SevenGo),
  )

g1 <- ggplot(accuracies, aes(y=p_acc)) + geom_boxplot() + ggtitle("Overall Accuracy")
g2 <- ggplot(accuracies, aes(y=p_acc_Go_OneGo)) + geom_boxplot() + ggtitle("Go_OneGo Accuracy")
g3 <- ggplot(accuracies, aes(y=p_acc_NoGo_OneGo)) + geom_boxplot() + ggtitle("NoGo_OneGo Accuracy")
g4 <- ggplot(accuracies, aes(y=p_acc_Go_ThreeGo)) + geom_boxplot() + ggtitle("Go_ThreeGo Accuracy")
g5 <- ggplot(accuracies, aes(y=p_acc_NoGo_ThreeGo)) + geom_boxplot() + ggtitle("NoGo_ThreeGo Accuracy")
g6 <- ggplot(accuracies, aes(y=p_acc_Go_FiveGo)) + geom_boxplot() + ggtitle("Go_FiveGo Accuracy")
g7 <- ggplot(accuracies, aes(y=p_acc_NoGo_FiveGo)) + geom_boxplot() + ggtitle("NoGo_FiveGo Accuracy")
g8 <- ggplot(accuracies, aes(y=p_acc_Go_SevenGo)) + geom_boxplot() + ggtitle("Go_SevenGo Accuracy")
g9 <- ggplot(accuracies, aes(y=p_acc_NoGo_SevenGo)) + geom_boxplot() + ggtitle("NoGo_SevenGo Accuracy")

plot_grid(g1, g2, g4,g6, g8, nrow=1)
plot_grid(g3, g5, g7, g9, nrow=1)
#exclusion subjects from dissertation
accuracies %>% filter(id %in% c(7, 14, 36, 38, 47, 55, 62, 83)) %>% select(id, starts_with("p_acc"))

accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_Go_OneGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_NoGo_OneGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_Go_ThreeGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_NoGo_ThreeGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_Go_FiveGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_NoGo_FiveGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_Go_SevenGo)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_NoGo_SevenGo)


## Subject exclusions

### Exclusion coding

# Generate a column called exclude_gng with the following scheme:
#   
# - 0: a good subject, no concerns
# - 1: a subject with minor problems that we should try to salvage
# - 2: a suspicious subject whose behavior is questionable, but not obviously bad
# - 3: a subject who should certainly be excluded.
# 
# Following this scheme, here are the combined recommendations of Nate and Michael. Here, I've built on the code above, plus Nate's gngPreprocessing R Markdown.

#### 1: questionable, but not clearly flawed

# - 15: Is considerably lower in accuracy than the rest of the group on go trials, however they perform worse as the length of the blocks decreases making me think these are true errors (e.g. experimental manipulation works in the proper direction)
# - 7: This person conversely is quite low in accuracy on all of the no-go trials, however their go performance is close to 1 making me think they are perhaps just strongly biased towards go.


gng <- gng %>% mutate(
  exclude_go_nogo=dplyr::recode(id,
                                      `15` = 1, `7` = 1,
                                      .default = 0
  ))


# RT problems -------------------------------------------------------------


gng %>% group_by(stim_cond) %>% dplyr::summarise(q5 = quantile(rt, .05), q1 = quantile(rt, .025))

ggplot(data = filter(gng, resp == 1), aes(x = rt)) + geom_histogram() + facet_wrap(~stim_cond, scales = "free")

#Number of trials with RTs < 200.

nrow(gng %>% filter(rt < 200 & rt != 0))
nrow(gng %>% filter(rt < 200 & rt != 0))/nrow(gng)
# 116 RTs. This is about 0.4% of the data, which is a very minor loss. Increasing to 250 ms includes ~650 more RTs, which I feel motivated to keep in the data

#Set accuracy and RT to NA when RT < 200ms.


gng <- gng %>% mutate(rt=if_else(rt < 200, NA_real_, rt), correct=if_else(rt < 200, NA_real_, correct))

# Winsorizing long RTs:


# What about Winsorizing the top 1% of inv-transformed RTs by condition and block?
test <- gng %>% filter(resp == 1) %>% #only look at trials where subjects emmitted go response
  group_by(stim_cond) %>% 
  mutate(
    rt_log=log10(rt), 
    rt_log_winsor=DescTools::Winsorize(rt_log, probs=c(0, .99), na.rm = T),
    rt_winsor=DescTools::Winsorize(rt, probs=c(0, .99), na.rm = T) 
  ) %>%
  ungroup() %>% group_by(id, stim_cond) %>% #use per-subject condition and block RT quantiles to trim
  mutate(
    rt_inv_trim=if_else(rt_inv < quantile(rt_inv, 0.01, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_trim=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_log_trim=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup() %>%
  group_by(stim_cond) %>% #use overall condition and block RT quantiles to trim
  mutate(
    rt_inv_trim_grp=if_else(rt_inv < quantile(rt_inv, 0.01, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_trim_grp=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_log_trim_grp=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup()


#the problem with per-subject trimming is that many plausible values are removed since they are unlikely compared to the subject's distribution
test %>% filter(is.na(rt_inv_trim)) %>% pull(rt) %>% hist(main="RTs dropped by subject-specific trim")
test %>% filter(is.na(rt_inv_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")
test %>% filter(is.na(rt_inv_trim_grp)) %>% nrow()
test %>% filter(is.na(rt)) %>% nrow()
test %>% filter(is.na(rt_trim_grp)) %>% nrow()
test %>% filter(is.na(rt_winsor)) %>% nrow()
test %>% filter(is.na(rt_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")

# raw RT distributions
trim <- ggplot(data = test, aes(x = rt_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = test, aes(x = rt_trim_grp)) + geom_density()
winsor <- ggplot(data = test, aes(x = rt_winsor)) + geom_histogram()
winsor_d <- ggplot(data = test, aes(x = rt_winsor)) + geom_density()

sub <- ggplot(data = test, aes(x = rt_trim)) + geom_histogram()
sub_d <- ggplot(data = test, aes(x = rt_trim)) + geom_density()

raw <- ggplot(data = test, aes(x = rt)) + geom_histogram()
raw_d <- ggplot(data = test, aes(x = rt)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, winsor, winsor_d, sub, sub_d, ncol = 2)

# Inverse RT distributions
trim <- ggplot(data = test, aes(x = rt_inv_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = test, aes(x = rt_inv_trim_grp)) + geom_density()
# winsor <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_histogram()
# winsor_d <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_density()

sub <- ggplot(data = test, aes(x = rt_inv_trim)) + geom_histogram()
sub_d <- ggplot(data = test, aes(x = rt_inv_trim)) + geom_density()

raw <- ggplot(data = test, aes(x = rt_inv)) + geom_histogram()
raw_d <- ggplot(data = test, aes(x = rt_inv)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, sub, sub_d, ncol = 2)

# log RT distributions
trim <- ggplot(data = test, aes(x = rt_log_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = test, aes(x = rt_log_trim_grp)) + geom_density()
# winsor <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_histogram()
# winsor_d <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_density()

sub <- ggplot(data = test, aes(x = rt_log_trim)) + geom_histogram()
sub_d <- ggplot(data = test, aes(x = rt_log_trim)) + geom_density()

raw <- ggplot(data = test, aes(x = rt_log)) + geom_histogram()
raw_d <- ggplot(data = test, aes(x = rt_log)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, sub, sub_d, ncol = 2)

g1 <- ggplot(test, aes(sample=rt_log_winsor)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q Winsorize log(RT)")
g2 <- ggplot(test, aes(sample=rt_log_trim)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q Trimmed log(RT)")
g3 <- ggplot(test, aes(sample=rt_inv_trim)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q Trimmed 1/RT")
g4 <- ggplot(test, aes(sample=rt_log_trim_grp)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q Group-Trimmed log(RT)")
g5 <- ggplot(test, aes(sample=rt_inv_trim_grp)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q Group-Trimmed 1/RT")
g6 <- ggplot(test, aes(sample=rt_inv)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
  facet_wrap(~stim) + ggtitle("Normal Q-Q 1/RT")

plot_grid(g1, g2, g3, g4, g5, g6, nrow=3)

## stick with group-trimmed log-RT, doesnt quite get nogo trials but they are not of imminent interest.

### look for high-NA subjects and unusually high mean RTs for potential exclusion

rt_sums <- test %>% group_by(id) %>% dplyr::summarise(mean_rt = harmonic.mean(rt_log_trim_grp), p_na = length(is.na(rt_log_trim_grp)))
hist(rt_sums$mean_rt) # looks fine to me
# View(rt_sums)
# View(test)
# test

nas <- test[which(is.na(test$rt_log_trim_grp)),]
# hist(nas$id)
table(nas$id) %>% sort() #%>% hist()
# yeah, lets get 15 out of here, they have waay more NAs than anyone else, clearly making them an outlier.
gng <- gng %>% mutate(
  exclude_go_nogo=dplyr::recode(id,
                                `15` = 3, 
                                `7` = 1,
                                .default = 0
  ))

# good news is that subjects with a high number of RAs have already been flagged above, I dont see the need to flag anyone else as problematic.
gng <- test %>% select(id, Block, stim_cond, trial, rt_inv_trim_grp) %>% right_join(gng, by = c("id", "Block", "stim_cond", "trial")) %>% 
  mutate(rt_inv_trim_grp = ifelse(stim == "NoGo", rt, rt_inv_trim_grp))

# RT distributions for clean and full sample  -----------------------------

gng_SNAP <- gng %>% left_join(select(SNAP_all_scored, id, exclude_SNAP))

gng_full <- dplyr::filter(gng_SNAP, exclude_SNAP ==0, exclude_go_nogo != 3, stim == "Go")
gng_mean_ci_conditions_full <- summarySEwithin(data = gng_full, measurevar = "rt", withinvars = "cond", idvar = "id", na.rm = TRUE)
full_n  <- length(unique(gng_full$id))

gng_point_full <- ggplot(gng_mean_ci_conditions_full, aes(x = rt_norm, y = cond)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt_norm-ci, xmax=rt_norm+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("gng task full sample: N=", full_n))


gng_clean <- dplyr::filter(gng_SNAP, exclude_SNAP ==0, exclude_go_nogo ==0)
gng_mean_ci_conditions_clean <- summarySEwithin(data = gng_clean, measurevar = "rt", withinvars = "cond", idvar = "id", na.rm = TRUE)

clean_n <- length(unique(gng_clean$id))
gng_point_clean <- ggplot(gng_mean_ci_conditions_clean, aes(x = rt_norm, y = cond)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt_norm-ci, xmax=rt_norm+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("gng task clean sample: N=", clean_n))

plot_grid(gng_point_clean, gng_point_full, ncol = 1)


# 
# 
# # display accuracies per subject split by condition -----------------------
# 
# gng_accs <- gng %>% group_by(id, cond, stim) %>% dplyr::summarise(as.numeric(table(correct)[2])/(as.numeric(table(correct)[2])+as.numeric(table(correct)[1])))
# # colnames(gng_accs)[3] <- "ACC"
# 
# gng_accs %>% arrange(ACC)
# 
# 
# gng_accs_proc <- gng_all %>% group_by(Subject, TrialType, Procedure_Trial_) %>% dplyr::summarise(as.numeric(table(ACC)[2])/(as.numeric(table(ACC)[2])+as.numeric(table(ACC)[1])))
# colnames(gng_accs_proc)[4] <- "ACC"
# 
# gng_accs_proc %>% arrange(ACC) 
# 
# 
# # drop 7: 33.3% on nogo trials is substantially lower than the rest. 
# 
# gng_all$exclude <- ifelse(gng_all$Subject == 7,1,0)
#   
# # look at RTs across go trials --------------------------------------------
# 
# go <- gng_all %>% filter(TrialType == "Go") #%>% hist(.$RT)
# 
# hist(go$RT)
# quantile(go$RT, seq(.1,1,.1))
# 
# ##according to MH diss the letters were only presented for 500 ms and therefore, anything beyond that should be investigated.
# 
# go %>% arrange(-RT)
# 
# 
# # scratch below -----------------------------------------------------------
# 
# 
# # sample_sub <- data %>%  filter(Subject == 1)
# # 
# # attr(sample_sub, "column.info") <- NULL
# # str(sample_sub)
# # 
# # nrow(sample_sub) #should equal 64*3 = 192
# # unique(table(data$Subject)) #perfect.
# # 
# # head(sample_sub)
# # 
# # for(i in colnames(sample_sub)){
# #   cat(i, "\n")
# #   print(table(sample_sub[,i]))
# # }
# # 
# # # ## not entirely intuitive, maybe pull one block of 64 trials and see what we can glean from its structure.
# # # 
# # # b1 <- sample_sub[1:64,]
# # # head(b1, 32)
# # # 
# # # ## still strange.
# # # 
# # # data_scores <- read.sas7bdat(paste0(datadir, "/gonogoscores.sas7bdat"))
# # # hist(data_scores$T_MeanGoRT)
# # 
# # 
# # 
# # data <- read.sas7bdat(paste0(datadir, "/gonogo.sas7bdat"))
# # ng <- data %>% group_by(Subject, Block) %>% filter(row_number()==1) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("NoGoDisplay")) %>% ungroup()
# # names(ng) <- sub("NoGoDisplay_", "", names(ng))
# # ng$TrialType <- "NoGo"
# # go <- data %>% group_by(Subject, Block) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("GoDisplay")) %>% ungroup()
# # names(go) <- sub("GoDisplay_", "", names(go))
# # go$TrialType <- "Go"
# # # 
# # # 
# # # ## questions for MH
# # # # like, so many. Better to just talk in person.
# # # 
# # # godisplay resp = space. correct response
# # # godisplay resp = empty. ommision error
# # # nogodisplay 
# # # 
# # # 
# # # 
# # # 
# # # 
