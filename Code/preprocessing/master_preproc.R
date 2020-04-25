
# Task cleaning: Master ---------------------------------------------------

### combine MNH and NTH cleaning procedures for a final take on who to drop and who to put in the group of folks that will be removed from analysis for sensitivity analysis
## also retains MNH's variable naming scheme, will be good to keep this the same across analyses.

basedir <- "~/github_repos/PD_Inhibition_DDM/"

if (!require(pacman)) { install.packages("pacman"); library(pacman) }
p_load(car, brms, nlme, lme4, loo, readr,tidyverse, emmeans, cowplot, glmmTMB, bbmle, broom, psych)


load(file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData")) 
SNAP_all_scored <- tibble(SNAP_all_scored) %>% dplyr::rename(id = subject)# load SNAP data. 


R.utils::sourceDirectory(paste0(basedir, "Code/Functions"))

# Data cleaning and setup: Flanker -------------------------------------------------


#raw flanker data, all subjects
load(paste0(basedir, "Data/R_Originals/flanker_raw.RData"))

flanker <- flanker %>%
  select(Subject, TrialSlide_ACC, TrialSlide_RT, Incongruent, CongruentBlock, Block) %>%
  dplyr::rename(id=Subject, correct=TrialSlide_ACC, rt=TrialSlide_RT, 
                block=CongruentBlock, trial=Block, cond=Incongruent) %>%
  mutate(
    trial = trial - 8,  #1:8 are practice trials
    block = factor(block, levels=c(0, 1), labels=c("most_incon", "most_con")),
    cond = factor(cond, levels=c(0,1), labels=c("congruent", "incongruent")),
    block_number = case_when(
      trial <= 40 ~ 1,
      trial > 40 & trial <= 80 ~ 2,
      trial > 80 & trial <= 120 ~ 3,
      trial > 120 & trial <= 160 ~ 4
    ),
    trial_z = as.vector(scale(trial)),
    rt_inv=-1000/rt, #puts Rts on the -4 -- -1 range (makes parameter estimates a bit easier to see)
    # cond_block = factor(paste0(cond, "_", block), levels = c("congruent_most_con", "congruent_most_incon", "incongruent_most_con", "incongruent_most_incon"))
    cond_block = paste0(cond, "_", block)
  ) %>% group_by(id, block_number) %>% arrange(trial) %>%
  mutate(
    run_trial = 1:40, #NB. This is a risky solution in general, but I've verified that all subjects show the correct ordering and don't have gaps
    run_trial_z = as.vector(scale(run_trial)),
    prev_rt=dplyr::lag(rt, 1, order_by=run_trial)/1000, #convert to seconds to put it on a smaller scale, roughly comparable to the -1000/rt rt_inv
    prev_rt_inv=dplyr::lag(rt_inv, 1, order_by=run_trial)
  ) %>% arrange(id, trial) %>% ungroup()

#sanity checks on structure
xtabs(~run_trial + block_number, flanker)
xtabs(~block + block_number, flanker)
xtabs(~run_trial + block, flanker)

str(flanker) #after basic processing



# sanity checks -----------------------------------------------------------


# Looks like the number of congruent and incongruent trials by block matches expectation
# xtabs(~cond + block, flanker)
prop.table(xtabs(~cond + block, flanker), 2)


# A few super-short RTs, and then a long RT tail

ggplot(flanker, aes(x=cond, y=rt)) + geom_boxplot() + facet_wrap(~block) + theme_cowplot()

#missed trials
#flanker %>% filter(is.na(rt))



accuracies <- flanker %>% group_by(id, cond) %>%
  summarise(nacc=sum(correct==1, na.rm=T), nerr=sum(correct==0, na.rm=T)) %>%
  pivot_wider(names_from="cond", values_from=c("nerr", "nacc")) %>%
  mutate(
    p_acc=(nacc_congruent + nacc_incongruent)/(nacc_congruent + nacc_incongruent + nerr_congruent + nerr_incongruent),
    p_acc_congruent=(nacc_congruent)/(nacc_congruent + nerr_congruent),
    p_acc_incongruent=(nacc_incongruent)/(nacc_incongruent + nerr_incongruent)
  )

g1 <- ggplot(accuracies, aes(y=p_acc)) + geom_boxplot() + ggtitle("Overall Accuracy")
g2 <- ggplot(accuracies, aes(y=p_acc_incongruent)) + geom_boxplot() + ggtitle("Incongruent Accuracy")
g3 <- ggplot(accuracies, aes(y=p_acc_congruent)) + geom_boxplot() + ggtitle("Congruent Accuracy")
plot_grid(g1, g2, g3, nrow=1)

#exclusion subjects from dissertation
accuracies %>% filter(id %in% c(7, 14, 36, 38, 47, 55, 62, 83))

## Unlikely values

## Look at harmonic means

## Subject exclusions

# Looking over exclusions from dissertation, I might advocate a more conservative approach (keep more data). But there are clearly folks who are systematically wrong on incongruent (< 50%), suggesting that they did not understand the task.

### Exclusion coding

# Generate a column called exclude_flanker with the following scheme:
#   
# - 0: a good subject, no concerns
# - 1: a subject with minor problems that we should try to salvage
# - 2: a suspicious subject whose behavior is questionable, but not obviously bad
# - 3: a subject who should certainly be excluded.
# 
# Following this scheme, here are the combined recommendations of Nate and Michael. Here, I've built on the code above, plus Nate's FlankerPreprocessing R Markdown.

#### 3: obviously bad subjects

# Systematically wrong on incongruent trials. Suggests that they simply didn't understand the task.
# 
# - 7: 11.1% accuracy on incongruent; 97.5% accuracy congruent
# - 47: 18.0% accuracy on incongruent; 93.2% accuracy congruent
# - 55: 6.4% accuracy on incongruent; 98.8% congruent
# - 83: 11.8% accuracy on incongruent; 93.8% congruent

#### 2: probably bad

# At chance, or somewhat below, on incongruent trials
# 
# - 14: 30.4% accuracy on incongruent; 98.8% accuracy congruent
# - 38: 51.2% accuracy on incongruent; 97.5% congruent
# - 62: 47.8% accuracy on incongruent; 90.3% congruent

#### 1: questionable, but not clearly flawed

# - 36: 88.3% accuracy on incongruent; 72.2% accuracy congruent
# - 15: mean RT > 3SDs above sample harmonic mean (per Nate)
# - 111: mean RT > 3SDs above sample harmonic mean (per Nate)
# - 90: strange behavior on incon (per Nate)?
# 
# I’d say we’re justified in dropping these folks given that both 15 and 111 are greater than 3SDs above the mean RT for the crossed stimulus type and block type cells (even after winsorizing). 90 is whacky bc they are so far out on the incong-incong combination. I’d say this warrants us booting them in order to be conservative.


flanker <- flanker %>% mutate(
  exclude_flanker=dplyr::recode(id,
                                `7`  = 3, `47` = 3, `55` = 3, `83` = 3,
                                `14` = 2, `38` = 2, `62` = 2, 
                                `36` = 1, `15` = 1, `111` = 1, `90` = 1,
                                .default = 0
  ))


# RT problems -------------------------------------------------------------

# Following Nate's document, let's treat RTs < 250ms as implausible. I started out with < 200ms as the target, but then in looking at transformed RTs, keeping the
# 200-250ms range led to a slightly long left tail in the distribution after the 1/RT transformation. I solved this by trimming the bottom 1% of RTs by condition and
# block (at group level), but ironically, this specifically targeted RTs between 201ms and 250ms. So, we land on the 250ms threshold either way we slice things, and
# the data are beautifully normal under a 1/RT transformation with this threshold.

#Number of trials with RTs < 250.


nrow(flanker %>% filter(rt < 250))

#This is about 1.3% of the data, which is a very minor loss.

#Set accuracy and RT to NA when RT < 250ms.


flanker <- flanker %>% mutate(rt=if_else(rt < 250, NA_real_, rt), correct=if_else(rt < 250, NA_real_, correct))

# Winsorizing long RTs:


# What about Winsorizing the top 1% of inv-transformed RTs by condition and block?
flanker <- flanker %>% group_by(cond, block) %>% 
  mutate(
    rt_log=log10(rt), 
    rt_log_winsor=DescTools::Winsorize(rt_log, probs=c(0, .99), na.rm = T),
    rt_winsor=DescTools::Winsorize(rt, probs=c(0, .99), na.rm = T) 
  ) %>%
  ungroup() %>% group_by(id, cond, block) %>% #use per-subject condition and block RT quantiles to trim
  mutate(
    rt_inv_trim=if_else(rt_inv < quantile(rt_inv, 0.01, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_trim=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_log_trim=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup() %>%
  group_by(cond, block) %>% #use overall condition and block RT quantiles to trim
  mutate(
    rt_inv_trim_grp=if_else(rt_inv < quantile(rt_inv, 0.01, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_trim_grp=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_log_trim_grp=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup()


#the problem with per-subject trimming is that many plausible values are removed since they are unlikely compared to the subject's distribution
# flanker %>% filter(is.na(rt_inv_trim)) %>% pull(rt) %>% hist(main="RTs dropped by subject-specific trim")
# flanker %>% filter(is.na(rt_inv_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")
# flanker %>% filter(is.na(rt_inv_trim_grp)) %>% nrow()
flanker %>% filter(is.na(rt)) %>% nrow()
flanker %>% filter(is.na(rt_trim_grp)) %>% nrow()
flanker %>% filter(is.na(rt_winsor)) %>% nrow()
flanker %>% filter(is.na(rt_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")

trim <- ggplot(data = flanker, aes(x = rt_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = flanker, aes(x = rt_trim_grp)) + geom_density()
winsor <- ggplot(data = flanker, aes(x = rt_winsor)) + geom_histogram()
winsor_d <- ggplot(data = flanker, aes(x = rt_winsor)) + geom_density()

plot_grid(trim, trim_d, winsor, winsor_d)
314/nrow(flanker)

## in terms of keeping distributions smooth, lets stick with group-level trimming


# RT distributions for clean and full sample  -----------------------------

flanker_SNAP <- flanker %>% left_join(select(SNAP_all_scored, id, exclude_SNAP))

flanker_full <- dplyr::filter(flanker_SNAP, exclude_SNAP ==0, exclude_flanker != 3)
flanker_mean_ci_conditions_full <- summarySEwithin(data = flanker_full, measurevar = "rt", withinvars = "cond_block", idvar = "id", na.rm = TRUE)
full_n  <- length(unique(flanker_full$id))

flanker_point_full <- ggplot(flanker_mean_ci_conditions_full, aes(x = rt, y = cond_block)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt-ci, xmax=rt+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("Flanker task full sample: N=", full_n))


flanker_clean <- dplyr::filter(flanker_SNAP, exclude_SNAP ==0, exclude_flanker ==0)
flanker_mean_ci_conditions_clean <- summarySEwithin(data = flanker_clean, measurevar = "rt", withinvars = "cond_block", idvar = "id", na.rm = TRUE)

clean_n <- length(unique(flanker_clean$id))
flanker_point_clean <- ggplot(flanker_mean_ci_conditions_clean, aes(x = rt, y = cond_block)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt-ci, xmax=rt+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("Flanker task clean sample: N=", clean_n))

plot_grid(flanker_point_clean, flanker_point_full, ncol = 1)

# Data cleaning and setup: Recent Probes -------------------------------------------------


#raw recent probes data, all subjects
load(paste0(basedir, "Data/R_Originals/recent_probes_raw.RData")) #loads in as rp

rp <- rp %>% select(Subject, ProbeDisplay_ACC, ProbeDisplay_RT, Condition,  Trial) %>%
  dplyr::rename(id=Subject, correct=ProbeDisplay_ACC, rt=ProbeDisplay_RT, 
                trial=Trial, stim=Condition) %>%
  mutate(
    # trial = trial -8,  #1:10 are practice trials
    stim = factor(stim, levels=c("Y", "Nfam0", "Nfam1", "Nfam2", "Nfam1RI"), labels=c("positive", "negative_unfamiliar", "negative_familiar", "negative_highly_familiar", "negative_rc")),
    cond = factor(ifelse(stim %in% c("positive", "negative_unfamiliar"), 0, 1), levels = c(0,1), labels = c("no_conflict", "conflict")),
    trial_z = as.vector(scale(trial)),
    rt_inv=-1000/rt #puts Rts on the -4 -- -1 range (makes parameter estimates a bit easier to see)
  ) %>% group_by(id) %>% arrange(id,trial)


id_82s <- rp %>% filter(trial == 82) %>%  dplyr::select(id) # fix trial numbering differences between groups

rp <- rp %>% mutate(trial = ifelse(id %in% id_82s$id, trial -10, trial -8)) %>% mutate(
  prev_rt=dplyr::lag(rt, 1)/1000, #convert to seconds to put it on a smaller scale, roughly comparable to the -1000/rt rt_inv
  prev_rt_inv=dplyr::lag(rt_inv, 1)
) %>% arrange(id, trial) %>% ungroup()

#sanity checks on structure
xtabs(~trial , rp)

str(rp) #after basic processing



# sanity checks -----------------------------------------------------------


# A few super-short RTs, and then a long RT tail

ggplot(rp, aes(x=stim, y=rt)) + geom_boxplot() + theme_cowplot() + theme(axis.text.x = element_text(angle = 45))

#missed trials
miss <- rp %>% filter(is.na(rt)) # 59 total
table(miss$id) %>% sort() %>% hist() # id 36 has 11 missed trials (15% of trials)


accuracies <- rp %>% group_by(id, stim) %>%
  dplyr::summarise(nacc=sum(correct==1, na.rm=T), nerr=sum(correct==0, na.rm=T)) %>%
  pivot_wider(names_from="stim", values_from=c("nerr", "nacc")) %>%
  mutate(
    p_acc=(nacc_positive + nacc_negative_unfamiliar + nacc_negative_familiar + nacc_negative_highly_familiar + nacc_negative_rc)/(nacc_positive + nacc_negative_unfamiliar + nacc_negative_familiar + nacc_negative_highly_familiar + nacc_negative_rc + nerr_positive + nerr_negative_unfamiliar + nerr_negative_familiar + nerr_negative_highly_familiar + nerr_negative_rc),
    p_acc_positive=(nacc_positive)/(nacc_positive + nerr_positive),
    p_acc_negative_unfamiliar=(nacc_negative_unfamiliar)/(nacc_negative_unfamiliar + nerr_negative_unfamiliar),
    p_acc_negative_familiar=(nacc_negative_familiar)/(nacc_negative_familiar + nerr_negative_familiar),
    p_acc_negative_highly_familiar=(nacc_negative_highly_familiar)/(nacc_negative_highly_familiar + nerr_negative_highly_familiar),
    p_acc_negative_rc=(nacc_negative_rc)/(nacc_negative_rc + nerr_negative_rc),
  )

g1 <- ggplot(accuracies, aes(y=p_acc)) + geom_boxplot() + ggtitle("Overall Accuracy")
g2 <- ggplot(accuracies, aes(y=p_acc_positive)) + geom_boxplot() + ggtitle("Positive Accuracy")
g3 <- ggplot(accuracies, aes(y=p_acc_negative_unfamiliar)) + geom_boxplot() + ggtitle("Negative Unfamiliar Accuracy")
g4 <- ggplot(accuracies, aes(y=p_acc_negative_familiar)) + geom_boxplot() + ggtitle("Negative Familiar Accuracy")
g5 <- ggplot(accuracies, aes(y=p_acc_negative_highly_familiar)) + geom_boxplot() + ggtitle("Negative Highly familiar Accuracy")
g6 <- ggplot(accuracies, aes(y=p_acc_negative_rc)) + geom_boxplot() + ggtitle("Negative Response Conflict Accuracy")
plot_grid(g1, g2, g3, g4,g5,g6, nrow=1)

#exclusion subjects from dissertation
accuracies %>% filter(id %in% c(7, 14, 36, 38, 47, 55, 62, 83)) %>% select(id, starts_with("p_acc"))

accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_positive)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_negative_unfamiliar)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_negative_familiar)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_negative_highly_familiar)
accuracies %>% select(id, starts_with("p_acc")) %>% arrange(p_acc_negative_rc)



## Subject exclusions

### Exclusion coding

# Generate a column called exclude_rp with the following scheme:
#   
# - 0: a good subject, no concerns
# - 1: a subject with minor problems that we should try to salvage
# - 2: a suspicious subject whose behavior is questionable, but not obviously bad
# - 3: a subject who should certainly be excluded.
# 
# Following this scheme, here are the combined recommendations of Nate and Michael. Here, I've built on the code above, plus Nate's rpPreprocessing R Markdown.

#### 2: probably bad

# - 65: 52% acc total, 41% positive, 89% negative unfamiliar, 33% negative familiar, 89% negative highly familiar, 44% negative response conflict. 

#### 1: questionable, but not clearly flawed

# - 40: 54% acc total, 61% positive, 56% negative unfamiliar, 33% negative familiar, 56% negative highly familiar, 44% negative response conflict. 
# - 100: 60% acc total, 67% positive, 78% negative unfamiliar, 78% negative familiar, 22% negative highly familiar, 56% negative response conflict. 
# - 36: 65% acc total, 58% positive, 67% negative unfamiliar, 78% negative familiar, 78% negative highly familiar, 67% negative response conflict. 
# - 38: 61% acc total, 50% positive, 67% negative unfamiliar, 67% negative familiar, 78% negative highly familiar, 78% negative response conflict. 


# Honestly, none of these are super crazy concerning on the basis of accuracy, maybe with the exception of 65 who performed quite low on positive trials, though this could certainly just be indicated by a nay-saying bias.
# 36 and 38 are repeat offenders that should likely be excluded from all analyses. Likewise, I've flagged 100 for an invalid SNAP profile so it doesnt seem that we are losing much at all.

rp <- rp %>% mutate(
  exclude_recent_probes=dplyr::recode(id,
                                      `65` = 2, 
                                      `40` = 1, `36` = 1, `100` = 1, `38` = 1,
                                      .default = 0
  ))


# RT problems -------------------------------------------------------------

# In line with documentation in the RecentProbesPreprocessing.Rmd we will exclude trials that are faster than 350 ms and above 3SDs above 
#the group harmonic mean, calculated separately for each condition.

#Number of trials with RTs < 350.

nrow(rp %>% filter(rt < 350))

#75 RTs. This is about 0.9% of the data, which is a very minor loss.

#Set accuracy and RT to NA when RT < 350ms.


rp <- rp %>% mutate(rt=if_else(rt < 350, NA_real_, rt), correct=if_else(rt < 350, NA_real_, correct))

# Winsorizing long RTs:


# What about Winsorizing the top 1% of inv-transformed RTs by condition and block?
test <- rp %>% group_by(stim) %>% 
  mutate(
    rt_log=log10(rt), 
    rt_log_winsor=DescTools::Winsorize(rt_log, probs=c(0, .99), na.rm = T),
    rt_winsor=DescTools::Winsorize(rt, probs=c(0, .99), na.rm = T) 
  ) %>%
  ungroup() %>% group_by(id, stim) %>% #use per-subject condition and block RT quantiles to trim
  mutate(
    rt_inv_trim=if_else(rt_inv < quantile(rt_inv, 0.01, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_trim=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_log_trim=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup() %>%
  group_by(stim) %>% #use overall condition and block RT quantiles to trim
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

## stick with group-trimmed log-RT

### look for high-NA subjects and unusually high mean RTs for potential exclusion

rt_sums <- test %>% group_by(id) %>% dplyr::summarise(mean_rt = harmonic.mean(rt_log_trim_grp), p_na = length(is.na(rt_log_trim_grp)))
hist(rt_sums$mean_rt) # looks fine to me
# View(rt_sums)
# View(test)
# test

nas <- test[which(is.na(test$rt_log_trim_grp)),]
# hist(nas$id)
table(nas$id) %>% sort() #%>% hist()

# good news is that subjects with a high number of RAs have already been flagged above, I dont see the need to flag anyone else as problematic.
recent_probes <- test

# RT distributions for clean and full sample  -----------------------------

recent_probes_SNAP <- recent_probes %>% left_join(select(SNAP_all_scored, id, exclude_SNAP))

recent_probes_full <- dplyr::filter(recent_probes_SNAP, exclude_SNAP ==0, exclude_recent_probes != 3)
recent_probes_mean_ci_conditions_full <- summarySEwithin(data = recent_probes_full, measurevar = "rt", withinvars = "stim", idvar = "id", na.rm = TRUE)
full_n  <- length(unique(recent_probes_full$id))

recent_probes_point_full <- ggplot(recent_probes_mean_ci_conditions_full, aes(x = rt, y = stim)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt-ci, xmax=rt+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("recent_probes task full sample: N=", full_n))


recent_probes_clean <- dplyr::filter(recent_probes_SNAP, exclude_SNAP ==0, exclude_recent_probes ==0)
recent_probes_mean_ci_conditions_clean <- summarySEwithin(data = recent_probes_clean, measurevar = "rt", withinvars = "stim", idvar = "id", na.rm = TRUE)

clean_n <- length(unique(recent_probes_clean$id))
recent_probes_point_clean <- ggplot(recent_probes_mean_ci_conditions_clean, aes(x = rt, y = stim)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt-ci, xmax=rt+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("recent_probes task clean sample: N=", clean_n))

plot_grid(recent_probes_point_clean, recent_probes_point_full, ncol = 1)

# Data cleaning and setup: Go/NoGo -------------------------------------------------

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
    # rt_inv=-1000/rt #puts Rts on the -4 -- -1 range (makes parameter estimates a bit easier to see). This should be estimated after trimming RTs
  ) %>% group_by(id) %>% arrange(id,trial) %>% mutate(
    prev_rt=dplyr::lag(rt, 1, order_by = trial)/1000, #convert to seconds to put it on a smaller scale, roughly comparable to the -1000/rt rt_inv
    # prev_rt_inv=dplyr::lag(rt_inv, 1, order_by = trial)
  ) %>% arrange(id, trial) %>% dplyr::group_by(id, cond) %>%
  dplyr::mutate(block_trial = dplyr::row_number()) %>% #in case we want to allow for the potential for the number of trials since the last nogo to influence RT rather than simply blocking by "goproc" or adding overall trial number (blind to sequencing of go and nogo trials) as covariate
  ungroup() 

for(i in 1:8){gng %>% filter(trial ==i) %>% print()}

#sanity checks on structure
xtabs(~cond +id, gng) # 12 of each condition.
xtabs(~stim +cond , gng)
# xtabs(~Block+id , gng)


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
# accuracies %>% filter(id %in% c(7, 14, 36, 38, 47, 55, 62, 83)) %>% select(id, starts_with("p_acc"))

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
# - 7: This person conversely is quite low in accuracy on all of the no-go trials, however their go performance is close to 1 on go trials making me think they are perhaps just strongly biased towards go.


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


gng <- gng %>% mutate(rt=if_else(rt < 200 & rt != 0, NA_real_, rt), correct=if_else(rt < 200 & rt != 0, NA_real_, correct)) #since 0's are valid no-gos


# winsorize and trim long RTs ---------------------------------------------


test <- gng  %>% 
  mutate(rt_inv=-1000/rt, #puts Rts on the -4 -- -1 range (makes parameter estimates a bit easier to see).
         rt_log=log10(rt)) %>%
  group_by(stim_cond) %>% 
  mutate( #use overall condition and RT quantiles. 
    rt_winsor_grp= DescTools::Winsorize(rt, probs=c(0, .99), na.rm = T), 
    rt_inv_winsor_grp= DescTools::Winsorize(rt_inv, probs=c(0, .99), na.rm = T),
    rt_log_winsor_grp= DescTools::Winsorize(rt_log, probs=c(0, .99), na.rm = T),
    
    rt_trim_grp=if_else(rt > quantile(rt, 0.98, na.rm=T), NA_real_, rt),
    rt_inv_trim_grp=if_else(rt_inv > quantile(rt_inv, 0.98, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_log_trim_grp=if_else(rt_log > quantile(rt_log, 0.98, na.rm=T), NA_real_, rt_log)
  ) %>%
  ungroup() %>% group_by(id, stim_cond) %>% #use per-subject condition and block RT quantiles to trim
  mutate(
    rt_trim=if_else(rt > quantile(rt, 0.99, na.rm=T), NA_real_, rt),
    rt_inv_trim=if_else(rt_inv < quantile(rt_inv, 0.99, na.rm=T), NA_real_, rt_inv), #rt_inv > quantile(rt_inv, 0.99) | 
    rt_log_trim=if_else(rt_log > quantile(rt_log, 0.99, na.rm=T), NA_real_, rt_log)
  ) %>% ungroup() 


#the problem with per-subject trimming is that many plausible values are removed since they are unlikely compared to the subject's distribution
# test %>% filter(is.na(rt_inv_trim)) %>% pull(rt) %>% hist(main="RTs dropped by subject-specific trim")
# test %>% filter(is.na(rt_inv_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")

hist(test[which(is.na(test$rt_trim_grp)),"rt"]$rt, main = "RTs dropped by group trim: raw")
drop_rt <- test[which(is.na(test$rt_trim_grp)),c("rt", "id", "stim_cond")]; xtabs(~id+stim_cond, drop_rt)
test[which(is.na(test$rt_trim_grp)),"rt"]$rt %>% length()


hist(test[which(is.na(test$rt_log_trim_grp)),"rt"]$rt, main = "RTs dropped by group trim: log")
drop_rt <- test[which(is.na(test$rt_log_trim_grp)),c("rt_log", "id", "stim_cond")]; xtabs(~id+stim_cond, drop_rt)
test[which(is.na(test$rt_log_trim_grp)),"rt"]$rt %>% length()

hist(test[which(is.na(test$rt_inv_trim_grp)),"rt"]$rt, main = "RTs dropped by group trim: inv")
drop_rt <- test[which(is.na(test$rt_trim_grp)),c("rt", "id", "stim_cond")]; xtabs(~id+stim_cond, drop_rt)
test[which(is.na(test$rt_trim_grp)),"rt"]$rt %>% length()




hist(test[which(is.na(test$rt_inv_trim_grp)),"rt"]$rt)
test[which(is.na(test$rt_inv_trim_grp)),c("rt", "rt_inv_trim_grp", "stim_cond", "correct")]
test %>% filter(is.na(rt_inv_trim_grp)) %>% nrow()
test %>% filter(is.na(rt)) %>% nrow()
test %>% filter(is.na(rt_trim_grp)) %>% nrow()
test %>% filter(is.na(rt_winsor_grp)) %>% nrow()
test %>% filter(is.na(rt_trim_grp)) %>% pull(rt) %>% hist(main="RTs dropped by group trim")

# raw RT distributions
trim <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_trim_grp)) + geom_density()
winsor <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_winsor_grp)) + geom_histogram()
winsor_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_winsor_grp)) + geom_density()

sub <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_trim)) + geom_histogram()
sub_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_trim)) + geom_density()

raw <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt)) + geom_histogram()
raw_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, winsor, winsor_d, sub, sub_d, ncol = 2)


# Inverse RT distributions
trim <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv_trim_grp)) + geom_density()
# winsor <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_histogram()
# winsor_d <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_density()

sub <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv_trim)) + geom_histogram()
sub_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv_trim)) + geom_density()

raw <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv)) + geom_histogram()
raw_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_inv)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, sub, sub_d, ncol = 2)

# log RT distributions
trim <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log_trim_grp)) + geom_histogram() 
trim_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log_trim_grp)) + geom_density()
# winsor <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_histogram()
# winsor_d <- ggplot(data = test, aes(x = rt_inv_winsor)) + geom_density()

sub <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log_trim)) + geom_histogram()
sub_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log_trim)) + geom_density()

raw <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log)) + geom_histogram()
raw_d <- ggplot(data = filter(test, stim == "Go" & resp == 1), aes(x = rt_log)) + geom_density()

plot_grid(raw, raw_d, trim, trim_d, sub, sub_d, ncol = 2)

g1 <- ggplot(test, aes(sample=rt_log_winsor_grp)) + stat_qq(distribution=qnorm) + stat_qq_line(color="blue") +
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

## go with group-trimmed log transformation

### look for high-NA subjects and unusually high mean RTs for potential exclusion

rt_sums <- test %>% group_by(id) %>% dplyr::summarise(mean_rt = harmonic.mean(rt_log_trim_grp), p_na = sum(is.na(rt_log_trim_grp)))
# test %>% filter(trial %in% 1:8) %>% select(id, stim_cond, resp, correct, trial, rt_log_trim_grp, rt) 
hist(rt_sums$mean_rt) 
hist(rt_sums$p_na) 
rt_sums %>% arrange(-mean_rt)
rt_sums %>% arrange(-p_na)

# yeah, lets get 15 out of here, they have waay more NAs than anyone else, clearly making them an outlier.

gng <- gng %>% mutate(
  exclude_go_nogo=dplyr::recode(id,
                                `15` = 3, 
                                `7` = 1,
                                .default = 0
  ))

# good news is that subjects with a high number of RAs have already been flagged above, I dont see the need to flag anyone else as problematic.
gng <- test %>% select(id, Block, stim_cond, trial, rt_log_trim_grp, rt_trim_grp) %>% right_join(gng, by = c("id", "Block", "stim_cond", "trial")) %>% 
  mutate(rt_log_trim_grp = ifelse(rt == 0,0, rt_log_trim_grp)) #this simply transforms -Inf values for nogo trials back to 0. These are "true" no-gos.


# looking big picture at the "cleanliness of our data"  ------------------------------------------------------

# detach("package:plyr", unload=TRUE, force = TRUE, character.only = TRUE)
# ssubs <- gng %>% dplyr::filter(trial == 1) 
# xtabs(~exclude_go_nogo, ssubs)

gng_SNAP <- gng %>% left_join(select(SNAP_all_scored, id, exclude_SNAP))


# 2/24/20. MPQ cleaning is done, load and include these folks.
load(file.path(basedir,"Data/preprocessed/MPQ_all_scored_final.RData"))

MPQ_all_scored <- dplyr::rename(MPQ_all_scored, id = subject)
flanker_SNAP <- flanker_SNAP %>% left_join(select(MPQ_all_scored, id, exclude_MPQ), by = "id")
recent_probes_SNAP <- recent_probes_SNAP %>% left_join(select(MPQ_all_scored, id, exclude_MPQ), by = "id")
gng_SNAP <- gng_SNAP %>% left_join(select(MPQ_all_scored, id, exclude_MPQ), by = "id")


# combine all exclusion information and see if folks need to be dropped listwise?
exclude_all <- flanker_SNAP %>% filter(trial == 1) %>% select(id, exclude_flanker, exclude_SNAP,exclude_MPQ) %>% full_join(select(filter(recent_probes_SNAP, trial ==1), id, exclude_recent_probes, exclude_SNAP, exclude_MPQ)) %>% full_join(select(filter(gng_SNAP, trial ==1), id, exclude_go_nogo)) %>% arrange(id) %>% print(n = 112)


exclude_all %>% filter(is.na(exclude_MPQ) | exclude_MPQ ==1 )
MPQ_all_scored %>% filter(id %in% c(7,44,48,82,85))

# if SNAP profile is bad, dont include in any estimations -----------------


exclude_all %>% filter(exclude_SNAP != 0 | is.na(exclude_SNAP) | is.na(exclude_MPQ) | exclude_MPQ == 1)
exclude_all %>% filter(!(exclude_SNAP != 0 | is.na(exclude_SNAP) | is.na(exclude_MPQ) | exclude_MPQ == 1)) %>% print(n = 112)

# Invalidity indices for all personality (e.g. listwise) dropped subjects. I'm generally okay with this.
SNAP_all_scored %>% full_join(MPQ_all_scored, by = "id") %>% select(id, contains("VRIN"), contains("TRIN"), Z_InvalidityIndex ) %>% select(!starts_with("T_")) %>% 
  filter(id %in% c(7,23,38,48,51,82,85,100)) %>%
  dplyr::rename(Z_VRIN_SNAP = Z_VRIN.x,
                Z_VRIN_MPQ = Z_VRIN.y,
                Z_TRIN_SNAP = Z_TRIN.x,
                Z_TRIN_MPQ = Z_TRIN.y,) %>% select(-VRIN, -TRIN) %>% arrange(id)

exclude_all <- exclude_all %>% mutate(
  exclude_flanker=dplyr::recode(id,
                                `7` = 3, 
                                `23` = 3, 
                                `38` = 3,
                                `44` = 3, # this is only bc it is missing, do not drop for other tasks.
                                `48` = 3, 
                                `51` = 3,
                                `82` = 3, 
                                `85` = 3, 
                                `100` = 3, 
                                .default = exclude_flanker
  ),
  exclude_recent_probes=dplyr::recode(id,
                                      `7` = 3, 
                                      `23` = 3, 
                                      `38` = 3,
                                      `48` = 3, 
                                      `51` = 3,
                                      `82` = 3, 
                                      `85` = 3, 
                                      `100` = 3, 
                                      .default = exclude_recent_probes
  ),
  
  exclude_go_nogo=dplyr::recode(id,
                                `7` = 3, 
                                `23` = 3, 
                                `38` = 3,
                                `48` = 3, 
                                `51` = 3,
                                `82` = 3, 
                                `85` = 3, 
                                `100` = 3,
                                .default = exclude_go_nogo
  ))

##### terrible form, but apply this to individual dfs

flanker <- flanker_SNAP %>% mutate(
  exclude_flanker=dplyr::recode(id,
                                `7` = 3, 
                                `23` = 3, 
                                `38` = 3,
                                `48` = 3, 
                                `51` = 3,
                                `82` = 3, 
                                `85` = 3, 
                                `100` = 3,
                                .default = exclude_flanker
  ))

recent_probes <- recent_probes_SNAP %>% mutate(
  exclude_recent_probes=dplyr::recode(id,
                                      `7` = 3, 
                                      `23` = 3, 
                                      `38` = 3,
                                      `48` = 3, 
                                      `51` = 3,
                                      `82` = 3, 
                                      `85` = 3, 
                                      `100` = 3,
                                      .default = exclude_recent_probes
  ))
gng <- go_nogo <- gng_SNAP %>% mutate(
  exclude_go_nogo =dplyr::recode(id,
                                 `7` = 3, 
                                 `23` = 3, 
                                 `38` = 3,
                                 `48` = 3, 
                                 `51` = 3,
                                 `82` = 3, 
                                 `85` = 3, 
                                 `100` = 3,
                                 .default = exclude_go_nogo
  ))



# RT distributions for clean and full sample: GNG  -----------------------------



gng_full <- dplyr::filter(gng, exclude_SNAP ==0, !is.na(exclude_SNAP), exclude_go_nogo != 3, stim == "Go")
gng_mean_ci_conditions_full <- summarySEwithin(data = gng_full, measurevar = "rt", withinvars = "cond", idvar = "id", na.rm = TRUE)
full_n  <- length(unique(gng_full$id))

gng_point_full <- ggplot(gng_mean_ci_conditions_full, aes(x = rt_norm, y = cond)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt_norm-ci, xmax=rt_norm+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("gng task full sample: N=", full_n))


gng_clean <- dplyr::filter(gng, exclude_SNAP ==0, exclude_go_nogo ==0)
gng_mean_ci_conditions_clean <- summarySEwithin(data = gng_clean, measurevar = "rt", withinvars = "cond", idvar = "id", na.rm = TRUE)

clean_n <- length(unique(gng_clean$id))
gng_point_clean <- ggplot(gng_mean_ci_conditions_clean, aes(x = rt_norm, y = cond)) +
  geom_point(shape = 21, size = 2, fill = "black") +
  geom_errorbar(width = 0, aes(xmin = rt_norm-ci, xmax=rt_norm+ci)) + 
  labs(y = "Condition", x = "Reaction Time (msec)",title = paste0("gng task clean sample: N=", clean_n))

plot_grid(gng_point_clean, gng_point_full, ncol = 1)
# not entirely clear what's going on at the bottom here, but fine not getting too bogged down for now

























# final look and rename vars if need be -----------------------------------

##### flanker

#only drop the folks that are definitely not usable.
flanker_full <- flanker %>% filter(exclude_flanker != 3) %>% select(id, correct, rt_trim_grp, cond, block, block_number,trial_z, rt_inv_trim_grp, block_number, run_trial, run_trial_z, prev_rt, prev_rt_inv) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp, stim = cond) %>% 
  mutate(stimblock = ifelse(as.character(stim) == "incongruent", paste0(stim, "_", block), as.character(stim))) # descriptives indicate that perhaps block only matters for incongruent stimuli, thus, create a contrast that lumps congruent together but separates stimulus by block for incongruent

write.csv(flanker_full, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/flanker_full_accCode.csv", row.names = FALSE)

# filter NA RTs as well
flanker_full_nafilt <- flanker %>% filter(exclude_flanker != 3) %>% select(id, correct, rt_trim_grp, cond, block, block_number,trial_z, rt_inv_trim_grp, block_number, run_trial, run_trial_z, prev_rt, prev_rt_inv) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp, stim = cond) %>% 
  mutate(stimblock = ifelse(as.character(stim) == "incongruent", paste0(stim, "_", block), as.character(stim))) %>% #descriptives indicate that perhaps block only matters for incongruent stimuli, thus, create a contrast that lumps congruent together but separates stimulus by block for incongruent
  filter(!is.na(rt)) # N.B. these need to be filtered out after changing the name of group-trimmed rt to rt for analysis.

write.csv(flanker_full_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/flanker_full_nafilt_accCode.csv", row.names = FALSE)

#drop the questionable folks too
flanker_clean <- flanker %>% filter(exclude_flanker == 0) %>% select(id, correct, rt_trim_grp, cond, block, block_number,trial_z, rt_inv_trim_grp, block_number, run_trial, run_trial_z, prev_rt, prev_rt_inv) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp, stim = cond) %>% 
  mutate(stimblock = ifelse(as.character(stim) == "incongruent", paste0(stim, "_", block), as.character(stim))) #descriptives indicate that perhaps block only matters for incongruent stimuli, thus, create a contrast that lumps congruent together but separates stimulus by block for incongruent
  
write.csv(flanker_clean, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/flanker_clean_accCode.csv", row.names = FALSE)

# filter NA RTs as well
flanker_clean_nafilt <- flanker %>% filter(exclude_flanker == 0) %>% select(id, correct, rt_trim_grp, cond, block, block_number,trial_z, rt_inv_trim_grp, block_number, run_trial, run_trial_z, prev_rt, prev_rt_inv) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp, stim = cond) %>%
  mutate(stimblock = ifelse(as.character(stim) == "incongruent", paste0(stim, "_", block), as.character(stim))) %>% #descriptives indicate that perhaps block only matters for incongruent stimuli, thus, create a contrast that lumps congruent together but separates stimulus by block for incongruent
  filter(!is.na(rt)) # N.B. these need to be filtered out after changing the name of group-trimmed rt to rt for analysis.
  
write.csv(flanker_clean_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/flanker_clean_nafilt_accCode.csv", row.names = FALSE)




##### recent_probes

#only drop the folks that are definitely not usable.
recent_probes_full <- recent_probes %>% filter(exclude_recent_probes != 3) %>% select(id, correct, rt_trim_grp, stim, cond, trial_z, rt_log_trim_grp, prev_rt) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) 

write.csv(recent_probes_full, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/recent_probes_full_accCode.csv", row.names = FALSE)

# filter NA RTs as well
recent_probes_full_nafilt <- recent_probes %>% filter(exclude_recent_probes != 3)%>% select(id, correct, rt_trim_grp, stim, cond, trial_z, rt_log_trim_grp, prev_rt) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) %>% 
  filter(!is.na(rt))
  
write.csv(recent_probes_full_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/recent_probes_full_nafilt_accCode.csv", row.names = FALSE)

#drop the questionable folks too
recent_probes_clean <- recent_probes %>% filter(exclude_recent_probes == 0) %>% select(id, correct, rt_trim_grp, stim, cond, trial_z, rt_log_trim_grp, prev_rt) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp)

write.csv(recent_probes_clean, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/recent_probes_clean_accCode.csv", row.names = FALSE)

# filter NA RTs as well
recent_probes_clean_nafilt <- recent_probes %>% filter(exclude_recent_probes == 0) %>% select(id, correct, rt_trim_grp, stim, cond, trial_z, rt_log_trim_grp, prev_rt) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) %>% filter(!is.na(rt))

write.csv(recent_probes_clean_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/recent_probes_clean_nafilt_accCode.csv", row.names = FALSE)

##### go_nogo

#only drop the folks that are definitely not usable.
go_nogo_full <- go_nogo %>% filter(exclude_go_nogo != 3) %>% select(id, correct, rt_trim_grp, cond, stim, cond, stim_cond, trial_z, rt_log_trim_grp, prev_rt, block_trial) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) 

write.csv(go_nogo_full, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/go_nogo_full_accCode.csv", row.names = FALSE)

# filter NA RTs as well
go_nogo_full_nafilt <- go_nogo %>% filter(exclude_go_nogo != 3) %>% select(id, correct, rt_trim_grp, cond, stim, cond, stim_cond, trial_z, rt_log_trim_grp, prev_rt, block_trial) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) %>%
  filter(!is.na(rt))

write.csv(go_nogo_full_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/go_nogo_full_nafilt_accCode.csv", row.names = FALSE)

#drop the questionable folks too
go_nogo_clean <- go_nogo %>% filter(exclude_go_nogo == 0) %>% select(id, correct, rt_trim_grp, cond, stim, cond, stim_cond, trial_z, rt_log_trim_grp, prev_rt, block_trial) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) 

write.csv(go_nogo_clean, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/go_nogo_clean_accCode.csv", row.names = FALSE)

# filter NA RTs as well
go_nogo_clean_nafilt <- go_nogo %>% filter(exclude_go_nogo == 0) %>% select(id, correct, rt_trim_grp, cond, stim, cond, stim_cond, trial_z, rt_log_trim_grp, prev_rt, block_trial) %>%
  dplyr::rename(subj_idx = id, response = correct, rt = rt_trim_grp) %>%
  filter(!is.na(rt))

write.csv(go_nogo_clean_nafilt, file = "~/github_repos/PD_Inhibition_DDM/Data/preprocessed/go_nogo_clean_nafilt_accCode.csv", row.names = FALSE)

