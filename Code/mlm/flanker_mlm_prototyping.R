#MLMs on Flanker Data
library(readr)
library(nlme)
library(brms)
library(tidyverse)
library(emmeans)

pdddm_home <- "~/Data_Analysis/PD_Inhibition_DDM"

flanker <- read_csv(file.path(pdddm_home, "Data/preprocessed", "flank_stimCode.csv")) %>%
  dplyr::rename(acc=TrialSlide_ACC, cond=stim, block=CongruentBlock) %>% select(-response) %>%
  mutate(rt=rt*1000, rt_sec=rt/1000, 
         cond=factor(cond, levels=c(0,1), labels=c("congruent", "incongruent")),
         block=factor(block, levels=c(0,1), labels=c("most_incon", "most_con"))) %>%
  group_by(subj_idx) %>%
  mutate(trial=1:n(), trial_z = as.vector(scale(trial)), prev_cond=dplyr::lag(cond, 1, order_by=trial)) %>% ungroup() %>%
  mutate(prev_cond=ifelse(trial_block==1, NA_character_, as.character(prev_cond))) %>% #reset prev_cond by trial
  na.omit() #drop the first trial of every block so that models can be compared directly


# table(flanker$TrialSlide_ACC)
# table(flanker$block)
# table(flanker$response)
# table(flanker$stim)
xtabs(~cond + block, flanker)
xtabs(~trial, flanker)

#starting point in ML: random intercept of subject, CS structure wrt trial
m1 <- lme(rt ~ trial + cond*block, random = ~ 1 | subj_idx,
          correlation=corCompSymm(form=~trial|subj_idx), na.action = na.exclude, 
          data = flanker, method='ML')

summary(m1)
emmeans(m1, ~cond | block)

#compound symmetry is already captured above by a simple variance component (G-side)
#which is why the CS correlation is zero in a model that has both
#go back to a simple random intercept model
m2 <- lme(rt ~ trial + cond*block, random = ~ 1 | subj_idx,
          na.action = na.exclude, 
          data = flanker, method='ML', control = lmeControl(opt="optim"))

summary(m2)
anova(m1, m2)
emmeans(m2, ~cond | block)

#same log likelihood!
logLik(m1) - logLik(m2)

#blocked compound symmetry approach
#random=list(Block=pdCompSymm(~Variety-1))

#Here's the R-side only model with compound symmetry (no random intercept)
m3 <- gls(rt ~ trial + cond*block, correlation=corCompSymm(form=~trial|subj_idx), na.action = na.exclude, 
          data = flanker, method='ML')

#this should also yield identical fit statistics -- just a reparameterization of the same approach
anova(m2, m3)
logLik(m2) - logLik(m3)

#look at model-predicted residual correlation structure

#R.a <- getVarCov(dental.lme.a,type="conditional",individual=1)
rcov <- getVarCov(m3, type="conditional", individual=1) #covariance matrix
rcorr <- cov2cor(getVarCov(m3, type="conditional",individual=1)) #residual correlation matrix for subj 1

#V_i <- getVarCov(m3, type="marginal", individual=1)

#how about random intercept + random slope of trial
m4 <- lme(rt ~ trial + cond*block, random = ~ 1 +trial | subj_idx,
          na.action = na.exclude, 
          data = flanker, method='ML')

summary(m4)

#AR(1) model on the R-side with random intercept and slope on the G-side
m5 <- lme(rt ~ trial + cond*block, 
          random = ~ trial | subj_idx,
          correlation=corAR1(form = ~ trial_z | subj_idx),
          na.action = na.exclude, 
          data = flanker, method='ML')

anova(m4, m5) #suggests that an AR(1) structure on residuals helps -- ~0.3 correlation for adjacent trials
rcov <- getVarCov(m5, type="conditional", individual=1) #covariance matrix
rcorr <- cov2cor(getVarCov(m5, type="conditional",individual=1)[[1]]) #residual correlation matrix for subj 1


###          


m4 <- lme(rt ~ trial_z + cond*block, random = list(subj_idx=pdCompSymm(~block - 1)),
          na.action = na.exclude, #correlation=corCompSymm(form=~trial_z|subj_idx),
          data = flanker, method='REML')

summary(m3)

#I believe this is a model where intercepts are allowed to vary by condition
m5 <- lme(rt ~ trial_z + cond*block, random = list(subj_idx=pdDiag(form=~block)),
          na.action = na.exclude, #correlation=corCompSymm(form=~trial_z|subj_idx),
          data = flanker, method='REML')

VarCorr(m4)
getVarCov(m4, "random.effects")

#https://m-clark.github.io/mixed-models-with-R/extensions.html
#heterogeneous intercept by block
m6 <- lme(rt ~ trial_z + cond*block, weights=varIdent(form=~1|block),
          random=~1|subj_idx,
          na.action = na.exclude,
          data = flanker, method='ML')

getVarCov(m6, "random.effects")

library(glmmTMB)


m7 <- lme(rt ~ trial_z + cond*block, random=~1|subj_idx/block,
          na.action = na.exclude,
          data = flanker, method='ML')

summary(m7)
getVarCov(m7, "random.effects")

anova(m6, m7)

library(glmmTMB)
test <- glmmTMB(rt ~ trial_z + cond*block + cs(trial_z + 0 | subj_idx),
          na.action = na.exclude, dispformula = ~0,
          data = flanker, REML = TRUE)


test <- glmmTMB(rt ~ trial_z + cond*block + (1|subj_idx) + diag(0 + block|subj_idx),
                na.action = na.exclude,
                data = flanker, REML = FALSE)


AICtab(m6, m7, test)
AIC(m7)
AIC(test)
library(bbmle)

#what about allowing for a random intercept for condition within block?
m8 <- lme(rt ~ trial_z + cond*block, random=~1|subj_idx/block/cond,
                na.action = na.exclude,
                data = flanker, method='ML')

anova(m7, m8) #no improvement

#what about effect of trial condition as nesting (no block)
m9 <- lme(rt ~ trial_z + cond*block, random=~1|subj_idx/cond,
          na.action = na.exclude,
          data = flanker, method='ML')

anova(m7, m9) #m9 is way worse

#add random slope of linear trial_z to m7
m10 <- lme(rt ~ trial_z + cond*block, random=~1+trial_z|subj_idx/block,
          na.action = na.exclude,
          data = flanker, method='ML')

anova(m7, m10) #looks like a bit improvement to let trial_z slope vary

#add ar1 residual correlation structure on top of block within subject and linear trial_z m10
#N.B. The form must use trial, not trial_z for lme to see the trial structure.
#  If we use trial_z in the corAR1(), lme sees this as a continuous predictor.
#  It's fine to use trial_z on the fixed side of the model, and even to have a random slope
#  of trial_z, which captures individual differences in the linear drift in RTs over the experiment
m11 <- lme(rt ~ trial_z + cond*block, 
          random = ~ 1 + trial_z | subj_idx/block,
          correlation=corAR1(form = ~ trial | subj_idx/block),
          na.action = na.exclude, 
          data = flanker, method='ML')

anova(m10, m11) #huge improvement -- adjacent trials correlate

#ARMA(1,1) even better
m12 <- lme(rt ~ trial_z + cond*block, 
           random = ~ 1 + trial_z | subj_idx/block,
           correlation=corARMA(form = ~ trial | subj_idx/block, p=1, q=1),
           na.action = na.exclude, 
           data = flanker, method='ML')

summary(m12)
anova(m11, m12) #yes, huge improvement here, too

emmeans(m12, ~cond | block)
emmip(m12, ~cond | block, CIs = TRUE)
pairs(emmeans(m12, ~cond * block))


#add previous congruency to model as fixed, in addition to ar1 residual correlation
m13 <- lme(rt ~ trial_z + cond*block*prev_cond, 
           random = ~ 1 + trial_z | subj_idx/block,
           correlation=corAR1(form = ~ trial | subj_idx/block),
           na.action = na.exclude, 
           data = flanker, method='ML')

#simpler version of m7 with prev_cond included, but no R-side error correlation
m14 <- lme(rt ~ trial_z + cond*block*prev_cond, 
           random = ~ 1 + trial_z | subj_idx/block, #random intercept and slope for block nested in subject
           na.action = na.exclude, 
           data = flanker, method='ML')

#ARMA(1,1) and prev_cond
m15 <- lme(rt ~ trial_z + cond*block*prev_cond, 
           random = ~ 1 + trial_z | subj_idx/block,
           correlation=corARMA(form = ~ trial | subj_idx/block, p=1, q=1),
           na.action = na.exclude, 
           data = flanker, method='ML')

m15_brms <- brm(rt_sec ~ 1 + trial_z + cond*block*prev_cond + (1 + trial_z | subj_idx/block) +
                  data=flanker, chains=4, cores=4, iter=3000)
m15_brms <- add_criterion(m15_brms, "waic")
m15_brms <- add_criterion(m15_brms, "loo")

m15_brms2 <- brm(rt_sec ~ 1 + trial_z + cond*block*prev_cond + (1 + trial_z | subj_idx) +
                  (1 | subj_idx:block), data=flanker, chains=4, cores=4, iter=3000)
m15_brms2 <- add_criterion(m15_brms2, "waic")
m15_brms2 <- add_criterion(m15_brms2, "loo")

model_weights(m15_brms, m15_brms2, weights = "loo")
model_weights(m15_brms, m15_brms2, weights = "waic")
loo_model_weights(m15_brms, m15_brms2, weights = "waic")


save(m15_brms, m15_brms2, file="/Users/mnh5174/Data_Analysis/PD_Inhibition_DDM/Outputs/brms/m15_brms.RData")

cor_ar(formula = ~ trial | subj_idx/block, p=1)


brm_2 <- make_stancode(rt ~ 1 + trial + cond*block + cosy(time = trial, gr=subj_idx) + (1 | subj_idx),
                       iter=400, data=flanker, cores = 4)

#AR(2) and prev_cond
m16 <- lme(rt ~ trial_z + cond*block*prev_cond, 
           random = ~ 1 + trial_z | subj_idx/block,
           correlation=corARMA(form = ~ trial | subj_idx/block, p=2, q=0),
           na.action = na.exclude, 
           data = flanker, method='ML')


#anova(m10, m11, m14, m13)
AICtab(m10, m11, m12, m13, m14, m15, m16, sort=TRUE, weights=TRUE)

#so far, m15 is the winner, though it is slow to estimate!

#simpler prototype for a trait
m15 <- lme(rt ~ trial_z + cond*block*prev_cond + cond*block*Z_Impulsivity, 
           random = ~ 1 + trial_z | subj_idx/block,
           #correlation=corARMA(form = ~ trial | subj_idx/block, p=1, q=1),
           na.action = na.exclude, 
           data = flanker, method='ML')

m15 <- lme(rt ~ trial_z + cond*block + cond*block*Z_Entitlement, 
           random = ~ 1 + trial_z | subj_idx/block,
           correlation=corARMA(form = ~ trial | subj_idx/block, p=1, q=0),
           na.action = na.exclude, 
           data = flanker, method='ML')

library(lme4)
mxx <- glmer(acc ~ trial_z + cond*block + (1 | subj_idx),
           #correlation=corARMA(form = ~ trial | subj_idx/block, p=1, q=0),
           na.action = na.exclude, 
           data = flanker, family="binomial")


####
#I believe this is a model where intercepts are allowed to vary by condition
m5 <- gls(rt ~ trial_z + cond*block, weights=varIdent(form=~1|trial_z),
          na.action = na.exclude, correlation=corCompSymm(form=~1|subj_idx),
          data = flanker, method='REML')

xx <- groupedData(rt ~ as.numeric(block) * as.numeric(trial_z) | subj_idx, data=flanker)
m6 <- gls(rt ~ trial_z + cond*block, weights=varIdent(form=~1|block),
          na.action = na.exclude, correlation=corCompSymm(form=~1|subj_idx),
          data = xx, method='REML')



brm_1_mod <- brm(rt ~ 1 + trial_z + cond*block + (1 | subj_idx),
             iter=1200, data=flanker, cores = 4,
             prior=set_prior("normal(200, 200)", class="Intercept")) #improves model tremendously


brm_1_ln <- brm(rt ~ 1 + trial_z + cond*block + (1 | subj_idx),
                iter=1200, data=flanker, cores = 4, family=shifted_lognormal(),
                prior=set_prior("normal(5, 1)", class="Intercept"), #exp(4-6) means intercept is in the 55-403ms range
                file="log_normal_basic_prior")
                 
brm_1_sec <- brm(rt_sec ~ 1 + trial_z + cond*block + (1 | subj_idx),
                iter=1200, data=flanker, cores = 4, family=shifted_lognormal(),
                prior=set_prior("normal(-1.5, 1)", class="Intercept")) #exp(-1.5) yields ~0.22 (s) intercept
                

brm_2_sec <- brm(rt_sec ~ 1 + trial_z + cond*block + ar(time = trial_z, gr=subj_idx, cov=TRUE),
                 iter=1200, data=flanker, cores = 4, family=shifted_lognormal(),
                 prior=set_prior("normal(-1.5, 1)", class="Intercept") #exp(-1.5) yields ~0.22 (s) intercept
                 )



#equivalent brm with compound symmetry
brm_2 <- brm(rt ~ 1 + trial_z + cond*block + cosy(time = trial, gr=subj_idx) + (1 | subj_idx),
             iter=400, data=flanker, cores = 4)


brm_2 <- make_stancode(rt ~ 1 + trial + cond*block + cosy(time = trial, gr=subj_idx) + (1 | subj_idx),
             iter=400, data=flanker, cores = 4)

save(brm_1, brm_2, file="Flanker_brms.RData")

