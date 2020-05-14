
# 4/30/2020 NH edits to latent structure ----------------------------------

## this is mostly me running through the factor analysis checks with faulty subjects dropped. 

pacman::p_load(tidyverse, psych, ggcorrplot, lavaan, GPArotation)

# basedir <- "C:/Users/timot/Documents/GitHub/PD_Inhibition_DDM" # Tim
basedir <- "~/github_dirs/PD_Inhibition_DDM" # Nate

# alldat <- read_csv("/Data/alldat.csv")
mpq <- get(load(file.path(basedir, "Data/preprocessed/MPQ_all_scored_final.RData"))) %>% dplyr::filter(exclude_MPQ == 0)
snap <- get(load(file.path(basedir, "Data/preprocessed/SNAP_all_scored_final.RData"))) %>% dplyr::filter(exclude_SNAP == 0)
alldat <- mpq %>% inner_join(snap, by = "subject") %>% as_tibble()
# length(unique(alldat$subject)) # retain 104 "clean" subjects


#select just the SNAP and MPQ variables
#corrplot of the scales
self_reps <- alldat %>% dplyr::select(subject, starts_with("T_"), MPS_wbT:MPS_abT) %>% select(-T_VRIN, -T_TRIN, -T_DRIN, -T_RareVirtues, - T_Deviance, -T_BackDeviance, -T_InvalidityIndex, -T_NegativeTemperament, -T_PositiveTemperament)
# self_reps <- alldat %>% select(subject, NEGTEMP:SUICPRON, MPS_wbR:MPS_abR)
# self_reps <- self_reps %>% select(-NEGTEMP, -POSTEMP, -DISINH, -DISINHP, -SELFHARM) # I think this simply drops all composite scales
# self_reps <- self_reps %>% select(-NEGTEMP, -POSTEMP, -DISINH, -DISINHP, -SUICPRON, -LOSLFEST) # LOSLFEST tends to be really wonky, perhaps just taking the self harm index would elminate it's wonk-levels. nah I think we should probably drop it altogether.
write.csv(self_reps, file.path(basedir, "Data/preprocessed/all_personality_no_invalids.csv"), row.names = FALSE)
self_reps <- read.csv(file.path(basedir, "Data/preprocessed/all_personality_no_invalids.csv"))
# srcorrs <- cor(select(self_reps, -subject), use = "pairwise.complete.obs")
# ggcorrplot(srcorrs, type = "upper", hc.order = TRUE, tl.cex = 6, lab = TRUE, lab_size = 1)
# ggsave(file.path(basedir,"Outputs/self_report_correlations.pdf"))
# 
# #histograms
# sr_hist <- pivot_longer(self_reps, T_Mistrust:MPS_abT, names_to = "measure", values_to = "value")
# ggplot(sr_hist, aes(x = value)) + geom_histogram(binwidth = 10) + facet_wrap(~measure, scales = "free")
# ggsave(file.path(basedir,"/Outputs/self_report_histograms.pdf"))

# not interested (psychoticism) or messy scale
# sr_fa <- self_reps %>% select(-subject) %>% select(-T_EccentricPerceptions, -MPS_abT, -T_LowSelfEsteem, -T_SuicideProneness, -T_Disinhibition) #-T_SelfHarm)
sr_fa <- self_reps %>% select(-subject) %>% select(-T_EccentricPerceptions, -MPS_abT, -T_SelfHarm, -T_Disinhibition) #-T_SelfHarm)
#map test, parallel analysis, and factor analysis
map <- nfactors(sr_fa, n=10, rotate="oblimin", fm="ml")
map
par <- fa.parallel(sr_fa, fm = "ml", fa = "fa") # suggests 4 factors
sr.fa <- fa(sr_fa, 9, rotate="oblimin", fm = "ml") #suggests 9 factors are suffciient
print(sr.fa, cut=0, sort = T) ##looks interpretable to me; heywood case detected so look at cautiously 
#pa estimation
sr.fa <- fa(sr_fa, 9, rotate="oblimin", fm = "pa") 
print(sr.fa, cut=.3, sort = T) 
print(sr.fa, cut=0, sort = T) 


#8f Still considered sufficient
# sr.fa <- fa(sr_fa, 8, rotate="oblimin", fm = "ml") #ml
# print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 
sr.fa <- fa(sr_fa, 8, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 

#7f Still sufficient
sr.fa <- fa(sr_fa, 7, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 7, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T)

#6f
sr.fa <- fa(sr_fa, 6, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 6, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 
# print(sr.fa, cut=0, sort = T) 

#5f - at this point antagonism and mistrust get collapsed together ; still sufficicent
sr.fa <- fa(sr_fa, 5, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 5, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 

#4f; sufficient
sr.fa <- fa(sr_fa, 4, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 4, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.2, sort = T) 

#3f
sr.fa <- fa(sr_fa, 3, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 3, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 


#2f
sr.fa <- fa(sr_fa, 2, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(sr_fa, 2, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.3, sort = T) 


#what about a 3 factor scale with no neuroticism items
sr.fa <- fa(select(sr_fa, -SUICPRON, -LOSLFEST, -MPS_abR, -ECCPERC, -MPS_srR, -DEPEN, -MPS_haR), 3, rotate="oblimin", fm = "ml") #ml
sr.fa <- fa(select(sr_fa, -SUICPRON, -LOSLFEST, -MPS_abR, -ECCPERC, -MPS_srR, -DEPEN, -MPS_haR), 3, rotate="oblimin", fm = "pa") #pa
print(sr.fa, cut=.2, sort = T) 

#c scales -> MPS_clR not included in self-reports
cscales <- alldat %>% select(MPS_clR, MPS_acR, MPS_tdR, HDWK, PROPER, IMPUL)
c_cor <- cor(cscales, use = "pairwise.complete.obs")
par <- fa.parallel(cscales, fm = "ml", fa = "fa")
c.fa <- fa(cscales, 3, rotate="oblimin", fm = "ml")
print(c.fa, cut=.2, sort = T)

#e scales
escales <- alldat %>% select(EXHIB, ENTITL, DETACH, MPS_wbR, MPS_spR, MPS_scR)
e_cor <- cor(escales, use = "pairwise.complete.obs")
par <- fa.parallel(escales, fm = "ml", fa = "fa")
e.fa <- fa(escales, 2, rotate="oblimin", fm = "ml")
print(e.fa, cut=0.2, sort = T)

#n scales
nscales <- alldat %>% select(LOSLFEST, MPS_srR, SUICPRON, DEPEN, MPS_wbR)
n_cor <- cor(nscales, use = "pairwise.complete.obs")
par <- fa.parallel(nscales, fm = "ml", fa = "fa")
n.fa <- fa(nscales, 3, rotate="oblimin", fm = "ml")
print(n.fa, cut=0.2, sort = T)

#flanker variables
flankvars <- alldat %>% select(C_Error, I_Error, C_MeanRT, I_MeanRT, I_Error_Resid, I_MeanRT_Resid, MISTRUST:SUICPRON, -SELFHARM, -POSTEMP, -DISINHP, -DISINH, MPS_wbR:MPS_abR)
flankcor <- cor(flankvars, use = "pairwise.complete.obs")
ggcorrplot(flankcor, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/flanker_sr_correlations.pdf")

#gng
sr_gng <- alldat %>% select(MeanGoRT, sumMisses, sumFalseAlarms, faresid, MISTRUST:SUICPRON, -SELFHARM, -POSTEMP, -DISINHP, -DISINH, MPS_wbR:MPS_abR)
sr_gng <- cor(sr_gng, use = "pairwise.complete.obs")
ggcorrplot(sr_gng, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/gng_sr_correlations.pdf")

#recent probes
sr_rp <- alldat %>% select(RP_Familiar_ACC, RP_Familiar_RTResid, MISTRUST:SUICPRON, -SELFHARM, -POSTEMP, -DISINHP, -DISINH, MPS_wbR:MPS_abR)
sr_rp <- cor(sr_rp, use = "pairwise.complete.obs")
ggcorrplot(sr_rp, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/recentprobe_sr_correlations.pdf")


### saving factor scores from the 9-f model and looking at task correlations
sr.fa <- fa(sr_fa, 9, rotate="oblimin", fm = "ml") #corrplot looks like maybe we could get 9 factors? 
print(sr.fa, cut=0, sort = T) ##looks interpretable to me
#save factor scores and link to ID
pd_9fscores <- as.data.frame(sr.fa$scores)
pd_9fscores$Subject <- self_reps$Subject
#name the factors
pdfac_cols <- c("Detachment", "AGG", "Agentic_E", "Impulsivity", "Perfection", "Psychoticism", "Neuroticism", "Tradition", "Mistrust", "Subject")
colnames(pd_9fscores) <- pdfac_cols

#pull performance variables from each task
flank_vars <- alldat %>% select(Subject, C_MeanRT, I_MeanRT, C_Error, I_Error, I_Error_Resid, I_MeanRT_Resid)
gng_vars <- alldat %>% select(Subject, MeanGoRT, sumFalseAlarms, sumMisses, faresid)
rp_vars <- alldat %>% select(Subject, RP_Familiar_ACC, RP_Familiar_RTResid)

#merge factors with performance variables from each task
flank_factor_dat <- left_join(pd_9fscores, flank_vars, by = "Subject")
gng_factor_dat <- left_join(pd_9fscores, gng_vars, by = "Subject")
rp_factor_dat <- left_join(pd_9fscores, rp_vars, by = "Subject")

#corrplots for each task with the 9 factors
#flank
flank_fac <- cor(flank_factor_dat, use = "pairwise.complete.obs")
ggcorrplot(flank_fac, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/flanker_9factor_correlations.pdf")
#gng
gng_fac <- cor(gng_factor_dat, use = "pairwise.complete.obs")
ggcorrplot(gng_fac, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/gng_9factor_correlations.pdf")
#rp
rp_fac <- cor(rp_factor_dat, use = "pairwise.complete.obs")
ggcorrplot(rp_fac, type = "upper", tl.cex = 6, lab = TRUE, lab_size = 1)
ggsave("C:/Users/timot/Google Drive/Hallquist_Dombrovski/PD_DDM/Outputs/recentprobe_9factor_correlations.pdf")


# trying some CFAs
B5_model <- '
 O =~ ECCPERC + MPS_abR
 N =~ LOSLFEST + MPS_srR + SUICPRON + DEPEN + MPS_wbR
 A =~ MISTRUST + MPS_alR + MANIP + AGG + MPS_agR + MPS_haR
 C =~ MPS_clR + IMPUL + HDWK + MPS_acR + PROPER + MPS_tdR
 E =~ DETACH + MPS_scR + EXHIB + ENTITL + MPS_spR
 '
B5_fit <- cfa(B5_model, data = self_reps, missing = "ML", estimator = "MLR")
summary(B5_fit, fit.measures = TRUE, standardized = TRUE)
B5mi <- modindices(B5_fit) #modindices mostly want facet-specific error covariances

#add the residual covariances
B5_model_resid <- '
 O =~ ECCPERC + MPS_abR
 N =~ LOSLFEST + MPS_srR + SUICPRON + DEPEN + MPS_wbR
 A =~ MISTRUST + MPS_alR + MANIP + AGG + MPS_agR + MPS_haR
 C =~ MPS_clR + IMPUL + HDWK + MPS_acR + PROPER + MPS_tdR
 E =~ DETACH + MPS_scR + EXHIB + ENTITL + MPS_spR
 
 DETACH	~~	MPS_scR
 HDWK	~~	MPS_acR
 MISTRUST	~~	MPS_alR
 PROPER	~~	MPS_tdR
 MPS_clR	~~	IMPUL
 '
B5_fitr <- cfa(B5_model_resid, data = self_reps, missing = "ML", estimator = "MLR")
summary(B5_fitr, fit.measures = TRUE, standardized = TRUE)
B5mir <- modindices(B5_fitr) #fit not very good here, also doesn't like 2 indicators for O

#what about just going with all 9
B9_model <- '
 O =~ ECCPERC + MPS_abR
 N =~ LOSLFEST + MPS_srR + SUICPRON + DEPEN + MPS_wbR
 TRUST =~ MISTRUST + MPS_alR 
 IMPOLITE =~ MANIP + AGG + MPS_agR + MPS_haR
 CONTROL =~ MPS_clR + IMPUL 
 ACHIEVE =~ HDWK + MPS_acR 
 RULES =~ PROPER + MPS_tdR
 AFF =~ DETACH + MPS_scR 
 AGENT =~ EXHIB + ENTITL + MPS_spR
 '
B9_fit <- cfa(B9_model, data = self_reps, missing = "ML", estimator = "MLR")
summary(B9_fit, fit.measures = TRUE, standardized = TRUE)
B9mi <- modindices(B9_fit) #not gonna work with 2 indicator latents. 

# just a four factor model? 
B4_model <- '
 #O =~ ECCPERC + MPS_abR
 N =~ LOSLFEST + MPS_srR + SUICPRON + MPS_wbR #+ DEPEN
 A =~ MISTRUST + MANIP + MPS_agR #+ AGG  + MPS_haR + MPS_alR 
 C =~ IMPUL + HDWK + PROPER #+ MPS_tdR + MPS_acR + MPS_clR
 E =~ EXHIB + MPS_spR + MPS_scR # + ENTITL + DETACH
#EXHIB	~~	MPS_spR
'
B4_fit <- cfa(B4_model, data = self_reps, missing = "ML", estimator = "MLR")
summary(B4_fit, fit.measures = TRUE, standardized = TRUE)
B4mi <- modindices(B4_fit) #modindices mostly want facet-specific error covariances

#try the scales you feel most solid about
B5_model_nocross <- '
 O =~ A*ECCPERC + A*MPS_abR
 #N =~ LOSLFEST + MPS_srR + SUICPRON
 A =~ MANIP + MPS_agR + AGG 
 C =~ IMPUL + PROPER + MPS_tdR + MPS_clR + MPS_acR + HDWK
 E =~ EXHIB + MPS_spR + MPS_scR + DETACH
'
B5_model_nocrossfit <- cfa(B5_model_nocross, data = self_reps, missing = "ML", estimator = "MLR")
summary(B5_model_nocrossfit, fit.measures = TRUE, standardized = TRUE)
B5_model_nocrossmi <- modindices(B5_model_nocrossfit) 

#3f model with no N, O factor
B3_model <- '
 E =~ EXHIB + MPS_spR + ENTITL + DETACH + MPS_scR 
 C =~ IMPUL + PROPER + MPS_tdR + MPS_clR + MPS_acR + HDWK
 ANT =~ MANIP + AGG + MPS_agR + MISTRUST + MPS_alR
 
 MPS_acR ~~ HDWK
 MISTRUST ~~ MPS_alR
 EXHIB ~~ MPS_spR
 PROPER	~~ MPS_tdR
'
B3_fit <- cfa(B3_model, data = self_reps, missing = "ML", estimator = "MLR")
summary(B3_fit, fit.measures = TRUE, standardized = TRUE)
B3_mi <- modindices(B3_fit) 
