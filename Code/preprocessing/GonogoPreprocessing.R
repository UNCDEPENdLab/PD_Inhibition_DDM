
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

pacman::p_load(tidyverse, sas7bdat)

data <- read.sas7bdat(paste0(datadir, "/gonogo.sas7bdat"))


# pull one subject and figure out data structure ---------------------------



ng <- data %>% group_by(Subject, Block) %>% filter(row_number()==1) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("NoGoDisplay")) %>% ungroup()
names(ng) <- sub("NoGoDisplay_", "", names(ng))
ng$TrialType <- "NoGo"
ng$ACC <- ifelse(ng$RESP == "{SPACE}", 0, 1)


go <- data %>% group_by(Subject, Block) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("GoDisplay")) %>% ungroup()
names(go) <- sub("GoDisplay_", "", names(go))
go$TrialType <- "Go"

gng_all <- data.frame()

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



# display accuracies per subject split by condition -----------------------

gng_accs <- gng_all %>% group_by(Subject, TrialType) %>% dplyr::summarise(as.numeric(table(ACC)[2])/(as.numeric(table(ACC)[2])+as.numeric(table(ACC)[1])))
colnames(gng_accs)[3] <- "ACC"

gng_accs %>% arrange(ACC)


gng_accs_proc <- gng_all %>% group_by(Subject, TrialType, Procedure_Trial_) %>% dplyr::summarise(as.numeric(table(ACC)[2])/(as.numeric(table(ACC)[2])+as.numeric(table(ACC)[1])))
colnames(gng_accs)[4] <- "ACC"

gng_accs_proc %>% arrange(ACC)


# drop 7: 33.3% on nogo trials is substantially lower than the rest. 

gng_all$exclude <- ifelse(gng_all$Subject == 7,1,0)
  
# look at RTs across go trials --------------------------------------------

go <- gng_all %>% filter(TrialType == "Go") #%>% hist(.$RT)

hist(go$RT)
quantile(go$RT, seq(.1,1,.1))

##according to MH diss the letters were only presented for 500 ms and therefore, anything beyond that should be investigated.

go %>% arrange(-RT)


# scratch below -----------------------------------------------------------


# sample_sub <- data %>%  filter(Subject == 1)
# 
# attr(sample_sub, "column.info") <- NULL
# str(sample_sub)
# 
# nrow(sample_sub) #should equal 64*3 = 192
# unique(table(data$Subject)) #perfect.
# 
# head(sample_sub)
# 
# for(i in colnames(sample_sub)){
#   cat(i, "\n")
#   print(table(sample_sub[,i]))
# }
# 
# # ## not entirely intuitive, maybe pull one block of 64 trials and see what we can glean from its structure.
# # 
# # b1 <- sample_sub[1:64,]
# # head(b1, 32)
# # 
# # ## still strange.
# # 
# # data_scores <- read.sas7bdat(paste0(datadir, "/gonogoscores.sas7bdat"))
# # hist(data_scores$T_MeanGoRT)
# 
# 
# 
# data <- read.sas7bdat(paste0(datadir, "/gonogo.sas7bdat"))
# ng <- data %>% group_by(Subject, Block) %>% filter(row_number()==1) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("NoGoDisplay")) %>% ungroup()
# names(ng) <- sub("NoGoDisplay_", "", names(ng))
# ng$TrialType <- "NoGo"
# go <- data %>% group_by(Subject, Block) %>% select(Subject, Block, Procedure_Block_, Procedure_Trial_, starts_with("GoDisplay")) %>% ungroup()
# names(go) <- sub("GoDisplay_", "", names(go))
# go$TrialType <- "Go"
# # 
# # 
# # ## questions for MH
# # # like, so many. Better to just talk in person.
# # 
# # godisplay resp = space. correct response
# # godisplay resp = empty. ommision error
# # nogodisplay 
# # 
# # 
# # 
# # 
# # 
