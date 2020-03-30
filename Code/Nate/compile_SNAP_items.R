
# compile item-level SNAP data for Tim ------------------------------------


pacman::p_load(tidyverse, sas7bdat)

snap_1 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_1.sas7bdat")
snap_2 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_2.sas7bdat")
snap_3 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_3.sas7bdat")
snap_4 <- read.sas7bdat("~/ics/Nate/PD_Inhibition_DDM/MH_Diss_Analyses/Data/snap_4.sas7bdat")

SNAP_all <-  snap_1 %>% select(-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1)
SNAP_all <- left_join(SNAP_all, select(snap_2,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")
SNAP_all <- left_join(SNAP_all, select(snap_3,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")
SNAP_all <- left_join(SNAP_all, select(snap_4,-Completed, -Time_Submitted, -IP_Address, -Last_Name, -Token, -attr1) , by = "id")

write.csv(SNAP_all, file = "~/ics/Nate/PD_Inhibition_DDM/Data/R_Originals/SNAP_all.csv",row.names = FALSE)
