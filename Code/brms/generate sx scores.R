library(sas7bdat)
library(psych)
k10 <- read.sas7bdat("~/Downloads/k10.sas7bdat") #%>% select(subject, k10Total) %>% rename(id = subject, k10 = k10Total)
stai <- read.sas7bdat("~/Downloads/stai.sas7bdat") #%>% select(subject, staitotal) %>% rename(id = subject, stai = staitotal)

fa_k10 <- fa(select(k10, -subject, -k10Total), nfactors = 1)
# fa_k10 <- fa(select(k10, -subject, -k10Total, -k10_1), nfactors = 1)
k10_scores <- data.frame(id = k10$subject, k10_score = fa_k10$scores) %>% rename(k10_score = MR1)

fa_stai <- fa(select(stai, -subject, -staitotal), nfactors = 1)
# fa_k10 <- fa(select(k10, -subject, -k10Total, -k10_1), nfactors = 1)
stai_scores <- data.frame(id = stai$subject, stai_score = fa_stai$scores) %>% rename(stai_score = MR1)

symptom_scores <- left_join(k10_scores, stai_scores, by = "id")
write.csv(symptom_scores, file = "~/ics/Alison/PDDDM_sx_scores.csv", row.names = FALSE)
