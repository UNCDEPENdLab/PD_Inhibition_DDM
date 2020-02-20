##scratch

cong_check <- ifelse(flanker$Central == flanker$Flanker,0,1 )
table(cong_check == flanker$Incongruent)
table(flanker$Procedure)
table(flanker$CongruentBlock)
table(flanker$TrialSlide_RESP)

head(data.frame(flanker))

Cdiff <- flanker$TrialSlide_RTTime - flanker$TrialSlide_OnsetTime
head(diff)
head(flanker$TrialSlide_RT)

table(flanker$TrialSlide_DurationError) #not sure what this means

table(flanker$Block)
length(table(flanker$Block)) #trial number?

table(flanker$Subject) 
length(table(flanker$Subject)) #okay so these check out

acc_check <- ifelse(flanker$TrialSlide_CRESP == flanker$TrialSlide_RESP, 1,0)
table(acc_check)

table(flanker$TrialSlide_ACC)

length(acc_check)
length(flanker$TrialSlide_ACC)

flanker$Correct


data.frame(flanker[which(is.na(flanker$TrialSlide_ACC)),"Flanker"]) == data.frame(flanker[which(is.na(flanker$TrialSlide_RT)),"Flanker"]) #looks like these are the same.

x <- data.frame(flanker[which(is.na(flanker$TrialSlide_ACC)),])
table(x$Subject)
hist(table(x$Subject))


flank_na.rm <- flank[-which(is.na(flank$TrialSlide_ACC)),]
# which(is.na(flank_na.rm$TrialSlide_ACC))

table(flanker$CongruentTrials_Cycle)
length(table(flanker$CongruentTrials_Sample))
table(flanker$CongruentTrials)

table(flanker$IncongruentTrials_Cycle)
length(table(flanker$IncongruentTrials_Sample))
table(flanker$CongruentTrials)


table(flanker$computer) #no idea...




# check blocking structure ------------------------------------------------

for(sub in unique(flanker$Subject)){
  df <- flanker %>% filter(Subject == sub) %>% data.frame()
  stopifnot(nrow(df) == 160)
  
  con <- df$CongruentTrials_Cycle
  con <- con[!is.nan(con)]
  
  stopifnot(all(table(con) == 40))
  
  icon <- df$CongruentTrials_Cycle
  icon <- con[!is.nan(icon)]
  
  stopifnot(all(table(icon) == 40))
  
  
  for(cyc in 1:2){
    con1 <- df %>% filter(CongruentTrials_Cycle == cyc)
    t <- as.numeric(table(con1$Incongruent))
    stopifnot(t[1]/(sum(t)) == 0.7)
    
    con1 <- df %>% filter(IncongruentTrials_Cycle == cyc)
    t <- as.numeric(table(con1$Incongruent))
    stopifnot(t[1]/(sum(t)) == 0.3)
  }
  
            
  
}



# check accuracies --------------------------------------------------------

head(flank)

summaries_tot <- flank %>% group_by(Subject) %>% summarise(Subject_ACC = as.numeric(table(TrialSlide_ACC)[2])/(as.numeric(table(TrialSlide_ACC)[2])+as.numeric(table(TrialSlide_ACC)[1]))) 
summaries_con <- flank %>% group_by(Subject, Incongruent)%>% summarise(Subject_ACC_con = as.numeric(table(TrialSlide_ACC)[2])/(as.numeric(table(TrialSlide_ACC)[2])+as.numeric(table(TrialSlide_ACC)[1]))) 

# commented to not muck up the makrdown but if interested can uncomment and run.

for (i in which(is.na(summaries_tot$Subject_ACC))) {
  df <- summaries_tot[i,]
  
  sub <- df$Subject
  # con <- df$Incongruent
  
  s <- flank %>% filter(Subject == sub)
  print(table(s[,c("TrialSlide_ACC")]))
  
  
}

# commented to not muck up the makrdown but if interested can uncomment and run.

for (i in which(is.na(summaries$Subject_ACC_con))) {
  df <- summaries[i,]

  sub <- df$Subject
  con <- df$Incongruent

  s <- flank %>% filter(Subject == sub)
  print(table(s[,c("TrialSlide_ACC", "Incongruent")]))


}

s %>% #group_by(Subject) %>% 
  group_by(Subject, Incongruent)%>% summarise(Subject_ACC_con = as.numeric(table(TrialSlide_ACC)[2])/(as.numeric(table(TrialSlide_ACC)[2])+as.numeric(table(TrialSlide_ACC)[1]))) 

table(s$TrialSlide_ACC)

s <- flank_na.rm %>% filter(Subject == 4)
table(s[,c("TrialSlide_ACC", "Incongruent")])

# table(s$TrialSlide_ACC)
# 
# sort(unique(flank$Subject_ACC))
hist(unique(flank$Subject_ACC))


#create column that denotes ABBA blocking structure
data.frame(flanker) %>% mutate(ABBA = ifelse(CongruentTrials_Cycle == 1, "A1", 
                                 ifelse(CongruentTrials_Cycle == 2, "A2",
                                        ifelse(IncongruentTrials_Cycle == 1, "B1", "B2")))) 


ggplot(SNAPz, aes_string(x = s)) + geom_histogram() 



x <- seq(1,5)
x2 <- x^2
x2[2] <- NA
y <- x - 1
X <- data.frame(x,x2,y)
harmonic.mean(x)
harmonic.mean(x2)
harmonic.mean(X)
harmonic.mean(X,na.rm=FALSE)
harmonic.mean(X,zero=FALSE)

sdhm <- function(x){
  xl <- length(x) - length(which(is.na(x)))
  x <- 
  out <- sqrt(mean(1/x, na.rm = TRUE))^(-4)*var(1/x, na.rm = TRUE)/length(xl)
  return(out )
}

s <- 102

sdhm <- function(x){
  x <- as.numeric(na.omit(x))
  out <- sqrt((mean(1/x))^(-4)*var(1/x)/length(x))
}


# str(flank_acc_trim)
flank_acc_trim <- flank_acc_trim %>% rename(`rt` = `TrialSlide_RT`)

flank_acc_trim <- flank_acc_trim %>% group_by(Subject) %>% mutate(Amean_rt = mean(rt, na.rm = TRUE),
                                                                  Hmean_rt = harmonic.mean(rt, na.rm = TRUE),
                                                                  SD_rt = sd(rt, na.rm = TRUE),
                                                                  SDH_rt = sdhm(rt)) %>% mutate(
                                                                    SD3Amean_rt_high = Amean_rt + 3*SD_rt,
                                                                    SD3Amean_rt_low = Amean_rt - 3*SD_rt,
                                                                    SD3Hmean_rt_high = Hmean_rt + 3*SDH_rt,
                                                                    SD3Hmean_rt_low = Hmean_rt - 3*SDH_rt
                                                                  ) %>% group_by(Subject, Incongruent) %>% mutate(con_Amean_rt = mean(rt, na.rm = TRUE),
                                                                                                                  con_Hmean_rt = harmonic.mean(rt, na.rm = TRUE),
                                                                                                                  con_SD_rt = sd(rt, na.rm = TRUE),
                                                                                                                  con_SDH_rt = sdhm(rt)) %>% mutate(
                                                                                                                    con_SD3Amean_rt_high = con_Amean_rt + 3*con_SD_rt,
                                                                                                                    con_SD3Amean_rt_low = con_Amean_rt - 3*con_SD_rt,
                                                                                                                    con_SD3Hmean_rt_high = con_Hmean_rt + 3*con_SDH_rt,
                                                                                                                    con_SD3Hmean_rt_low = con_Hmean_rt - 3*con_SDH_rt)






flank_acc_trim <- flank_acc_trim %>% mutate(group_Amean_rt = mean(rt, na.rm = TRUE),
                                            group_Hmean_rt = harmonic.mean(rt, na.rm = TRUE),
                                            group_SD_rt = sd(rt, na.rm = TRUE),) %>% mutate(group_SD3Amean_rt_high = group_Amean_rt + 3*SD_rt,
                                                                                            group_SD3Amean_rt_low = group_Amean_rt - 3*SD_rt,
                                                                                            group_SD3AHmean_rt_high = group_Hmean_rt + 3*SD_rt,
                                                                                            group_SD3AHmean_rt_low = group_Hmean_rt - 3*SD_rt) %>% 
  group_by(Subject) %>% mutate(Amean_rt = mean(rt, na.rm = TRUE),
                               Hmean_rt = harmonic.mean(rt, na.rm = TRUE),
                               SD_rt = sd(rt, na.rm = TRUE),) %>% mutate(SD3Amean_rt_high = Amean_rt + 3*SD_rt,
                                                                         SD3Amean_rt_low = Amean_rt - 3*SD_rt,
                                                                         SD3AHmean_rt_high = Hmean_rt + 3*SD_rt,
                                                                         SD3AHmean_rt_low = Hmean_rt - 3*SD_rt) %>% 
  group_by(Subject, Incongruent) %>% mutate(con_Amean_rt = mean(rt, na.rm = TRUE),
                                            con_Hmean_rt = harmonic.mean(rt, na.rm = TRUE),
                                            con_SD_rt = sd(rt, na.rm = TRUE),) %>% mutate(con_SD3Amean_rt_high = con_Amean_rt + 3*con_SD_rt,
                                                                                          con_SD3Amean_rt_low = con_Amean_rt - 3*con_SD_rt,
                                                                                          con_SD3AHmean_rt_high = con_Hmean_rt + 3*con_SD_rt,
                                                                                          con_SD3AHmean_rt_low = con_Hmean_rt - 3*con_SD_rt)

      
                                                                                                                    
x <- ggplot(flank_acc_trim, aes(x = rt)) + geom_histogram() + 
  geom_vline(xintercept = unique(flank_acc_trim$group_Hmean_rt)) +
  geom_vline(xintercept = unique(flank_acc_trim$group_SD3AHmean_rt_low), color = "red") +
  geom_vline(xintercept = unique(flank_acc_trim$group_SD3AHmean_rt_high), color = "red") + xlim(0,1000) +
  labs(title = "Harmonic mean and 3 SD over/under, whole group")

y <- ggplot(flank_acc_trim, aes(x = rt)) + geom_histogram() + 
  geom_vline(xintercept = unique(flank_acc_trim$group_Amean_rt)) +
  geom_vline(xintercept = unique(flank_acc_trim$group_SD3Amean_rt_low), color = "red") +
  geom_vline(xintercept = unique(flank_acc_trim$group_SD3Amean_rt_high), color = "red") + xlim(0,1000) +
  labs(title = "Arithmetic mean and 3 SD over/under, whole group")

plot_grid(x,y, ncol = 1)
                                                                                                                                                                                                                   )

str(flank_acc_trim)
flank_acc_trim %>% select(ends_with("_rt")) %>% cor()
flank_acc_trim[,c(12:28)]  %>% cor()


# pdf("Figures/subject_rt_dists_raw.pdf", width = 11, height = 8)
x <- ggplot(flank_acc_trim, aes(x = rt)) + geom_histogram() + ggtitle("Group")
plot(x)
for(s in unique(flank_acc_trim$Subject)){
  df <- flank_acc_trim %>% filter(Subject == s)
  df
  x <- ggplot(df, aes(x = rt)) + geom_histogram() + ggtitle(paste0("Subject ",s))
  plot(x)
  
}
# dev.off()

dfi <- filter(df, Incongruent == 1); dfc <- filter(df, Incongruent == 0)

x_low <- max(0,min((unique(df$SD3AHmean_rt_low)-50), (unique(df$SD3Amean_rt_low) -50), (range(df$rt, na.rm = TRUE)[1]) - 50))
x_high <- min(1000,max((unique(df$SD3AHmean_rt_high)+50), (unique(df$SD3Amean_rt_high) +50), (range(df$rt, na.rm = TRUE)[2]) + 50))


by_sub_harm <- ggplot(df, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(df$Hmean_rt)) +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_low), color = "red") +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s),subtitle = "Harmonic mean and 3 SD over/under, per subject")
# plot(x)

by_sub_art <- ggplot(df, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(df$Amean_rt)) +
  geom_vline(xintercept = unique(df$SD3Amean_rt_low), color = "red") +
  geom_vline(xintercept = unique(df$SD3Amean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s),subtitle = "Arimetic mean and 3 SD over/under, per subject")

by_sub_incon_harm <- ggplot(dfi, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(dfi$Hmean_rt), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_low), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_high), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(dfi$con_Hmean_rt), color = "red") +
  geom_vline(xintercept = unique(dfi$con_SD3AHmean_rt_low), color = "red") +
  geom_vline(xintercept = unique(dfi$con_SD3AHmean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s, " incongruent trials"),subtitle = "Harmonic mean and 3 SD over/under, per subjectfor ONLY incongruent trials")

by_sub_incon_art <- ggplot(dfi, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(dfi$Amean_rt), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3Amean_rt_low), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3Amean_rt_high), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(dfi$con_Amean_rt), color = "red") +
  geom_vline(xintercept = unique(dfi$con_SD3Amean_rt_low), color = "red") +
  geom_vline(xintercept = unique(dfi$con_SD3Amean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s, " incongruent trials"),subtitle = "Arithmetic mean and 3 SD over/under, per subject for ONLY incongruent trials")


by_sub_con_harm <- ggplot(dfc, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(dfc$Hmean_rt), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_low), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3AHmean_rt_high), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(dfc$con_Hmean_rt), color = "red") +
  geom_vline(xintercept = unique(dfc$con_SD3AHmean_rt_low), color = "red") +
  geom_vline(xintercept = unique(dfc$con_SD3AHmean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s, " congruent trials"),subtitle = "Harmonic mean and 3 SD over/under, per subject for ONLY congruent trials")

by_sub_con_art <- ggplot(dfc, aes(x = rt)) + geom_histogram(bins = 40) + 
  geom_vline(xintercept = unique(dfc$Amean_rt), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3Amean_rt_low), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(df$SD3Amean_rt_high), color = "gray", linetype = "dotted") +
  geom_vline(xintercept = unique(dfc$con_Amean_rt), color = "red") +
  geom_vline(xintercept = unique(dfc$con_SD3Amean_rt_low), color = "red") +
  geom_vline(xintercept = unique(dfc$con_SD3Amean_rt_high), color = "red") + xlim(x_low, x_high) +
  labs(title = paste0("Subject ",s, " congruent trials"),subtitle = "Arithmetic mean and 3 SD over/under, per subject for ONLY congruent trials")


plot_grid(by_sub_harm,  by_sub_con_harm, by_sub_con_art, by_sub_art,by_sub_incon_harm, by_sub_incon_art)


# plot(by_sub_incon_harm)
# 
# 
# 
# # y <- ggplot(df, aes(x = rt)) + geom_histogram(bins = 40) + 
# #   geom_vline(xintercept = unique(df$Amean_rt)) +
# #   geom_vline(xintercept = unique(df$SD3Amean_rt_low), color = "red") +
# #   geom_vline(xintercept = unique(df$SD3Amean_rt_high), color = "red") + xlim(x_low, x_high) +
# #   labs(title = paste0("Subject ",s),subtitle = "Arimetic mean and 3 SD over/under, per subject")
# 
# 
# 
# 
# cowplot::plot_grid(x,y, ncol = 1)





flank_acc_trim <- data.frame(flank_acc_trim)
ggplot(data = flank_acc_trim, aes(x = con_SD3AHmean_rt_high)) + geom_histogram()


