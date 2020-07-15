# convert flanker summaries back to initial raw coefficients

df <- read.csv("/Users/natehall/github_repos/PD_Inhibition_DDM/Code/brms/posterior_summaries.csv")

for(i in which(df$beta_contrast == "incongruent")){
  stopifnot(df[i+1, "beta_contrast"] == "intercept")
  df[i, "v_MAP"] <- df[i, "v_MAP"] - df[i+1, "v_MAP"]
}
  
write.csv(df, file = "/Users/natehall/github_repos/PD_Inhibition_DDM/Code/brms/posterior_summaries_fix_incong_contrast.csv", row.names = FALSE) 
  
  
