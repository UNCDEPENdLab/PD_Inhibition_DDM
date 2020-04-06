## abandoned function that looks to identify subjects that endorse uncommon responses often. Lower scores for each subject denote that they frequently chose a response that others did not choose.


percentage_checks <- function(dff, inds, plot_hist = FALSE){
  df <- dff %>% select(-inds) #strip indices for easier use.
  
  perc_df <- data.frame(subj = dff[,inds])
  # perc_df <- NA
  # perc.rank <- function(x, xo)  length(x[x <= xo])/length(x)*100
  for (i in 1:ncol(df)) {
    
    perc_1 <- table(df[,i])[1]/nrow(df)
    perc_2 <- table(df[,i])[2]/nrow(df)
    
    # perc_1 <- perc.rank(df[,i],1)
    # perc_2 <- perc.rank(df[,i],2)
    # 
    percs <- data.frame(ifelse(df[,i] == 1, perc_1, perc_2))
    colnames(percs) <- names(df)[i]
    perc_df <- cbind(perc_df, percs)
  }
  
  
  # range(perc_df[,-1])
  
  perc_df$perc_sums <- rowSums(perc_df[,-1])
  
  perc_df <- perc_df %>% arrange(perc_sums)
  # hist(perc_df$perc_sums)
  
  if(plot_hist){
    gg <- ggplot(perc_df, aes(x = perc_sums)) + geom_histogram(); print(gg)
  }
  
  
  return(perc_df)
}





#scrarch below, including more advanced mahalanobis distance checks from various packages 


# cov(df)
# dff <- mpq %>% select(subject, starts_with("MPQ_"))
# 
# # inds <- unique(df$subject)
# inds <- "subject"
# 
# 
# 
# 
# 
# 
# #mahalanobis distance wrapper function
# mah <- function(x, cx = NULL) {
#   browser()
#   if(is.null(cx)) cx <- cov(x)
#   out <- lapply(1:nrow(x), function(i) {
#     mahalanobis(x = x, 
#                 center = do.call("c", x[i, ]),
#                 cov = cx)
#   })
#   return(as.dist(do.call("rbind", out)))
# }
# 
# 
# x <- data.frame(X = c(rnorm(10, 0), rnorm(10, 5)), 
#                 Y = c(rnorm(10, 0), rnorm(10, 7)), 
#                 Z = c(rnorm(10, 0), rnorm(10, 12)))
# rownames(x) <- LETTERS[1:20]
# plot(x, pch = LETTERS[1:20])
# 
# d <- mah(x)
# d
# 
# 
# 
# library(assertr)
# 
# x <- maha_dist(mtcars)
# 
# y <- maha_dist(iris)#, robust=TRUE)
# z <- maha_dist(iris, robust=TRUE)
# cor(y,z)
# 
# hist(y)
# hist(z)
# 
# library(magrittr)            # for piping operator
# library(dplyr)               # for "everything()" function
# 
# # using every column from mtcars, compute mahalanobis distance
# # for each observation, and ensure that each distance is within 10
# # median absolute deviations from the median
# mtcars %>%
#   insist_rows(maha_dist, within_n_mads(10), everything())
# 
# maha_dist(df[,1:75], robust = TRUE)
# 
# 
# 
# quantile()
# 
# # require(graphics)
# # 
# # ma <- cbind(1:6, 1:3)
# # (S <-  var(ma))
# # mahalanobis(c(0, 0), 1:2, S)
# # 
# # x <- matrix(rnorm(100*3), ncol = 3)
# # stopifnot(mahalanobis(x, 0, diag(ncol(x))) == rowSums(x*x))
# # ##- Here, D^2 = usual squared Euclidean distances
# # 
# # Sx <- cov(x)
# # D2 <- mahalanobis(x, colMeans(x), Sx)
# # plot(density(D2, bw = 0.5),
# #      main="Squared Mahalanobis distances, n=100, p=3") ; rug(D2)
# # qqplot(qchisq(ppoints(100), df = 3), D2,
# #        main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
# #                            " vs. quantiles of" * ~ chi[3]^2))
# # abline(0, 1, col = 'gray')