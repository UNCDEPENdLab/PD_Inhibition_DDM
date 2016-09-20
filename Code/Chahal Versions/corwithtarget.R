corwithtarget <- function(df, omit=NULL, target, withvars=NULL, pmin=NULL) {
  if (!is.null(omit)) {
    dnames <- which(names(df) %in% omit)
    df <- df[,-1*dnames]
  }



  if (is.null(withvars)) {
    withvars <- names(df)[which(!names(df) %in% target)]
  }



  res <- sapply(target, function(tv) {
        cvec <- sapply(withvars, function(wv) {
              tryCatch(rc <- Hmisc::rcorr(df[,tv], df[,wv]), error=function(e) { print(e); browser() } )
              list(r=plyr::round_any(rc$r[1,2], .001), p=plyr::round_any(rc$P[1,2], .001))
            }
        )



        if (!is.null(pmin)) {
          sigr <- which(unlist(cvec["p",]) <= pmin)
          if (length(sigr) == 0L) { cvec <- c()
          } else { cvec <- cvec[,sigr, drop=FALSE] }
        }
        return(cvec)



        #print(cvec)
      }, simplify=FALSE)



  return(res)
}