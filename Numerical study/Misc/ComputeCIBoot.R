load("Results/dfParamTheta7.Rda") 
load("Results/dfParamBoot.Rda")

originalSeed <- 17
kVec <- c(50,100,200,300,400)

for(k in kVec){
  dfParamTheta7k <- dfParamTheta7[1:(k-1)]
  
  
  theta7True <- -0.3
  variation <- dfParamTheta7k[originalSeed]-dfParam
  
  q99 <- (1-0.99)/2
  l99 <- quantile(variation,q99, na.rm=TRUE)
  r99 <- quantile(variation,1-q99, na.rm=TRUE)
  
  q95 <- (1-0.95)/2
  l95 <- quantile(variation,q95, na.rm=TRUE)
  r95 <- quantile(variation,1-q95, na.rm=TRUE)
  
  
  q90 <- (1-0.90)/2
  l90 <- quantile(variation,q90, na.rm=TRUE)
  r90 <- quantile(variation,1-q90, na.rm=TRUE)
  
  
  C99 <- sum(dfParamTheta7k+l99 <= theta7True &  dfParamTheta7k+r99 >= theta7True, na.rm=TRUE) / length(na.omit(dfParamTheta7k))
  
  C95 <- sum(dfParamTheta7k+l95 <= theta7True &  dfParamTheta7k+r95 >= theta7True, na.rm=TRUE) / length(na.omit(dfParamTheta7k))
  
  C90 <- sum(dfParamTheta7k+l90 <= theta7True &  dfParamTheta7k+r90 >= theta7True, na.rm=TRUE) / length(na.omit(dfParamTheta7k))

  print(paste0("k=",k,", alpha=.99: ",round(C99*100,1) ))
  print(paste0("k=",k,", alpha=.95: ",round(C95*100,1)))
  print(paste0("k=",k,", alpha=.90: ",round(C90*100,1)))
}

