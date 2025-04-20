library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file="../Results/dfParamEst.Rda")
source("../Simulation/semi_markov_sim.R")

##Initialize parameters
tstop <- 5 #terminal time of trajectory
epsilon <- 10^(-7)
h <- 1/10
time.grid <- seq(from=0,to=tstop,by=h)
x.grid <- seq(from=-4,to=4,by=h)
nstepV <- tstop/h
nstepVSmall <- nstepV/5


namesPredErr <- c("i","MSE_Prop","MAE_Prop","MSE_Orac","MAE_Orac","MSE_App","MAE_App","MSE_Naiv1","MAE_Naiv1","MSE_Naiv2","MAE_Naiv2")
dfPredErr <- data.frame(matrix(ncol = length(namesPredErr), nrow = nrow(dfParam)))
colnames(dfPredErr) <- namesPredErr
dfPredErr$i <- 1:nrow(dfParam)

g <- c(0.8,-1.2)

mu <-  function(x, y, s, t, cov, paramVec){
  if(x == 1 && y==2){
    as.numeric(muTilde12(t,cov,paramVec[1],paramVec[2],paramVec[3]))
  }
  else if(x == 1 && y==3){
    as.numeric(muTilde13(t,cov,paramVec[4],paramVec[5],paramVec[6])) 
  }
  else if(x==2 && y==3){
    
    pX <- (1-exp(-g[1]*cov^2/2))*(1-exp(1/g[2]))
    num <- pX*exp((1-exp(paramVec[7]*(t-s)*cov^2))/(paramVec[7]*cov^2))*muTilde23(t,s,cov,paramVec[7])
    denom <- 1-pX*(1-exp( (1-exp(paramVec[7]*(t-s)*cov^2))/(paramVec[7]*cov^2)   ))
    
    res <- as.numeric(num/denom)  
    
    ifelse(is.nan(res),0,res)
    }
  else{
    error("invalid transition")
  }
}

###Differential equation solver for estimand
f <- NULL

rk4 <- function(df, a, b, f0, n) {
  
  h = (b-a)/n
  
  f[1] = f0
  for (i in 1:n) {
    s = a + h * (i-1)
    k1 = h * df(s,f[i])
    k2 = h * df(s+0.5*h, f[i]+0.5*k1)
    k3 = h * df(s+0.5*h, f[i]+0.5*k2)
    k4 = h * df(s+h, f[i]+k3)
    
    f[i+1] = f[i]+1/6*(k1+2*k2+2*k3+k4)
  }
  return(f)
}

diffVi <- function(s,t,v,haz23,x){
  d1 = haz23(s,t,x)*v[1]-1
  return(c(d1))
}

Vi <- function(t0,s,haz23,x){
  V <- rk4(function(t, v) diffVi(t,v,haz23=haz23,s=s,x=x), a=tstop, b=t0+epsilon, f0=0, nstepVSmall)
  return(tail(V,n=1))
}

diffVa <- function(t,v,haz12,haz13,haz23, x){
  d1 = (haz12(0,t,x)+haz13(0,t,x))*v[1]-(t<=tstop)*haz12(0,t,x)*Vi(t,t,haz23,x)
  return(c(d1))
}

VaVec <- function(t0,haz12,haz13,haz23,x){
  V <- rk4(function(t,v) diffVa(t,v, haz12=haz12, haz13=haz13, haz23=haz23, x=x), a=tstop, b=t0+epsilon, f0=0, nstepV)
  return(V)
}

haz12 <- function(s,t,x){mu(1,2,s,t,x,params)}
haz13 <- function(s,t,x){mu(1,3,s,t,x,params)}
haz23 <- function(s,t,x){mu(2,3,s,t,x,params)}

##Compute Va oracle values over a grid
hVa <- (tstop-0)/nstepV
Va.time.grid <- hVa*(0:nstepV)
Va.value.grid <- expand.grid(t=Va.time.grid, X=x.grid)
Va.value.grid$VaTrue <- NA

for(m in 1:length(x.grid)){
  Va.value.grid$VaTrue[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12,haz13,haz23,x.grid[m])
}


#theta proposed (remove runs that did not converge)
idxProp <- which(!is.na(dfParam$theta1))
for(i in idxProp){
  paramsProp <- c(dfParam$theta1[i],dfParam$theta2[i],dfParam$theta3[i],dfParam$theta4[i],dfParam$theta5[i],dfParam$theta6[i],dfParam$theta7[i])

  haz12Prop <- function(s,t,x){mu(1,2,s,t,x,paramsProp)}
  haz13Prop <- function(s,t,x){mu(1,3,s,t,x,paramsProp)}
  haz23Prop <- function(s,t,x){mu(2,3,s,t,x,paramsProp)} 
  
  for(m in 1:length(x.grid)){
    Va.value.grid$VaPred[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12Prop,haz13Prop,haz23Prop,x.grid[m])
  }
  
  MSE <- h^2*sum((Va.value.grid$VaPred-Va.value.grid$VaTrue)^2,na.rm=TRUE)
  MAE <- h^2*sum(abs(Va.value.grid$VaPred-Va.value.grid$VaTrue),na.rm=TRUE)
    
  dfPredErr$MSE_Prop[i] <- MSE  
  dfPredErr$MAE_Prop[i] <- MAE   
  
}

gc()


#theta oracle (remove runs that did not converge)
idxOrac <- which(!(dfParam$theta1oracle > 0))

for(i in idxOrac){
  paramsOrac <- c(dfParam$theta1oracle[i],dfParam$theta2oracle[i],dfParam$theta3oracle[i],dfParam$theta4oracle[i],dfParam$theta5oracle[i],dfParam$theta6oracle[i],dfParam$theta7oracle[i])
  
  haz12Orac <- function(s,t,x){mu(1,2,s,t,x,paramsOrac)}
  haz13Orac <- function(s,t,x){mu(1,3,s,t,x,paramsOrac)}
  haz23Orac <- function(s,t,x){mu(2,3,s,t,x,paramsOrac)} 
  
  for(m in 1:length(x.grid)){
    Va.value.grid$VaPred[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12Orac,haz13Orac,haz23Orac,x.grid[m])
  }
  
  MSE <- h^2*sum((Va.value.grid$VaPred-Va.value.grid$VaTrue)^2,na.rm=TRUE)
  MAE <- h^2*sum(abs(Va.value.grid$VaPred-Va.value.grid$VaTrue),na.rm=TRUE)
  
  dfPredErr$MSE_Orac[i] <- MSE  
  dfPredErr$MAE_Orac[i] <- MAE   
  
}

gc()


#theta approx
for(i in 1:nrow(dfParam)){
  paramsApp <- c(dfParam$theta1app[i],dfParam$theta2app[i],dfParam$theta3app[i],dfParam$theta4app[i],dfParam$theta5app[i],dfParam$theta6app[i],dfParam$theta7app[i])
  
  haz12App <- function(s,t,x){mu(1,2,s,t,x,paramsApp)}
  haz13App <- function(s,t,x){mu(1,3,s,t,x,paramsApp)}
  haz23App <- function(s,t,x){mu(2,3,s,t,x,paramsApp)} 
  
  for(m in 1:length(x.grid)){
    Va.value.grid$VaPred[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12App,haz13App,haz23App,x.grid[m])
  }
  
  MSE <- h^2*sum((Va.value.grid$VaPred-Va.value.grid$VaTrue)^2,na.rm=TRUE)
  MAE <- h^2*sum(abs(Va.value.grid$VaPred-Va.value.grid$VaTrue),na.rm=TRUE)
  
  dfPredErr$MSE_App[i] <- MSE  
  dfPredErr$MAE_App[i] <- MAE   
  
}

gc()


#theta Naive 1 (remove runs that did not converge)
idxCens0 <- which(!(dfParam$theta1cens0 > 0))

for(i in idxCens0){
  paramsCens0 <- c(dfParam$theta1cens0[i],dfParam$theta2cens0[i],dfParam$theta3cens0[i],dfParam$theta4cens0[i],dfParam$theta5cens0[i],dfParam$theta6cens0[i],dfParam$theta7cens0[i])
  
  haz12Cens0 <- function(s,t,x){mu(1,2,s,t,x,paramsCens0)}
  haz13Cens0 <- function(s,t,x){mu(1,3,s,t,x,paramsCens0)}
  haz23Cens0 <- function(s,t,x){mu(2,3,s,t,x,paramsCens0)} 
  
  for(m in 1:length(x.grid)){
    Va.value.grid$VaPred[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12Cens0,haz13Cens0,haz23Cens0,x.grid[m])
  }
  
  MSE <- h^2*sum((Va.value.grid$VaPred-Va.value.grid$VaTrue)^2,na.rm=TRUE)
  MAE <- h^2*sum(abs(Va.value.grid$VaPred-Va.value.grid$VaTrue),na.rm=TRUE)
  
  dfPredErr$MSE_Naiv1[i] <- MSE  
  dfPredErr$MAE_Naiv1[i] <- MAE   
  
}

gc()


#theta Naive 2 (remove runs that did not converge)
idxCens1 <- which(!(dfParam$theta1cens1 > 0))

for(i in idxCens1){
  paramsCens1 <- c(dfParam$theta1cens1[i],dfParam$theta2cens1[i],dfParam$theta3cens1[i],dfParam$theta4cens1[i],dfParam$theta5cens1[i],dfParam$theta6cens1[i],dfParam$theta7cens1[i])
  
  haz12Cens1 <- function(s,t,x){mu(1,2,s,t,x,paramsCens1)}
  haz13Cens1 <- function(s,t,x){mu(1,3,s,t,x,paramsCens1)}
  haz23Cens1 <- function(s,t,x){mu(2,3,s,t,x,paramsCens1)} 
  
  for(m in 1:length(x.grid)){
    Va.value.grid$VaPred[(1+(m-1)*length(Va.time.grid)):(m*length(Va.time.grid))] <- VaVec(0,haz12Cens1,haz13Cens1,haz23Cens1,x.grid[m])
  }
  
  MSE <- h^2*sum((Va.value.grid$VaPred-Va.value.grid$VaTrue)^2,na.rm=TRUE)
  MAE <- h^2*sum(abs(Va.value.grid$VaPred-Va.value.grid$VaTrue),na.rm=TRUE)
  
  dfPredErr$MSE_Naiv2[i] <- MSE  
  dfPredErr$MAE_Naiv2[i] <- MAE   
  
}

gc()

##save prediction errors
#save(dfPredErr, file = "../Results/dfPredErr.Rda") 

colMeans(dfPredErr %>% select(-c(i)), na.rm = TRUE)
