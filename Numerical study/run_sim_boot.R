source("Simulation/semi_markov_sim.R")
source("Simulation/semi_markov_sim_adj.R")
library("dplyr")
library("stats")
library("pracma")


##Initialize parameters
originalSeed <- 17
set.seed(originalSeed)
nsim <- 1500 #number of simulations 
X <- runif(nsim,min=-4,max=4) #simulate X
tstop <- 5 #terminal time of trajectory 
V <- runif(nsim,min=0,max=1)
C <- runif(nsim,min=V,max=tstop)
nmax <- 2 #maximum number of jumps (used for memory initialization)

namesData <- c("id","X","V","C","tn","tnNext","yn","ynNext","delayReport", "tnAdj1" , "tnAdj2", "tnAdj3")
dataRaw <- data.frame(matrix(ncol = length(namesData), nrow = nsim*nmax))
names(dataRaw) <- namesData

###Simulate original data
weibullTransform <- function(w,X,lambda,k,beta){
  return(1/lambda*(-log(1-w^{exp(-X*beta)}))^{1/k})
}

simReportingDelay <- function(state, X){
  W <- runif(n=1,min=0,max=1)
  if(state == 1){
    return(weibullTransform(W,X,2,0.5,0.1))
  }
  else if(state == 2){
    return(weibullTransform(W,X,1,1.5,0.2))
  }
}

k <- 1
for(i in 1:nsim){
  Xi <- X[i]
  Vi <- V[i]
  Ci <- C[i]
  simi <- simulator(1,Vi,Vi,Ci,Xi)
  
  for(j in 1:nmax){
    
    repDelay <- 0
    
    if(!is.na(simi$marks[j+1])){
      if(simi$marks[j] == 1 & simi$marks[j+1] == 3){
        repDelay <- simReportingDelay(1,Xi)
      }
      else if(simi$marks[j] == 2 & simi$marks[j+1] == 3){
        repDelay <- simReportingDelay(2,Xi)
      }
    }
    
    adjT0 <- NA
    adjT1 <- NA
    adjT2 <- NA
    
    if(!is.na(simi$marks[j+1]) & simi$marks[j]==2 & simi$marks[j+1]==3 & simi$times[j+1]+repDelay < tstop){
      adj_sim <- simulator_adj(1, simi$times[j+1]+repDelay, simi$times[j+1]+repDelay, tstop, Xi)
      adjT0 <- adj_sim$times[1]
      adjT1 <- adj_sim$times[2]
      adjT2 <- adj_sim$times[3]
    }
    
    dataRaw[k,] <- c(i,Xi,Vi,Ci,
                     simi$times[j],
                     simi$times[j+1],
                     simi$marks[j],
                     simi$marks[j+1],
                     repDelay,
                     adjT0,adjT1,adjT2) 
    k <- k+1
  }
}
dataRaw <- dataRaw %>% filter(!is.na(tn))  #remove "all-NA" rows

data <- dataRaw

data <- data %>% mutate(thinned = (tnNext+delayReport > tstop) ) #flag those jumps that are deleted due to delays

dataOrig <- data


###Functions

##Calculate absorption probability
absProb <- function(state,dur,X,p1adj,p2adj){
  if(state == 1){
    return((1-exp(-p1adj*X^2/(dur+2)))*(1-exp(1/p2adj)))
  } else if(state == 2){
    return(1-exp(1/p2adj*exp(p2adj*dur)))
  } else if(state == 3){
    return(1)
  }
}
absProbVec <- Vectorize(absProb)


##Estimate adjudication
loglikeAdj12 <- function(data,params){
  p <- params[1]
  
  loglike <- data %>% 
    filter(yn == 2, ynNext == 3, thinned == FALSE) %>% 
    mutate(tnAdj1 = tnNext+delayReport, tend = ifelse(is.na(tnAdj2),tstop,tnAdj2) ) %>% 
    mutate(dur = tend-tnAdj1, 
           occTerm = ifelse(is.na(tnAdj2),0,log(p*(X/(dur+2))^2)), 
           expoTerm = p*dur*X^2/(2*dur+4) ) %>%
    summarise(sum(occTerm-expoTerm)) %>%
    as.numeric()
  
  return(loglike)
}

loglikeAdj23 <- function(data,params){
  p <- params[1]
  
  loglike <- data %>% 
    filter(yn == 2, ynNext == 3, thinned == FALSE, !is.na(tnAdj2)) %>% 
    mutate(tend = ifelse(is.na(tnAdj3),tstop,tnAdj3) ) %>% 
    mutate(durEnd = tend-tnAdj2, 
           occTerm = ifelse(is.na(tnAdj3),0,log(exp(p*durEnd)) ), 
           expoTerm = 1/p*(exp(p*durEnd)-exp(p*0) ) ) %>% 
    summarise(sum(occTerm-expoTerm)) %>%
    as.numeric()
  
  return(loglike)
}

##Estimate reporting delay and event hazards
alphaDelay <- function(t,X,lambda,k,beta){
  alpha0 = k*lambda^k*t^(k-1)/(exp((lambda*t)^k)-1)
  alpha = alpha0*exp(X*beta)
  return(alpha)
}

loglikeRep13 <- function(data,params){
  lambda <- params[1]
  k <- params[2]
  beta <- params[3]
  
  loglike <- data %>% 
    filter(yn == 1, ynNext == 3, thinned == FALSE) %>% 
    mutate(occTerm = log(alphaDelay(delayReport,X,lambda,k,beta)), 
           expoTerm = -exp(X*beta)*log((1-exp(-(lambda*delayReport)^k))/(1-exp(-(lambda*(tstop-tnNext))^k))),
           absProb = 1) %>%
    summarise(sum((occTerm-expoTerm)*absProb)) %>% 
    as.numeric()
  
  return(loglike)
}

loglikeRep23 <- function(data,p1adj,p2adj,params){
  lambda <- params[1]
  k <- params[2]
  beta <- params[3]
  
  loglike <- data %>% 
    filter(yn == 2, ynNext == 3, thinned == FALSE) %>% 
    mutate(occTerm = log(alphaDelay(delayReport,X,lambda,k,beta)), 
           expoTerm = -exp(X*beta)*log((1-exp(-(lambda*delayReport)^k))/(1-exp(-(lambda*(tstop-tnNext))^k))),
           adjState = ifelse(!is.na(tnAdj3),3,ifelse(!is.na(tnAdj2),2,1)),
           adjDur = ifelse(adjState==3,0,ifelse(adjState==2,tstop-tnAdj2,tstop-tnAdj1)),
           absProb = absProbVec(adjState,adjDur,X,p1adj,p2adj)) %>% 
    summarise(sum((occTerm-expoTerm)*absProb)) %>% 
    as.numeric()
  
  return(loglike)
}


loglike12app <- function(data,params){
  for(i in 1:length(params)){assign(paste0("p",i),params[i])}
  
  loglike <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned13 = sum(yn==1 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned13 = ifelse(is.na(anyThinned13),0,anyThinned13)) %>%
    filter(yn == 1) %>% 
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned13 == 0, tnNext, C), 
           occTerm = ifelse(!is.na(ynNext) & ynNext==2, log(muTilde12(tnNext,X,p1,p2,p3)), 0), 
           expoTerm = 1/p2*(muTilde12(tend,X,p1,p2,p3)-muTilde12(tn,X,p1,p2,p3)) ) %>%
    summarise(sum(occTerm-expoTerm)) %>%
    as.numeric()
  
  return(loglike)
}


delayDist13 <- function(t,X,lambda,k,beta) (1-exp(-(lambda*(tstop-t))^k))^(exp(X*beta))

nu13 <- function(t,X,p1,p2,p3,lambda,k,beta){
  delayDist <- delayDist13(t,X,lambda,k,beta)
  return(muTilde13(t,X,p1,p2,p3)*delayDist)
}

loglike13app <- function(data,p13delay,params){
  for(i in 1:length(params)){assign(paste0("p",i),params[i])}
  
  lambda <- p13delay[1]
  k <- p13delay[2]
  beta <- p13delay[3]
  
  loglike <- data %>%
    mutate(occTerm = ifelse(!is.na(ynNext) & anyThinned13 == 0 & ynNext==3, log(muTilde13(tnNext,X,p1,p2,p3)), 0)) %>%
    rowwise() %>%
    mutate(expoTerm = integrate(f=nu13, lower=tn, upper=tend, X=X, p1=p1,p2=p2,p3=p3,lambda=lambda,k=k,beta=beta)$value ) %>% 
    ungroup() %>%
    summarise(sum(occTerm-expoTerm)) %>%
    as.numeric()
  
  return(loglike)
}

delayDist23 <- function(t,X,lambda,k,beta) (1-exp(-(lambda*(tstop-t))^k))^(exp(X*beta))

mu23Ast <- function(t,s,X,p1){
  pX <- (1-exp(-0.8*X^2/2))*(1-exp(-1/1.2))
  Vt <- t-s
  denom <- 1-pX*(1-exp(1/(p1*X^2)*(1-exp(p1*Vt*X^2))))
  num <- pX*muTilde23(t,s,X,p1)*exp(1/(p1*X^2)*(1-exp(p1*Vt*X^2)))
  muAst23 <- num/denom
  return(muAst23)
}

nu23 <- function(t,s,X,p1,lambda,k,beta){
  muAst23 <- mu23Ast(t,s,X,p1)
  delayDist <- delayDist23(t,X,lambda,k,beta)
  return(muAst23*delayDist)
}

loglike23app <- function(data,p23delay,p1adj,p2adj,params){
  for(i in 1:length(params)){assign(paste0("p",i),params[i])}
  
  lambda <- p23delay[1]
  k <- p23delay[2]
  beta <- p23delay[3]
  
  loglike <- data %>% 
    mutate(occTermW1 = ifelse(!is.na(ynNext) & anyThinned23 == 0, log(nu23(tnNext,tn,X,p1,lambda,k,beta)), 0),
           occTermW0 = 0) %>%
    rowwise() %>%
    mutate(expoTermW1 = integrate(f=nu23, lower=tn, upper=tend, s=tn, X=X, p1=p1,lambda=lambda,k=k,beta=beta)$value ) %>% 
    mutate(expoTermW0 = integrate(f=nu23, lower=tn, upper=C, s=tn, X=X, p1=p1,lambda=lambda,k=k,beta=beta)$value ) %>%
    ungroup() %>%
    mutate(adjState = ifelse(!is.na(tnAdj3),3,ifelse(!is.na(tnAdj2),2,1)),
           adjDur = ifelse(adjState==3,0,ifelse(adjState==2,tstop-tnAdj2,tstop-tnAdj1)),
           absProb = ifelse(!is.na(tnAdj1),absProbVec(adjState,adjDur,X,p1adj,p2adj),0)) %>%
    summarise(sum( (occTermW1-expoTermW1)*absProb+(occTermW0-expoTermW0)*(1-absProb) )) %>%
    as.numeric()
  
  return(loglike)
}

P22 <- function(t,s,X,p1){
  integral <- sapply(t, function(tt) integrate(f=mu23Ast, lower=s, upper=tt, s=s, X=X, p1=p1)$value)
  return(exp(-integral))
}

denomIntegrand2 <- function(t,s,X,p23delay,p1){
  return(P22(t,s,X,p1)*nu23(t,s,X,p1,p23delay[1],p23delay[2],p23delay[3]))
}

mu2Factor <- function(t,s,X,p23delay,p1){
  numerator <- P22(t,s,X,p1) 
  denominator <- 1-sapply(t, function(tt) integrate(f=denomIntegrand2, lower=s, upper=tt, s=s, X=X, p23delay=p23delay, p1=p1)$value) 
  
  result <- numerator/denominator 
  result <- min(result,1) #for numerical stability
  result <- max(result,numerator) #for numerical stability
  
  return(numerator/denominator)
}

mu23 <- function(t,s,X,p23delay,p1){
  nu <- nu23(t,s,X,p1,p23delay[1],p23delay[2],p23delay[3])
  factor <- mu2Factor(t,s,X,p23delay,p1)
  
  return(factor*nu)
}

loglike2 <- function(data,p23delay,p1adj,p2adj,params){
  for(i in 1:length(params)){assign(paste0("p",i),params[i])}
  
  #restrict search space for numerical stability
  if(p1 > 0 | p1 < -1.5 ){return(-Inf)}
  
  loglike <- data %>% 
    rowwise() %>%
    mutate(occTermW1 = ifelse(!is.na(ynNext) & anyThinned23 == 0, log(mu23(tnNext,tn,X,p23delay,p1)), 0),
           occTermW0 = 0,
           expoTermW1 = integrate(f=mu23, lower=tn, upper=tend, s=tn, X=X, p1=p1,p23delay=p23delay)$value, 
           expoTermW0 = integrate(f=mu23, lower=tn, upper=C, s=tn, X=X, p1=p1,p23delay=p23delay)$value ) %>%
    ungroup() %>%
    mutate(adjState = ifelse(!is.na(tnAdj3),3,ifelse(!is.na(tnAdj2),2,1)),
           adjDur = ifelse(adjState==3,0,ifelse(adjState==2,tstop-tnAdj2,tstop-tnAdj1)),
           absProb = ifelse(!is.na(tnAdj1),absProbVec(adjState,adjDur,X,p1adj,p2adj),0)) %>%
    summarise(sum( (occTermW1-expoTermW1)*absProb+(occTermW0-expoTermW0)*(1-absProb) )) %>%
    as.numeric()
  
  return(loglike)
}



runEstimation <- function(data) {
  
  #params
  p1adj <- NA
  p2adj <- NA
  p23delay <- NA
  theta2 <- NA
  
  tryCatch({p1adj <- optim(par = c(1), fn = loglikeAdj12, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  tryCatch({p2adj <- optim(par = c(-1), fn = loglikeAdj23, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  tryCatch({p23delay <- optim(par = c(1,1,1), fn = loglikeRep23, data=data, p1adj=p1adj, p2adj=p2adj, control=list(fnscale=-1))$par}, error = function(e){})
  
  dataloglike2 <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned23 = sum(yn==2 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned23 = ifelse(is.na(anyThinned23),0,anyThinned23)) %>%
    filter(yn == 2) %>%
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned23 == 0, tnNext, C)) %>%
    mutate(occTerm = NA, expoTerm = NA)
  
  tryCatch({theta2 <- optim(par = c(-0.25), fn = loglike2, data=dataloglike2, p23delay=p23delay, p1adj=p1adj, p2adj=p2adj, control=list(fnscale=-1), method="BFGS")$par}, error = function(e){})
  
  return(theta2)
  
}


runSim <- function(){
  totalSim <- 1000
  startSim <- 1
  
  dataParam <- c()
  for(i in startSim:totalSim){
    set.seed(originalSeed+i)
    resulti <- NA
    
    data <- data.frame(matrix(ncol = ncol(dataOrig), nrow = 0))
    names(data) <- names(dataOrig)
    
    ids <- sample(1:nsim,size=nsim,replace=TRUE)
    
    k <- 1
    for(j in ids){
      data <- rbind(data,dataOrig %>% filter(id==j) %>% mutate(id=k))
      k <- k+1
    }
    
    tryCatch({resulti <- runEstimation(data)}, error = function(e){})
    dataParam[i] <- resulti
    
    ##To save partial results, uncomment here:
    #if( (i <= startSim+9) | (i %% 20 == 0) ){
    #  save(dataParam, file = paste0("Results/dfParamBoot",i,".Rda")) 
    #}
    
  }
  
  return(dataParam)
}

dataParam <- runSim()
save(dataParam, file = "Results/dfParamBoot.Rda") 