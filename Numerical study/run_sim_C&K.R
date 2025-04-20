source("Simulation/semi_markov_sim.R")
source("Simulation/semi_markov_sim_adj.R")
library("dplyr")
library("stats")
library("pracma")

runSimSeed <- function(seed) {
  
  ##Initialize parameters
  set.seed(seed)
  nsim <- 1500 #number of simulations 
  X <- runif(nsim,min=-4,max=4) #simulate X
  tstop <- 5 #terminal time of trajectory 
  V <- runif(nsim,min=0,max=1)
  C <- runif(nsim,min=V,max=tstop)
  nmax <- 2 #maximum number of jumps (used for memory initialization)
  
  namesData <- c("id","X","V","C","tn","tnNext","yn","ynNext","delayReport", "tnAdj1" , "tnAdj2", "tnAdj3")
  dataRaw <- data.frame(matrix(ncol = length(namesData), nrow = nsim*nmax))
  names(dataRaw) <- namesData
  
  #params
  p23delay <- NA
  theta23app <- NA
  
  ##Simulate and save observations
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
  
  ##Estimate adjudication Cook&Kosorok
  
  dataAdj <- data %>% filter(!is.na(tnAdj1)) %>% mutate(adjOutc=ifelse(is.na(tnAdj3),0,1))
  
  tryCatch({fitAdj <- glm(adjOutc ~ X + tnNext + delayReport, data = dataAdj, family = "binomial")$coefficients}, error = function(e){})
  
  absProb <- function(x,tnnext,delayreport){
    linpred <- fitAdj[1]+fitAdj[2]*x+fitAdj[3]*tnnext+fitAdj[4]*delayreport
    return(exp(linpred)/(1+exp(linpred)))
  }
  
  absProbVec <- Vectorize(absProb)
  
  ##Estimate reporting delay
  alphaDelay <- function(t,X,lambda,k,beta){
    alpha0 = k*lambda^k*t^(k-1)/(exp((lambda*t)^k)-1)
    alpha = alpha0*exp(X*beta)
    return(alpha)
  }
  
  loglikeRep23 <- function(data,params){
    lambda <- params[1]
    k <- params[2]
    beta <- params[3]
    
    loglike <- data %>% 
      filter(yn == 2, ynNext == 3, thinned == FALSE) %>% 
      mutate(occTerm = log(alphaDelay(delayReport,X,lambda,k,beta)), 
             expoTerm = -exp(X*beta)*log((1-exp(-(lambda*delayReport)^k))/(1-exp(-(lambda*(tstop-tnNext))^k))),
             adjState = ifelse(!is.na(tnAdj3),3,ifelse(!is.na(tnAdj2),2,1)),
             adjDur = ifelse(adjState==3,0,ifelse(adjState==2,tstop-tnAdj2,tstop-tnAdj1)),
             absProb = absProbVec(X,tnNext,delayReport)) %>% 
      summarise(sum((occTerm-expoTerm)*absProb)) %>% 
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({p23delay <- optim(par = c(1,1,1), fn = loglikeRep23, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  ##Estimate hazard
  
  #model based
  
  #Poisson approximation
  
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
  
  dataloglike23 <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned23 = sum(yn==2 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned23 = ifelse(is.na(anyThinned23),0,anyThinned23)) %>%
    filter(yn == 2) %>%
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned23 == 0, tnNext, C))
  
  loglike23app <- function(data,p23delay,params){
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
             absProb = ifelse(!is.na(tnAdj1),absProbVec(X,tnNext,delayReport),0)) %>%
      summarise(sum( (occTermW1-expoTermW1)*absProb+(occTermW0-expoTermW0)*(1-absProb) )) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta23app <- optim(par = c(-0.1), fn = loglike23app, data=dataloglike23, p23delay=p23delay, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  #results
  paramRes <- c(p23delay[1],
                p23delay[2],
                p23delay[3],
                theta23app)
  
  return(paramRes)
  
}


runSim <- function(){
  totalSim <- 400
  startSim <- 1
  dataParam <- data.frame(f4=NA,
                          f5=NA,
                          f6=NA,
                          theta7app=NA)
  for(i in startSim:totalSim){
    resulti <- rep(NA,ncol(dataParam))
    tryCatch({resulti <- runSimSeed(i)}, error = function(e){})
    dataParam[i,] <- resulti
    print(i)
    
  }
  
  return(dataParam)
}

dataParam <- runSim()
save(dataParam, file = "Results/dfParamEstC&K.Rda") 
