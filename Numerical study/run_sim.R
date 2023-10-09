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
  p1adj <- NA
  p2adj <- NA
  p13delay <- NA
  p23delay <- NA
  theta1 <- NA
  theta2 <- NA
  theta12app <- NA
  theta13app <- NA
  theta23app <- NA
  theta12oracle <- NA
  theta13oracle <- NA
  theta23oracle <- NA
  theta12backcensor0 <- NA
  theta13backcensor0 <- NA
  theta23backcensor0 <- NA
  theta12backcensor1 <- NA
  theta13backcensor1 <- NA
  theta23backcensor1 <- NA
  
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
  
  tryCatch({p1adj <- optim(par = c(1), fn = loglikeAdj12, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
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
  
  tryCatch({p2adj <- optim(par = c(-1), fn = loglikeAdj23, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  ##Estimate reporting delay
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
  
  tryCatch({p13delay <- optim(par = c(1,1,1), fn = loglikeRep13, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
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
  
  tryCatch({p23delay <- optim(par = c(1,1,1), fn = loglikeRep23, data=data, p1adj=p1adj, p2adj=p2adj, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  ##Estimate hazard
  
  #model based
  
  #Poisson approximation
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
  
  tryCatch({theta12app <- optim(par = c(log(0.15),0.1,0.3), fn = loglike12app, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  delayDist13 <- function(t,X,lambda,k,beta) (1-exp(-(lambda*(tstop-t))^k))^(exp(X*beta))
  
  nu13 <- function(t,X,p1,p2,p3,lambda,k,beta){
    delayDist <- delayDist13(t,X,lambda,k,beta)
    return(muTilde13(t,X,p1,p2,p3)*delayDist)
  }
  
  
  dataloglike13 <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned13 = sum(yn==1 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned13 = ifelse(is.na(anyThinned13),0,anyThinned13)) %>%
    filter(yn == 1) %>% 
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned13 == 0, tnNext, C))
  
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
  
  tryCatch({theta13app <- optim(par = c(log(0.1),0.03,-0.1), fn = loglike13app, data=dataloglike13, p13delay=p13delay, control=list(fnscale=-1))$par}, error = function(e){})
  
  
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
  
  tryCatch({theta23app <- optim(par = c(-0.1), fn = loglike23app, data=dataloglike23, p23delay=p23delay, p1adj=p1adj, p2adj=p2adj, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  
  #Exact
  erfiApp <- function(z) pi^(-1/2)*(2*z+2/3*z^3+1/5*z^5+1/21*z^7+1/108*z^9+1/660*z^11+1/4680*z^13+1/37800*z^15) #Taylor approx around 0
  
  P11 <- function(t,s,X,p1,p2,p3,p4,p5,p6, type="approx"){
    if(type=="approx"){
      integral <- exp(p4+p6*cos(0.5*pi*X))/sqrt(p5)*sqrt(pi)/2*(erfiApp(sqrt(p5)*t)-erfiApp(sqrt(p5)*s))+exp(p1+p2*X+p3*sin(0.5*pi*X))/p2*(exp(p2*t)-exp(p2*s))
    }
    else{
      integral <- exp(p4+p6*cos(0.5*pi*X))/sqrt(p5)*sqrt(pi)/2*(erfi(sqrt(p5)*t)-erfi(sqrt(p5)*s))+exp(p1+p2*X+p3*sin(0.5*pi*X))/p2*(exp(p2*t)-exp(p2*s))
    }
    return(exp(-integral))
  }
  
  denomIntegrand1 <- function(t,s,X,p13delay,p1,p2,p3,p4,p5,p6){
    return(P11(t,s,X,p1,p2,p3,p4,p5,p6)*(muTilde12(t,X,p1,p2,p3)+nu13(t,X,p4,p5,p6,p13delay[1],p13delay[2],p13delay[3])))
  }
  
  mu1Factor <- function(t,s,k,X,p13delay,p1,p2,p3,p4,p5,p6){
    numerator <- P11(t,s,X,p1,p2,p3,p4,p5,p6) 
    denominator <- 1-sapply(t, function(tt) integrate(f=denomIntegrand1, lower=s, upper=tt, s=s, X=X, p13delay=p13delay, p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6)$value) 
    
    result <- numerator/denominator 
    result <- min(result,1) #for numerical stability
    result <- max(result,numerator) #for numerical stability
    
    return(numerator/denominator)
  }
  
  mu1k <- function(t,s,k,X,p13delay,p1,p2,p3,p4,p5,p6){
    if(k == 2){
      nu <- muTilde12(t,X,p1,p2,p3)
    }
    else if(k == 3){
      nu <- nu13(t,X,p4,p5,p6,p13delay[1],p13delay[2],p13delay[3])
    }
    factor <- mu1Factor(t,s,k,X,p13delay,p1,p2,p3,p4,p5,p6)
    
    return(factor*nu)
  }
  
  mu1dot <- function(t,s,X,p13delay,p1,p2,p3,p4,p5,p6){
    return( mu1Factor(t,s,k,X,p13delay,p1,p2,p3,p4,p5,p6)*(muTilde12(t,X,p1,p2,p3)+nu13(t,X,p4,p5,p6,p13delay[1],p13delay[2],p13delay[3])) )
  }
  
  dataloglike1 <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned13 = sum(yn==1 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned13 = ifelse(is.na(anyThinned13),0,anyThinned13)) %>%
    filter(yn == 1) %>%
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned13 == 0, tnNext, C)) %>%
    mutate(occTerm2 = NA, occTerm3 = NA, expoTerm = NA)
  
  loglike1 <- function(data,p13delay,params){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    #restrict search space for numerical stability
    if(p1 > -1 | p1 < -4 ){return(-Inf)}
    if(p2 > 1 | p2 < 0 ){return(-Inf)}
    if(p3 > 1 | p3 < -1 ){return(-Inf)}
    if(p4 > -1 | p4 < -4 ){return(-Inf)}
    if(p5 > 0.06 | p5 < 0 ){return(-Inf)}
    if(p6 > 1 | p6 < -1){return(-Inf)}
    
    
    loglike <- data %>% 
      rowwise() %>%
      mutate(occTerm2 = ifelse(!is.na(ynNext) & ynNext==2, log(mu1k(tnNext,tn,2,X,p13delay,p1,p2,p3,p4,p5,p6)), 0), 
             occTerm3 = ifelse(!is.na(ynNext) & ynNext==3 & anyThinned13 == 0, log(mu1k(tnNext,tn,3,X,p13delay,p1,p2,p3,p4,p5,p6)), 0),
             expoTerm = integrate(f=mu1dot, lower=tn, upper=tend, s=tn, X=X, p13delay=p13delay, p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6)$value) %>% 
      ungroup() %>%
      summarise(sum(occTerm2+occTerm3-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta1 <- optim(par = c(theta12app,theta13app), fn = loglike1, data=dataloglike1, p13delay=p13delay, control=list(fnscale=-1), method="BFGS")$par}, error = function(e){})
  
  
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
  
  dataloglike2 <- data %>% 
    group_by(id) %>% 
    mutate(anyThinned23 = sum(yn==2 & ynNext == 3 & thinned == TRUE) ) %>%
    ungroup() %>% 
    mutate(anyThinned23 = ifelse(is.na(anyThinned23),0,anyThinned23)) %>%
    filter(yn == 2) %>%
    mutate(tend = ifelse(!is.na(ynNext) & anyThinned23 == 0, tnNext, C)) %>%
    mutate(occTerm = NA, expoTerm = NA)
  
  loglike2 <- function(data,p23delay,params){
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
  
  tryCatch({theta2 <- optim(par = c(theta23app), fn = loglike2, data=dataloglike2, p23delay=p23delay, control=list(fnscale=-1), method="BFGS")$par}, error = function(e){})
  
  
  
  #oracle
  loglike12oracle <- function(data,params){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    loglike <- data %>% 
      filter(yn == 1) %>% 
      mutate(tend = ifelse(!is.na(ynNext), tnNext, C), 
             occTerm = ifelse(!is.na(ynNext) & ynNext==2, log(muTilde12(tnNext,X,p1,p2,p3)), 0), 
             expoTerm = 1/p2*(muTilde12(tend,X,p1,p2,p3)-muTilde12(tn,X,p1,p2,p3)) ) %>% 
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta12oracle <- optim(par = c(0.1,0.1,0.1), fn = loglike12oracle, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  mu13Integral <- function(t,s,X,p1,p2,p3,type="approx"){
    if(type=="approx"){
      integral <-  sqrt(pi/(4*p2))*exp(p1+p3*cos(0.5*pi*X))*(erfiApp(sqrt(p2)*t)-erfiApp(sqrt(p2)*s))
    }
    else{
      integral <- sqrt(pi/(4*p2))*exp(p1+p3*cos(0.5*pi*X))*(erfi(sqrt(p2)*t)-erfi(sqrt(p2)*s))
    }
    return(integral)
  }
  
  loglike13oracle <- function(data,params){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    #restrict search space for numerical stability
    if(p1 > -1 | p1 < -4 ){return(-Inf)}
    if(p2 > 1 | p2 < 0 ){return(-Inf)}
    if(p3 > 0 | p3 < -1 ){return(-Inf)}
    
    loglike <- data %>% 
      filter(yn == 1) %>% 
      mutate(tend = ifelse(!is.na(ynNext), tnNext, C), 
             occTerm = ifelse(!is.na(ynNext) & ynNext==3, log(muTilde13(tnNext,X,p1,p2,p3)), 0), 
             expoTerm = mu13Integral(tend,tn,X,params[1],params[2],params[3]) ) %>% 
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta13oracle <- optim(par = c(-1.5,0.1,-0.1), fn = loglike13oracle, data=data, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  
  dataOracle <- data %>% mutate(adjState = ifelse(!is.na(tnAdj3),3,ifelse(!is.na(tnAdj2),2,1)),
                                adjDur = ifelse(adjState==3,0,ifelse(adjState==2,tstop-tnAdj2,ifelse(!is.na(tnAdj1),tstop-tnAdj1,0))),
                                absProb = absProbVec(adjState,adjDur,X,0.8,-1.2)) %>%
    rowwise() %>%
    mutate(xi= ifelse(!is.na(absProb),rbinom(1,1,absProb),NA)) %>%
    ungroup()
  
  loglike23oracle <- function(data,params){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    if(p1 > 0 | p1 < -0.6 ){return(-Inf)}
    
    loglike <- data %>% 
      filter(yn == 2) %>%
      mutate(tend = ifelse(!is.na(ynNext) & xi==1, tnNext, C)) %>%
      mutate(occTerm = ifelse(!is.na(ynNext) & xi==1, log(mu23Ast(tnNext,tn,X,p1)), 0)) %>%
      rowwise() %>%
      mutate(expoTerm = integrate(f=mu23Ast, lower=tn, upper=tend, s=tn, X=X, p1=p1)$value ) %>% 
      ungroup() %>%
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta23oracle <- optim(par = c(-0.1), fn = loglike23oracle, data=dataOracle, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  #backcensoring 0
  backcensor0 <- tstop
  
  loglike12backcensor <- function(data,params,backcensor){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    loglike <- data %>% 
      group_by(id) %>% 
      mutate(anyThinned13 = sum(yn==1 & ynNext == 3 & thinned == TRUE) ) %>%
      ungroup() %>% 
      mutate(anyThinned13 = ifelse(is.na(anyThinned13),0,anyThinned13)) %>% 
      filter(yn == 1) %>% 
      mutate(Cback = pmin(C,backcensor),
             tend = ifelse(!is.na(ynNext) & tnNext <= Cback & anyThinned13==0, tnNext, Cback), 
             occTerm = ifelse(!is.na(ynNext) & ynNext==2 & tnNext <= Cback, log(muTilde12(tnNext,X,p1,p2,p3)), 0), 
             expoTerm = 1/p2*(muTilde12(tend,X,p1,p2,p3)-muTilde12(tn,X,p1,p2,p3)) ) %>% 
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  tryCatch({theta12backcensor0 <- optim(par = c(0.1,0.1,0.1), fn = loglike12backcensor, data=data, backcensor=backcensor0, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  loglike13backcensor <- function(data,params,backcensor){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    #restrict search space for numerical stability
    if(p1 > -1 | p1 < -4 ){return(-Inf)}
    if(p2 > 1 | p2 < 0 ){return(-Inf)}
    if(p3 > 0 | p3 < -1 ){return(-Inf)}
    
    loglike <- data %>% 
      group_by(id) %>% 
      mutate(anyThinned13 = sum(yn==1 & ynNext == 3 & thinned == TRUE) ) %>%
      ungroup() %>% 
      mutate(anyThinned13 = ifelse(is.na(anyThinned13),0,anyThinned13)) %>% 
      filter(yn == 1) %>% 
      mutate(Cback = pmin(C,backcensor),
             tend = ifelse(!is.na(ynNext) & tnNext <= Cback & anyThinned13==0, tnNext, Cback), 
             occTerm = ifelse(!is.na(ynNext) & ynNext==3 & tnNext <= Cback & anyThinned13==0, log(muTilde13(tnNext,X,p1,p2,p3)), 0), 
             expoTerm = mu13Integral(tend,tn,X,params[1],params[2],params[3]) ) %>% 
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta13backcensor0 <- optim(par = c(-1.5,0.1,-0.1), fn = loglike13backcensor, data=data, backcensor=backcensor0, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  loglike23backcensor <- function(data,params,backcensor){
    for(i in 1:length(params)){assign(paste0("p",i),params[i])}
    
    if(p1 > 0 | p1 < -0.6 ){return(-Inf)}
    
    loglike <- data %>% 
      group_by(id) %>% 
      mutate(anyThinned23 = sum(yn==2 & ynNext == 3 & thinned == TRUE) ) %>%
      ungroup() %>% 
      mutate(anyThinned23 = ifelse(is.na(anyThinned23),0,anyThinned23)) %>%
      filter(yn == 2) %>%
      mutate(Cback = pmin(C,backcensor),
             TooLongAdj = ((backcensor==tstop-2) & (is.na(tnAdj3)) & (!is.na(tnAdj1)) & (tstop-tnAdj1 > 2)),
             tend = ifelse(!is.na(ynNext) & tnNext <= Cback & anyThinned23 == 0 & !TooLongAdj, tnNext, Cback),
             occTerm = ifelse(!is.na(ynNext) & tnNext <= Cback & anyThinned23 == 0 & !TooLongAdj, log(mu23Ast(tnNext,tn,X,p1)), 0)) %>%
      rowwise() %>%
      mutate(expoTerm = integrate(f=mu23Ast, lower=tn, upper=tend, s=tn, X=X, p1=p1)$value ) %>% 
      ungroup() %>%
      summarise(sum(occTerm-expoTerm)) %>%
      as.numeric()
    
    return(loglike)
  }
  
  tryCatch({theta23backcensor0 <- optim(par = c(-0.1), fn = loglike23backcensor, data=data, backcensor=backcensor0, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  #backcensoring 1
  
  backcensor1 <- tstop-1
  
  tryCatch({theta12backcensor1 <- optim(par = c(0.1,0.1,0.1), fn = loglike12backcensor, data=data, backcensor=backcensor1, control=list(fnscale=-1))$par}, error = function(e){})
  tryCatch({theta13backcensor1 <- optim(par = c(-1.5,0.1,-0.1), fn = loglike13backcensor, data=data, backcensor=backcensor1, control=list(fnscale=-1))$par}, error = function(e){})
  tryCatch({theta23backcensor1 <- optim(par = c(-0.1), fn = loglike23backcensor, data=data, backcensor=backcensor1, control=list(fnscale=-1))$par}, error = function(e){})
  
  
  #results
  paramRes <- c(p1adj,
                p2adj,
                p13delay[1],
                p13delay[2],
                p13delay[3],
                p23delay[1],
                p23delay[2],
                p23delay[3],
                theta1[1],
                theta1[2],
                theta1[3],
                theta1[4],
                theta1[5],
                theta1[6],
                theta2,
                theta12app[1],
                theta12app[2],
                theta12app[3],
                theta13app[1],
                theta13app[2],
                theta13app[3],
                theta23app,
                theta12oracle[1],
                theta12oracle[2],
                theta12oracle[3],
                theta13oracle[1],
                theta13oracle[2],
                theta13oracle[3],
                theta23oracle,
                theta12backcensor0[1],
                theta12backcensor0[2],
                theta12backcensor0[3],
                theta13backcensor0[1],
                theta13backcensor0[2],
                theta13backcensor0[3],
                theta23backcensor0,
                theta12backcensor1[1],
                theta12backcensor1[2],
                theta12backcensor1[3],
                theta13backcensor1[1],
                theta13backcensor1[2],
                theta13backcensor1[3],
                theta23backcensor1)
  
  return(paramRes)
  
}


runSim <- function(){
  totalSim <- 400
  startSim <- 1
  dataParam <- data.frame(g1=NA,
                          g2=NA,
                          f1=NA,
                          f2=NA,
                          f3=NA,
                          f4=NA,
                          f5=NA,
                          f6=NA,
                          theta1=NA,
                          theta2=NA,
                          theta3=NA,
                          theta4=NA,
                          theta5=NA,
                          theta6=NA,
                          theta7=NA,
                          theta1app=NA,
                          theta2app=NA,
                          theta3app=NA,
                          theta4app=NA,
                          theta5app=NA,
                          theta6app=NA,
                          theta7app=NA,
                          theta1oracle=NA,
                          theta2oracle=NA,
                          theta3oracle=NA,
                          theta4oracle=NA,
                          theta5oracle=NA,
                          theta6oracle=NA,
                          theta7oracle=NA,
                          theta1cens0=NA,
                          theta2cens0=NA,
                          theta3cens0=NA,
                          theta4cens0=NA,
                          theta5cens0=NA,
                          theta6cens0=NA,
                          theta7cens0=NA,
                          theta1cens1=NA,
                          theta2cens1=NA,
                          theta3cens1=NA,
                          theta4cens1=NA,
                          theta5cens1=NA,
                          theta6cens1=NA,
                          theta7cens1=NA)
  for(i in startSim:totalSim){
    resulti <- rep(NA,ncol(dataParam))
    tryCatch({resulti <- runSimSeed(i)}, error = function(e){})
    dataParam[i,] <- resulti
    
    ##To save partial results, uncomment here:
    #if( (i <= startSim+9) | (i %% 20 == 0) ){
    #  save(dataParam, file = paste0("Results/dfParamEst",i,".Rda")) 
    #}
    
  }
  
  return(dataParam)
}

dataParam <- runSim()
save(dataParam, file = "Results/dfParamEst.Rda") 