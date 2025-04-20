library(tidyverse)
library(scales)
library(plotrix)

LECDK19 <- readRDS("LECDK19.Rdata")

allIds <- unique(c(LECDK19$adjReac$id,LECDK19$adjDisab$id,LECDK19$delay$id,LECDK19$disab$id,LECDK19$reac$id))

rk4 <- function(df, a, b, f0, n) {
  
  h = (b-a)/n
  
  f = f0
  for (i in 1:n) {
    s = a + h * (i-1)
    k1 = h * df(s,f)
    k2 = h * df(s+0.5*h, f+0.5*k1)
    k3 = h * df(s+0.5*h, f+0.5*k2)
    k4 = h * df(s+h, f+k3)
    
    f = f+1/6*(k1+2*k2+2*k3+k4)
  }
  return(f)
}

diffAdjReac = function(t, x, mu12, mu21, mu13, mu14, mu24) {
  d1 = -x[1]*(mu12(t)+mu13(t)+mu14(t))+x[2]*mu21(t)
  d2 = -x[2]*(mu21(t)+mu24(t))+x[1]*mu12(t)
  d3 = x[1]*mu13(t)
  d4= x[1]*mu14(5)+x[2]*mu24(t)
  return( c(d1, d2, d3, d4) )
}

#note: model space is extended so "rejectedBefore" is captured in statespace - it is now state 4 and dead is state 5
diffAdjDisab = function(t, x, mu13, mu12, mu24, mu42, mu43, mu15, mu45, mu25) {
  d1 = -x[1]*(mu13(t)+mu12(t)+mu15(t))
  d2 = -x[2]*(mu24(t)+mu25(t))+x[4]*mu42(t)
  d3 = x[1]*mu13(t)+x[4]*mu43(t)
  d4 = -x[4]*(mu43(t)+mu42(t)+mu45(t))+x[2]*mu24(t)
  d5 = x[4]*mu45(t)+x[1]*mu15(t)+x[2]*mu25(t)
  return( c(d1, d2, d3, d4, d5) )
}


alpha <- function(t,ageAtDisab,gender,lambda,k,beta1,beta2,beta3){
  alpha0 = k*lambda^k*t^(k-1)/(exp((lambda*t)^k)-1)
  alpha = alpha0*exp(beta1*(gender=="M")+beta2*(gender=="F")+beta3*ageAtDisab)
  return(alpha)
}


adjStateDisab <- function(adjState,rejectedBefore){c((adjState==1 & rejectedBefore==0),(adjState==2),(adjState==3),(adjState==1 & rejectedBefore==1),(adjState==5))}

loglikeRep <- function(data,params){
  lambda <- params[1]
  k <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  beta3 <- params[5]
  
  if(lambda <= 0){return(Inf)}
  if(k <= 0 | k > 2){return(Inf)}
  
  loglike <- data %>% 
    rowwise() %>%
    mutate(occTerm = log(alpha(repDelayDisab,age-repDelayDisab,gender,lambda,k,beta1,beta2,beta3)), 
           expoTerm = -exp(beta1*(gender=="M")+beta2*(gender=="F")+beta3*(age-repDelayDisab))*log((1-exp(-(lambda*repDelayDisab)^k))/(1-exp(-(lambda*repDelayMax)^k))) ) %>% 
    ungroup() %>%
    summarise(sum((occTerm-expoTerm)*weight*multiplicity)) %>% 
    as.numeric()
  
  return(loglike)
}


delayDist <- function(t,ageAtDisab,gender,lambda,k,beta1,beta2,beta3){
  return( (1-exp(-(lambda*t)^k))^(exp( beta1*(gender=="M")+beta2*(gender=="F")+beta3*ageAtDisab )) )
}

tmax=8
nstep=100
eta <- "2019-09-01"
t0 <- "2015-01-31"
stepsize <- 1/12
dayYear <- 365.25
monthYear <- 12
halfMonthInDays <- 15 

#initialize parameter lists
paramList <- list()

est <- function(data, plot=FALSE, naive=FALSE){

  
  #adjudication for reactivations
  dataAdjReacActive <- data$adjReac %>% 
    filter(adjState==1) %>% 
    group_by(adjState,gender,age,durReac,durDisab) %>%
    summarise(expo = sum(expo*multiplicity), 
              awardOcc = sum(awardOcc*multiplicity), 
              rejectOcc=sum(rejectOcc*multiplicity), 
              deadOcc=sum(deadOcc*multiplicity)) %>% 
    filter(expo > 0)
  
  dataAdjReacInactive <- data$adjReac %>% 
    filter(adjState==2) %>% 
    group_by(adjState,gender,age,durReac,durDisab) %>%
    summarise(expo = sum(expo*multiplicity), 
              reapplyOcc = sum(reapplyOcc*multiplicity), 
              deadOcc=sum(deadOcc*multiplicity)) %>% 
    filter(expo > 0)
  
  modelReacReapply <- glm(data = dataAdjReacInactive, family = poisson(link = log), reapplyOcc ~ age + gender +durReac + durDisab - 1, offset = log(expo))
  modelReacAward <- glm(data = dataAdjReacActive, family = poisson(link = log), awardOcc ~ age + gender+ durReac + durDisab  - 1, offset = log(expo))
  modelReacReject <- glm(data = dataAdjReacActive, family = poisson(link = log), rejectOcc ~ age + gender+ durReac + durDisab - 1, offset = log(expo))
  modelReacDeadIBNR <- glm(data = dataAdjReacInactive, family = poisson(link = log), deadOcc ~ age + gender+ durReac + durDisab - 1, offset = log(expo))
  modelReacDeadRBNS <- glm(data = dataAdjReacActive, family = poisson(link = log), deadOcc ~ age + gender+ durReac + durDisab - 1, offset = log(expo))
  
  absProbReac <- function(y0,age0,gender0,durDisab0,durReac0){
    mu12 <- function(t){if(modelReacReject$converged == FALSE){0}else{exp(modelReacReject$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReac0+t,durDisab0+t))}}
    mu21 <- function(t){if(modelReacReapply$converged == FALSE){0}else{exp(modelReacReapply$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReac0+t,durDisab0+t))}}
    mu13 <- function(t){if(modelReacAward$converged == FALSE){0}else{exp(modelReacAward$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReac0+t,durDisab0+t))}}
    mu14 <- function(t){if(modelReacDeadRBNS$converged == FALSE){0}else{exp(modelReacDeadRBNS$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReac0+t,durDisab0+t))}}
    mu24 <- function(t){if(modelReacDeadIBNR$converged == FALSE){0}else{exp(modelReacDeadIBNR$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReac0+t,durDisab0+t))}}
    
    distrib <- rk4(function(t, x) diffAdjReac(t, x, mu12=mu12, mu21=mu21, mu13=mu13, mu14=mu14, mu24=mu24), a=0, b=tmax, f0=y0, nstep)
    return(distrib[3])
  } 
  
  
  #adjudication for disabilities
  dataAdjDisabActive <- data$adjDisab %>% 
    filter(adjState==1) %>%
    group_by(gender,age,durDisab,durDisabReport,rejectedBefore) %>%
    summarise(expo = sum(expo*multiplicity), 
              awardOcc = sum(awardOcc*multiplicity), 
              rejectOcc = sum(rejectOcc*multiplicity), 
              deadOcc=sum(deadOcc*multiplicity)) %>% 
    filter(expo > 0)
  
  dataAdjDisabInactive <- data$adjDisab  %>% 
    filter(adjState==2) %>%
    group_by(gender,age,durDisab,durDisabReport,rejectedBefore) %>%
    summarise(expo = sum(expo*multiplicity), 
              reapplyOcc = sum(reapplyOcc*multiplicity), 
              deadOcc=sum(deadOcc*multiplicity)) %>% 
    filter(expo > 0)
  
  modelDisabReapply <- glm(data = dataAdjDisabInactive, family = poisson(link = log), reapplyOcc ~ age+ gender +durDisabReport + durDisab  - 1, offset = log(expo))
  modelDisabAward <- glm(data = dataAdjDisabActive, family = poisson(link = log), awardOcc ~ age + gender + durDisabReport  + durDisab   + rejectedBefore  - 1, offset = log(expo))
  modelDisabRejec <- glm(data = dataAdjDisabActive, family = poisson(link = log), rejectOcc ~ age + gender + durDisabReport + durDisab + rejectedBefore  - 1, offset = log(expo))
  modelDeadRBNS <- glm(data = dataAdjDisabActive, family = poisson(link = log), deadOcc ~ age + gender + durDisabReport + durDisab + rejectedBefore  - 1, offset = log(expo))
  modelDeadIBNR <- glm(data = dataAdjDisabInactive, family = poisson(link = log), deadOcc ~ age + gender + durDisabReport + durDisab  - 1, offset = log(expo))
  
  absProbDisab <- function(y0,age0,gender0,durDisab0,durReport0){
    mu13 <- function(t){if(modelDisabAward$converged == FALSE){0}else{exp(modelDisabAward$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,0))}}
    mu12 <- function(t){if(modelDisabRejec$converged == FALSE){0}else{exp(modelDisabRejec$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,0))}}
    mu24 <- function(t){if(modelDisabReapply$converged == FALSE){0}else{exp(modelDisabReapply$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t))}}
    mu42 <- function(t){if(modelDisabRejec$converged == FALSE){0}else{exp(modelDisabRejec$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,1))}}
    mu43 <- function(t){if(modelDisabAward$converged == FALSE){0}else{exp(modelDisabAward$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,1))}} 
    mu15 <- function(t){if(modelDeadRBNS$converged == FALSE){0}else{exp(modelDeadRBNS$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,0))}}
    mu45 <-function(t){if(modelDeadRBNS$converged == FALSE){0}else{exp(modelDeadRBNS$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t,1))}}
    mu25 <- function(t){if(modelDeadIBNR$converged == FALSE){0}else{exp(modelDeadIBNR$coefficients %*% c(age0+t,(gender0=="F"),(gender0=="M"),durReport0+t,durDisab0+t))}}
    
    distrib <- rk4(function(t, x) diffAdjDisab(t, x, mu13=mu13, mu12=mu12, mu24=mu24, mu42=mu42, mu43=mu43, mu15=mu15, mu45=mu45, mu25=mu25), a=0, b=tmax, f0=y0, nstep)
    return(distrib[3])
  } 
  
  
  #reporting delays
  dataDelay <- data$delay %>%
    rowwise() %>%
    mutate(weight = absProbDisab(adjStateDisab(adjState,rejectedBefore),age,gender,durDisab,durDisabReport)) %>%
    ungroup() %>%
    mutate(weight = pmax(0,pmin(weight,1)))
  
  pDelay <- optim(par = c(0.5,1,0.1,2,0.01), fn = loglikeRep, data=dataDelay, control=list(fnscale=-1))$par
  
  
  #disability
  idDisab <-  data$disab %>% filter(disabOcc==1) %>% select(id) %>% unlist()
  
  dataNoDisab <- data$disab  %>% filter(!(id %in% idDisab))
  
  dataOneDisab <- data$disab %>% filter((id %in% idDisab))
  
  dataOneDisabWeight <- dataOneDisab  %>% 
    filter(disabOcc==1) %>%
    rowwise() %>%
    mutate(weight=absProbDisab(adjStateDisab(adjState,rejectedBefore),age,gender,durDisab,durDisabReport) ) %>%
    ungroup() %>%
    mutate(weight = pmax(0,pmin(weight,1))) %>%
    select(id,weight)
  
  dataOneDisabConfirm <- left_join(dataOneDisab,dataOneDisabWeight,by=c('id'))
  
  dataOneDisabReject <- dataOneDisabConfirm
  dataOneDisabReject$weight <- 1-dataOneDisabReject$weight
  
  dataOneDisabReject <- dataOneDisabReject %>% filter(weight > 0)
  
  for(i in unique(dataOneDisabReject$id)){
    idInv <- which(dataOneDisabReject$id == i & dataOneDisabReject$disabOcc==1)
    dataOneDisabReject[idInv,"disabOcc"] <- 0
    dates <- seq(as.Date(dataOneDisabReject[idInv,"dateStart"]), as.Date(eta), by="months")
    expos <- as.numeric(diff(seq(as.Date(dataOneDisabReject[idInv,"dateStart"]), as.Date(eta), by="months")))/dayYear
    dataOneDisabReject[idInv,"expo"] <- expos[c(1)]
    dates <- dates[-c(1)]
    dates <- head(dates,-1)
    if(length(dates)==0){
      next
    }
    expos <- expos[-c(1)]
    ages <- stepsize/2+dataOneDisabReject[idInv,"age"]+(seq(from=dataOneDisabReject[idInv,"age"],length.out=length(dates))-dataOneDisabReject[idInv,"age"])*stepsize
    dfId <- dataOneDisabReject[idInv,]
    dfId[1:length(dates),] <- dfId[1,]
    dfId$dateStart <- dates
    dfId$age <- ages
    dfId$expo <- expos
    dataOneDisabReject <- rbind(dataOneDisabReject,dfId)
  }
  
  dataNoDisab$weight <- 1
  dataOneDisabAdj <- dataNoDisab %>% union(dataOneDisabConfirm) %>% union(dataOneDisabReject)
  
  dataOneDisabAdjDelay <- dataOneDisabAdj
  
  dataOneDisabAdjDelay$expo <- dataOneDisabAdjDelay$expo*delayDist(as.numeric(difftime(eta,as.Date(dataOneDisabAdjDelay$dateStart)+halfMonthInDays,units="days"))/dayYear,
                                                                   dataOneDisabAdjDelay$age,dataOneDisabAdjDelay$gender,
                                                                   pDelay[1],pDelay[2],pDelay[3],pDelay[4],pDelay[5])
  dataOneDisabAdjDelay <- dataOneDisabAdjDelay %>% filter(expo > 0)
  
  dataOneDisabAdjDelayFit <- dataOneDisabAdjDelay %>%
    mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
    group_by(gender,age,time) %>%
    summarise(disabOccW=sum(disabOcc*weight*multiplicity), expoW=sum(expo*weight*multiplicity))%>%
    filter(expoW > 0)
  
  modelDisab <- glm(data = dataOneDisabAdjDelayFit, family = poisson(link = log), disabOccW ~ age + gender + time- 1, offset = log(expoW))
  
  
  #reactivation
  idReac <-  data$reac %>% filter(reacOcc==1) %>% select(id) %>% unlist()
  
  dataNoReac <-  data$reac %>% filter(!(id %in% idReac))
  
  dataOneReac <-  data$reac %>% filter((id %in% idReac))
  
  adjStateReac <- function(adjState) c((adjState==1),(adjState==2),(adjState==3),(adjState==4))
  
  dataOneReacWeight <- dataOneReac  %>% 
    filter(reacOcc==1) %>%
    rowwise() %>%
    mutate(weight=1-absProbReac(adjStateReac(adjState),age,gender,durDisab,durReac) ) %>%
    ungroup() %>%
    mutate(weight = pmax(0,pmin(weight,1))) %>%
    select(id,weight)
  
  dataOneReacConfirm <- left_join(dataOneReac,dataOneReacWeight,by=c('id'))
  
  dataOneReacReject <- dataOneReacConfirm
  dataOneReacReject$weight <- 1-dataOneReacReject$weight
  dataOneReacReject <- dataOneReacReject %>% filter(weight > 0)
  
  for(i in unique(dataOneReacReject$id)){
    idReac <- which(dataOneReacReject$id == i & dataOneReacReject$reacOcc==1)
    dataOneReacReject[idReac,"reacOcc"] <- 0
    dates <- seq(as.Date(dataOneReacReject[idReac,"dateStart"]), as.Date(eta), by="months")
    expos <- as.numeric(diff(seq(as.Date(dataOneReacReject[idReac,"dateStart"]), as.Date(eta), by="months")))/dayYear
    dataOneReacReject[idReac,"expo"] <- expos[c(1)]
    dates <- dates[-c(1)]
    dates <- head(dates,-1)
    if(length(dates)==0){
      next
    }
    expos <- expos[-c(1)]
    ages <- stepsize/2+dataOneReacReject[idReac,"age"]+(seq(from=dataOneReacReject[idReac,"age"],length.out=length(dates))-dataOneReacReject[idReac,"age"])*stepsize
    durs <- dataOneReacReject[idReac,"durDisab"]+stepsize/2+(seq(from=dataOneReacReject[idReac,"durDisab"],length.out=length(dates))-dataOneReacReject[idReac,"durDisab"])*stepsize
    dfId <- dataOneReacReject[idReac,]
    dfId[1:length(dates),] <- dfId[1,]
    dfId$dateStart <- dates
    dfId$age <- ages
    dfId$durDisab <- durs
    dfId$expo <- expos
    dataOneReacReject <- rbind(dataOneReacReject,dfId)
  }
  
  dataNoReac$weight <- 1
  dataOneReacAdj <- dataNoReac %>% union(dataOneReacConfirm) %>% union(dataOneReacReject)
  
  dataOneReacAdjFit <- dataOneReacAdj %>% 
    mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
    group_by(gender,age,durDisab,time) %>%
    summarise(occReacW=sum(reacOcc*weight*multiplicity), expoW=sum(expo*weight*multiplicity))%>%
    filter(expoW > 0)
  
  modelReac <- glm(data = dataOneReacAdjFit, family = poisson(link = log), occReacW ~ durDisab + age + gender + time - 1, offset = log(expoW))
  
  
  paramList$reacReapplyCoef <- if(modelReacReapply$converged==0){rep(NA,5)}else{modelReacReapply$coefficients}
  paramList$reacAwardCoef <- if(modelReacAward$converged==0){rep(NA,5)}else{modelReacAward$coefficients} 
  paramList$reacRejectCoef <-  if(modelReacReject$converged==0){rep(NA,5)}else{modelReacReject$coefficients}  
  paramList$reacDeadIBNRCoef <- if(modelReacDeadIBNR$converged==0){rep(NA,5)}else{modelReacDeadIBNR$coefficients}  
  paramList$reacDeadRBNSCoef <- if(modelReacDeadRBNS$converged==0){rep(NA,5)}else{modelReacDeadRBNS$coefficients}  
  
  paramList$disabReapplyCoef <- if(modelDisabReapply$converged==0){rep(NA,5)}else{modelDisabReapply$coefficients} 
  paramList$disabAwardCoef <- if(modelDisabAward$converged==0){rep(NA,6)}else{modelDisabAward$coefficients}  
  paramList$disabRejectCoef <- if(modelDisabRejec$converged==0){rep(NA,6)}else{modelDisabRejec$coefficients}  
  paramList$disabDeadIBNRCoef <- if(modelDeadIBNR$converged==0){rep(NA,5)}else{modelDeadIBNR$coefficients}   
  paramList$disabDeadRBNSCoef <- if(modelDeadRBNS$converged==0){rep(NA,6)}else{modelDeadRBNS$coefficients}  
  
  paramList$delayCoef <- pDelay
  
  paramList$disabCoef <- if(modelDisab$converged==0){rep(NA,4)}else{modelDisab$coefficients}  
  paramList$reacCoef <- if(modelReac$converged==0){rep(NA,5)}else{modelReac$coefficients}  
  
  if(plot){
    
    discrete <- 0.1
    
    #disability naive
    dataDisabNaive <- data$disab %>%     
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
      group_by(gender,age,time) %>%
      summarise(disabOcc=sum(disabOcc*multiplicity), expo=sum(expo*multiplicity))%>%
      filter(expo > 0)
    
    modelDisabNaive <- glm(data = dataDisabNaive, family = poisson(link = log), disabOcc ~ age + gender + time - 1, offset = log(expo))
    
    #disability backcensor 1 year
    dataDisabBackcens <- data$disab %>%
      filter(dateStart <= "2018-08-01 UTC") %>%
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
      group_by(gender,age,time) %>%
      summarise(disabOcc=sum(disabOcc*multiplicity), expo=sum(expo*multiplicity))%>%
      filter(expo > 0)
    
    modelDisabBackcens <- glm(data = dataDisabBackcens, family = poisson(link = log), disabOcc ~ age + gender + time - 1, offset = log(expo))
    
    #reactivation naive
    dataReacNaive <- data$reac %>% 
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
      group_by(gender,age,durDisab,time) %>%
      summarise(occReac=sum(reacOcc*multiplicity), expo=sum(expo*multiplicity))%>%
      filter(expo > 0)
    
    modelReacNaive <- glm(data = dataReacNaive, family = poisson(link = log), occReac ~ durDisab + age + gender + time - 1, offset = log(expo))
    
    #reactivation backcensor 1 year
    dataReacBackcens <- data$reac %>% 
      filter(dateStart <= "2018-08-01 UTC") %>%
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
      group_by(gender,age,durDisab,time) %>%
      summarise(occReac=sum(reacOcc*multiplicity), expo=sum(expo*multiplicity))%>%
      filter(expo > 0)
    
    modelReacBackcens <- glm(data = dataReacBackcens, family = poisson(link = log), occReac ~ durDisab + age + gender + time - 1, offset = log(expo))
    
    
    #disability
    dataDisabOE <- dataOneDisabAdjDelayFit
    dataDisabOE$pred <- predict(modelDisab,newdata=dataDisabOE, type="response")
    dataDisabOESum <- dataDisabOE %>% 
      mutate(time = round(time/discrete)*discrete) %>% 
      group_by(time) %>%
      summarise(disabOccW = sum(disabOccW), expoW=sum(expoW), pred=sum(pred)) %>%
      mutate(OE = disabOccW/expoW, OEpred = pred/expoW)
    
    dataDisabOENaive <- dataDisabNaive
    dataDisabOENaive$pred <- predict(modelDisabNaive,newdata=dataDisabOENaive, type="response")
    dataDisabOENaive$predBackcens <- predict(modelDisabBackcens,newdata=dataDisabOENaive, type="response")
    dataDisabOENaiveSum <- dataDisabOENaive %>% 
      mutate(time = round(time/discrete)*discrete) %>% 
      group_by(time) %>%
      summarise(disabOcc = sum(disabOcc), expo=sum(expo), pred=sum(pred), predBackcens=sum(predBackcens)) %>%
      mutate(OE = disabOcc/expo, OEpred = pred/expo, OEpredBackcens = predBackcens/expo)
    
    
    #reactivation
    dataReacOE <- dataOneReacAdj %>%
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear) %>%
      mutate(reacOccW=reacOcc*weight, expoW=expo*weight)
    dataReacOE$pred <- predict(modelReac, newdata=dataReacOE, type="response")
    dataReacOESum  <- dataReacOE %>% 
      mutate(time = round(time/discrete)*discrete) %>% 
      group_by(time) %>%
      summarise(occReacW=sum(reacOccW), expoW=sum(expoW), pred=sum(pred))%>%
      filter(expoW > 0) %>%
      mutate(OE = occReacW/expoW, OEpred = pred/expoW)
    
    dataReacOENaive <- data$reac  %>%
      mutate(time=as.numeric(difftime(as.Date(dateStart)+halfMonthInDays,t0, units="days"))/dayYear)
    dataReacOENaive$pred <- predict(modelReacNaive,newdata=dataReacOENaive, type="response")
    dataReacOENaive$predBackcens <- predict(modelReacBackcens,newdata=dataReacOENaive, type="response")
    dataReacNaiveOESum <- dataReacOENaive%>% 
      mutate(time = round(time/discrete)*discrete) %>% 
      group_by(time) %>%
      summarise(occReac=sum(reacOcc), expo=sum(expo), pred=sum(pred), predBackcens=sum(predBackcens))%>%
      filter(expo > 0) %>%
      mutate(OE = occReac/expo, OEpred = pred/expo, OEpredBackcens=predBackcens/expo)
    
    png("Figures/FittedRates.png", width = 10, height = 4, units = 'in', res = 300)
    par(mfrow=c(1,2))
    par(mar=c(4,4,2,2))
    
    idxNotCens <- which(dataDisabOENaiveSum$time <= max(dataDisabOENaiveSum$time)-1)
    idxCens <- which(dataDisabOENaiveSum$time >= max(dataDisabOENaiveSum$time)-1)
    
    plot(dataDisabOESum$time,dataDisabOESum$OE, 
         ylim = c(0,0.0051),
         xlab= "time (years)",
         ylab = "disability rate")
    points(dataDisabOENaiveSum$time,dataDisabOENaiveSum$OE, col=scales::alpha("black", 0.3))
    lines(dataDisabOESum$time,dataDisabOESum$OEpred, col=scales::alpha("black", 1))
    lines(dataDisabOENaiveSum$time,dataDisabOENaiveSum$OEpred, col=scales::alpha("black", 0.3))
    lines(dataDisabOENaiveSum$time[idxNotCens],dataDisabOENaiveSum$OEpredBackcens[idxNotCens], col=scales::alpha("black", 0.6))
    lines(dataDisabOENaiveSum$time[idxCens],dataDisabOENaiveSum$OEpredBackcens[idxCens], col=scales::alpha("black", 0.6), lty=2)
    
    idxNotCensReac <- which(dataReacNaiveOESum$time <= max(dataReacNaiveOESum$time)-1)
    idxCensReac <- which(dataReacNaiveOESum$time >= max(dataReacNaiveOESum$time)-1)
    
    plot(dataReacOESum$time,dataReacOESum$OE,
         xlab= "time (years)",
         ylab = "reactivation rate")
    points(dataReacNaiveOESum$time,dataReacNaiveOESum$OE, col=scales::alpha("black", 0.3))
    lines(dataReacOESum$time, dataReacOESum$OEpred)
    lines(dataReacNaiveOESum$time,dataReacNaiveOESum$OEpred, col=scales::alpha("black", 0.3))
    lines(dataReacNaiveOESum$time[idxNotCensReac],dataReacNaiveOESum$OEpredBackcens[idxNotCensReac], col=scales::alpha("black", 0.6))
    lines(dataReacNaiveOESum$time[idxCensReac],dataReacNaiveOESum$OEpredBackcens[idxCensReac], col=scales::alpha("black", 0.6), lty=2)
    axis.break(axis = 2, breakpos = 0.09, style="slash")
    
    dev.off()
    
    
    #delay
    dataDelayOENames <- c("time","age","gender","expo","occ")
    dataDelayOE <- data.frame(matrix(nrow = 0, ncol = length(dataDelayOENames)))
    colnames(dataDelayOE) <- dataDelayOENames
    
    for(i in 1:nrow(dataDelay)){
      dataDelayi <- dataDelay[i,]
      
      times <- seq(from=dataDelayi$durDisabReport,to=dataDelayi$repDelayMax,by=stepsize/2)
      ages <- dataDelayi$age - dataDelayi$durDisabReport
      genders = rep(dataDelayi$gender,times=length(times))
      expos = rep(stepsize/2,times=length(times))
      occs = c(1,rep(0,length(times)-1))
      absProbs <- rep(dataDelayi$weight,times=length(times))
      
      dataDelayOEi <- data.frame(matrix(nrow = length(times), ncol = length(dataDelayOENames)))
      colnames(dataDelayOEi) <- dataDelayOENames
      dataDelayOEi$time <- times
      dataDelayOEi$age <- ages
      dataDelayOEi$gender <- genders
      dataDelayOEi$expo <- expos*dataDelayi$weight
      dataDelayOEi$occ <- occs*dataDelayi$weight
      
      
      dataDelayOE <- rbind(dataDelayOE,dataDelayOEi)
    }
    
    dataDelayOEPred <- dataDelayOE %>%
      rowwise() %>%
      mutate(pred=alpha(time,age,gender,pDelay[1],pDelay[2],pDelay[3],pDelay[4],pDelay[5])*expo) %>%
      ungroup()
    
    dataDelayOEPlot <- dataDelayOEPred %>% 
      mutate(time = round(time/discrete)*discrete) %>%
      group_by(time) %>% 
      summarise(expo=sum(expo),occ=sum(occ),pred=sum(pred))  %>% 
      mutate(OE=occ/expo,OEpred=pred/expo)
    
    png("Figures/DelayRates.png", width = 10, height = 4, units = 'in', res = 300)
    par(mar=c(4,4,2,2))
    
    plot(dataDelayOEPlot$time,dataDelayOEPlot$OE,
         ylab = "reporting delay reverse time hazard",
         xlab = "time (years)",
         xlim = c(0,1),
         ylim = c(0,30))
    lines(dataDelayOEPlot$time,dataDelayOEPlot$OEpred)
    
    dev.off()
    
    
    #adjudication reactivation and disability
    dataAdjInactiveOE <- dataAdjReacInactive
    dataAdjInactiveOE$predReapply <- predict(modelReacReapply,newdata=dataAdjInactiveOE, type="response")
    dataAdjInactiveOE$predDeadIBNR <- predict(modelReacDeadIBNR,newdata=dataAdjInactiveOE, type="response")
    dataAdjInactiveOESum <- dataAdjInactiveOE %>% 
      mutate(durReac=round(durReac/discrete)*discrete) %>%
      group_by(durReac) %>%
      summarise(reapplyOcc = sum(reapplyOcc), 
                deadOcc = sum(deadOcc),
                expo=sum(expo), 
                predReapply=sum(predReapply),
                predDeadIBNR=sum(predDeadIBNR)) %>%
      mutate(OEReapply = reapplyOcc/expo, 
             OEDeadIBNR = deadOcc/expo,
             OEpredReapply = predReapply/expo,
             OEpredDeadIBNR = predDeadIBNR/expo)
    
    dataAdjActiveOE <- dataAdjReacActive
    dataAdjActiveOE$predAward <- predict(modelReacAward,newdata=dataAdjActiveOE, type="response")
    dataAdjActiveOE$predReject <- predict(modelReacReject,newdata=dataAdjActiveOE, type="response")
    dataAdjActiveOE$predDeadRBNS <- predict(modelReacDeadRBNS,newdata=dataAdjActiveOE, type="response")
    dataAdjActiveOESum <- dataAdjActiveOE %>% 
      mutate(durReac=round(durReac/discrete)*discrete) %>%
      group_by(durReac) %>%
      summarise(awardOcc = sum(awardOcc), 
                rejectOcc = sum(rejectOcc),
                deadOcc = sum(deadOcc),
                expo=sum(expo), 
                predAward=sum(predAward),
                predReject=sum(predReject),
                predDeadRBNS=sum(predDeadRBNS)) %>%
      mutate(OEAward = awardOcc/expo, 
             OEReject = rejectOcc/expo,
             OEDead = deadOcc/expo,
             OEpredAward = predAward/expo,
             OEpredReject = predReject/expo,
             OEpredDead = predDeadRBNS/expo)
    
    dataAdjDisabInactiveOE <- dataAdjDisabInactive
    dataAdjDisabInactiveOE$predReapply <- predict(modelDisabReapply,newdata=dataAdjDisabInactiveOE, type="response")
    dataAdjDisabInactiveOESum <- dataAdjDisabInactiveOE %>% 
      mutate(durDisab = round(durDisab/discrete)*discrete) %>%
      group_by(durDisab) %>%
      summarise(reapplyOcc = sum(reapplyOcc), 
                expo=sum(expo), 
                predReapply=sum(predReapply)) %>%
      mutate(OEReapply = reapplyOcc/expo, 
             OEpredReapply = predReapply/expo)
    
    
    dataAdjDisabActiveOE <- dataAdjDisabActive
    dataAdjDisabActiveOE$predAward <- predict(modelDisabAward,newdata=dataAdjDisabActiveOE, type="response")
    dataAdjDisabActiveOE$predReject <- predict(modelDisabRejec,newdata=dataAdjDisabActiveOE, type="response")
    dataAdjDisabActiveOE$predDeadRBNS <- predict(modelDeadRBNS,newdata=dataAdjDisabActiveOE, type="response")
    dataAdjDisabActiveOESum <- dataAdjDisabActiveOE %>% 
      mutate(durDisab = round(durDisab/discrete)*discrete) %>%
      group_by(durDisab) %>%
      summarise(awardOcc = sum(awardOcc), 
                rejectOcc = sum(rejectOcc),
                deadOcc = sum(deadOcc),
                expo=sum(expo), 
                predAward = sum(predAward), 
                predReject = sum(predReject),
                predDead = sum(predDeadRBNS)) %>%
      mutate(OEAward = awardOcc/expo, 
             OEReject = rejectOcc/expo,
             OEDead = deadOcc/expo,
             OEpredAward = predAward/expo,
             OEpredReject = predReject/expo,
             OEpredDead = predDead/expo)
    
    
    
    png("Figures/AdjRates.png", width = 6, height = 4, units = 'in', res = 300)
    par(mfrow=c(5,2))
    par(mar = c(1.6, 1.6, 0.1, 1))
    par(cex.lab=0.9, cex.axis=0.9)
    par(mgp=c(3, .5, 0))
    
    plot(dataAdjInactiveOESum$durReac, dataAdjInactiveOESum$OEReapply)
    legend("topright", legend=expression(paste('ω'[21],' reactivation')), bty = "n")
    lines(dataAdjInactiveOESum$durReac, dataAdjInactiveOESum$OEpredReapply)
    
    plot(dataAdjInactiveOESum$durReac, dataAdjInactiveOESum$OEDeadIBNR)
    legend("topright", legend=expression(paste('ω'[24],' reactivation')), bty = "n")
    lines(dataAdjInactiveOESum$durReac, dataAdjInactiveOESum$OEpredDeadIBNR)
    
    plot(dataAdjActiveOESum$durReac, dataAdjActiveOESum$OEAward)
    legend("topright", legend=expression(paste('ω'[13],' reactivation')), bty = "n")
    lines(dataAdjActiveOESum$durReac, dataAdjActiveOESum$OEpredAward)
    
    plot(dataAdjActiveOESum$durReac, dataAdjActiveOESum$OEReject)
    legend("topright", legend=expression(paste('ω'[12],' reactivation')), bty = "n")
    lines(dataAdjActiveOESum$durReac, dataAdjActiveOESum$OEpredReject)
    
    plot(dataAdjActiveOESum$durReac, rep(0,nrow(dataAdjActiveOESum)))
    legend("topright", legend=expression(paste('ω'[14],' reactivation')), bty = "n")
    lines(dataAdjActiveOESum$durReac, rep(0,nrow(dataAdjActiveOESum)))
    
    plot(dataAdjDisabInactiveOESum$durDisab, dataAdjDisabInactiveOESum$OEReapply)
    legend("topright", legend=expression(paste('ω'[21],' disability')), bty = "n")
    lines(dataAdjDisabInactiveOESum$durDisab, dataAdjDisabInactiveOESum$OEpredReapply)
    
    plot(dataAdjDisabInactiveOESum$durDisab, rep(0,nrow(dataAdjDisabInactiveOESum)))
    legend("topright", legend=expression(paste('ω'[24],' disability')), bty = "n")
    lines(dataAdjDisabInactiveOESum$durDisab, rep(0,nrow(dataAdjDisabInactiveOESum)))
    
    plot(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEAward)
    legend("topright", legend=expression(paste('ω'[13],' disability')), bty = "n")
    lines(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEpredAward)
    
    plot(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEReject)
    legend("topright", legend=expression(paste('ω'[12],' disability')), bty = "n")
    lines(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEpredReject)
    
    plot(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEDead)
    legend("topright", legend=expression(paste('ω'[14],' disability')), bty = "n")
    lines(dataAdjDisabActiveOESum$durDisab, dataAdjDisabActiveOESum$OEpredDead)
    
    dev.off()
    
  }
  
  if(naive){
    listParamAndNaive <- list()
    listParamAndNaive$paramList <- paramList
    listParamAndNaive$modelDisabNaive <- modelDisabNaive$coefficients
    listParamAndNaive$modelReacNaive <- modelReacNaive$coefficients
    listParamAndNaive$modelDisabBackcens <- modelDisabBackcens$coefficients
    listParamAndNaive$modelReacBackcens <- modelReacBackcens$coefficients
    return(listParamAndNaive)
  }
  
  return(paramList)
}



runSim <- function(){
  
  startSim <- 1
  totalSim <- 500
  
  paramNames <- c("reacReapplyAge","reacReapplyGenderF","reacReapplyGenderM","reacReapplyDurReac","reacReapplyDurDisab",
                  "reacAwardAge","reacAwardGenderF","reacAwardGenderM","reacAwardDurReac","reacAwardDurDisab",
                  "reacRejectAge","reacRejectGenderF","reacRejectGenderM","reacRejectDurReac","reacRejectDurDisab",
                  "reacDeadIBNRAge","reacDeadIBNRGenderF","reacDeadIBNRGenderM","reacDeadIBNRDurReac","reacDeadIBNRDurDisab",
                  "reacDeadRBNSAge","reacDeadRBNSGenderF","reacDeadRBNSGenderM","reacDeadRBNSDurReac","reacDeadRBNSDurDisab",
                  "disabReapplyAge","disabReapplyGenderF","disabReapplyGenderM","disabReapplyDurDisabReport","disabReapplyDurDisab",
                  "disabAwardAge","disabAwardGenderF","disabAwardGenderM","disabAwardDurDisabReport","disabAwardDurDisab","disabAwardRejectedBefore",
                  "disabRejectAge","disabRejectGenderF","disabRejectGenderM","disabRejectDurDisabReport","disabRejectDurDisab","disabRejectRejectedBefore",
                  "disabDeadIBNRAge","disabDeadIBNRGenderF","disabDeadIBNRGenderM","disabDeadIBNRDurDisabReport","disabDeadIBNRDurDisab",
                  "disabDeadRBNSAge","disabDeadRBNSGenderF","disabDeadRBNSGenderM","disabDeadRBNSDurDisabReport","disabDeadRBNSDurDisab","disabDeadRBNSRejectedBefore",
                  "delayLambda","delayK","delayGenderM","delayGenderF","delayAge",
                  "disabAge","disabGenderF","disabGenderM","disabTime",
                  "reacDurDisab","reacAge","reacGenderF","reacGenderM","reacTime"
  )
  paramBoot <- as.data.frame(matrix(ncol=length(paramNames)))
  colnames(paramBoot) <- paramNames
  
  dataOrig <- LECDK19
  
  for(m in startSim:totalSim){
    set.seed(m)
    resultm <- NA
    
    dataEst <- list()
    dataEst$adjReac <- dataOrig$adjReac[0,]
    dataEst$adjDisab <- dataOrig$adjDisab[0,]
    dataEst$delay <- dataOrig$delay[0,]
    dataEst$disab <- dataOrig$disab[0,]
    dataEst$reac <- dataOrig$reac[0,]
    
    ids <- sample(allIds,size=length(allIds),replace=TRUE)
    idCount <- as.data.frame(ids) %>% group_by(ids) %>% summarise(multiplicity=n())
    
    dataEst$adjReac <- dataOrig$adjReac %>% inner_join(idCount, by=c("id"="ids"))
    dataEst$adjDisab <- dataOrig$adjDisab %>% inner_join(idCount, by=c("id"="ids"))
    dataEst$delay <- dataOrig$delay %>% inner_join(idCount, by=c("id"="ids"))
    dataEst$disab <- dataOrig$disab %>% inner_join(idCount, by=c("id"="ids"))
    dataEst$reac <- dataOrig$reac %>% inner_join(idCount, by=c("id"="ids"))
    
    tryCatch({resultm <- est(dataEst)}, error = function(e){})
    
    tryCatch({paramBoot[m,1:5] <- resultm$reacReapplyCoef}, error = function(e){})
    tryCatch({paramBoot[m,6:10] <- resultm$reacAwardCoef}, error = function(e){})
    tryCatch({paramBoot[m,11:15] <- resultm$reacRejectCoef}, error = function(e){})
    tryCatch({paramBoot[m,16:20] <- resultm$reacDeadIBNRCoef}, error = function(e){})
    tryCatch({paramBoot[m,21:25] <- resultm$reacDeadRBNSCoef}, error = function(e){})
    tryCatch({paramBoot[m,26:30] <- resultm$disabReapplyCoef}, error = function(e){})
    tryCatch({paramBoot[m,31:36] <- resultm$disabAwardCoef}, error = function(e){})
    tryCatch({paramBoot[m,37:42] <- resultm$disabRejectCoef}, error = function(e){})
    tryCatch({paramBoot[m,43:47] <- resultm$disabDeadIBNRCoef}, error = function(e){})
    tryCatch({paramBoot[m,48:53] <- resultm$disabDeadRBNSCoef}, error = function(e){})
    tryCatch({paramBoot[m,54:58] <- resultm$delayCoef}, error = function(e){})
    tryCatch({paramBoot[m,59:62] <- resultm$disabCoef}, error = function(e){})
    tryCatch({paramBoot[m,63:67] <- resultm$reacCoef}, error = function(e){})
    
    ##To save partial results, uncomment here:
    #if( (m <= startSim+9) | (m %% 10 == 0) ){
    #  save(paramBoot, file = paste0("Results/paramBoot",m,".Rda")) 
    #}
    
    gc()
  }
  
  return(paramBoot)
}

##Ordinary estimation
LECDK19Mult <- LECDK19
LECDK19Mult$adjReac <- LECDK19Mult$adjReac %>% mutate(multiplicity=1)
LECDK19Mult$adjDisab <- LECDK19Mult$adjDisab %>% mutate(multiplicity=1)
LECDK19Mult$delay <- LECDK19Mult$delay %>% mutate(multiplicity=1)
LECDK19Mult$disab <- LECDK19Mult$disab %>% mutate(multiplicity=1)
LECDK19Mult$reac <- LECDK19Mult$reac %>% mutate(multiplicity=1)

paramEst <- est(LECDK19Mult,plot=TRUE,naive=TRUE)
#save(paramEst, file = "Results/paramEst.Rda") 

##Bootstrap estimation
paramBoot <- runSim()

paramBoot <- paramBoot %>% filter(!(is.na(disabAge))) %>% head(n=400)
#save(paramBoot, file = "Results/paramBoot.Rda") 


q95 <- (1-0.95)/2

#disability
quantile(paramBoot$disabAge,q95, na.rm=TRUE)
quantile(paramBoot$disabAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabTime,q95, na.rm=TRUE)
quantile(paramBoot$disabTime,1-q95, na.rm=TRUE)

#reactivation
quantile(paramBoot$reacAge,q95, na.rm=TRUE)
quantile(paramBoot$reacAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacTime,q95, na.rm=TRUE)
quantile(paramBoot$reacTime,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacDurDisab,1-q95, na.rm=TRUE)

#reactivation adjudication

quantile(paramBoot$reacRejectAge,q95, na.rm=TRUE)
quantile(paramBoot$reacRejectAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacRejectGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacRejectGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacRejectGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacRejectGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacRejectDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacRejectDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$reacRejectDurReac,q95, na.rm=TRUE)
quantile(paramBoot$reacRejectDurReac,1-q95, na.rm=TRUE)


quantile(paramBoot$reacAwardAge,q95, na.rm=TRUE)
quantile(paramBoot$reacAwardAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacAwardGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacAwardGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacAwardGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacAwardGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacAwardDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacAwardDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$reacAwardDurReac,q95, na.rm=TRUE)
quantile(paramBoot$reacAwardDurReac,1-q95, na.rm=TRUE)


quantile(paramBoot$reacDeadRBNSAge,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadRBNSAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadRBNSGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadRBNSGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadRBNSGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadRBNSGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadRBNSDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadRBNSDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadRBNSDurReac,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadRBNSDurReac,1-q95, na.rm=TRUE)


quantile(paramBoot$reacReapplyAge,q95, na.rm=TRUE)
quantile(paramBoot$reacReapplyAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacReapplyGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacReapplyGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacReapplyGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacReapplyGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacReapplyDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacReapplyDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$reacReapplyDurReac,q95, na.rm=TRUE)
quantile(paramBoot$reacReapplyDurReac,1-q95, na.rm=TRUE)


quantile(paramBoot$reacDeadIBNRAge,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadIBNRAge,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadIBNRGenderM,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadIBNRGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadIBNRGenderF,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadIBNRGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadIBNRDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadIBNRDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$reacDeadIBNRDurReac,q95, na.rm=TRUE)
quantile(paramBoot$reacDeadIBNRDurReac,1-q95, na.rm=TRUE)

#disability adjudication

quantile(paramBoot$disabRejectAge,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabRejectGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabRejectGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabRejectDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$disabRejectDurDisabReport,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectDurDisabReport,1-q95, na.rm=TRUE)

quantile(paramBoot$disabRejectRejectedBefore,q95, na.rm=TRUE)
quantile(paramBoot$disabRejectRejectedBefore,1-q95, na.rm=TRUE)


quantile(paramBoot$disabAwardAge,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabAwardGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabAwardGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabAwardDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$disabAwardDurDisabReport,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardDurDisabReport,1-q95, na.rm=TRUE)

quantile(paramBoot$disabAwardRejectedBefore,q95, na.rm=TRUE)
quantile(paramBoot$disabAwardRejectedBefore,1-q95, na.rm=TRUE)


quantile(paramBoot$disabDeadRBNSAge,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadRBNSGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadRBNSGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadRBNSDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadRBNSDurDisabReport,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSDurDisabReport,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadRBNSRejectedBefore,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadRBNSRejectedBefore,1-q95, na.rm=TRUE)


quantile(paramBoot$disabReapplyAge,q95, na.rm=TRUE)
quantile(paramBoot$disabReapplyAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabReapplyGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabReapplyGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabReapplyGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabReapplyGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabReapplyDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$disabReapplyDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$disabReapplyDurDisabReport,q95, na.rm=TRUE)
quantile(paramBoot$disabReapplyDurDisabReport,1-q95, na.rm=TRUE)


quantile(paramBoot$disabDeadIBNRAge,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadIBNRAge,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadIBNRGenderM,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadIBNRGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadIBNRGenderF,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadIBNRGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadIBNRDurDisab,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadIBNRDurDisab,1-q95, na.rm=TRUE)

quantile(paramBoot$disabDeadIBNRDurDisabReport,q95, na.rm=TRUE)
quantile(paramBoot$disabDeadIBNRDurDisabReport,1-q95, na.rm=TRUE)

#delay

quantile(paramBoot$delayLambda,q95, na.rm=TRUE)
quantile(paramBoot$delayLambda,1-q95, na.rm=TRUE)

quantile(paramBoot$delayK,q95, na.rm=TRUE)
quantile(paramBoot$delayK,1-q95, na.rm=TRUE)

quantile(paramBoot$delayGenderM,q95, na.rm=TRUE)
quantile(paramBoot$delayGenderM,1-q95, na.rm=TRUE)

quantile(paramBoot$delayGenderF,q95, na.rm=TRUE)
quantile(paramBoot$delayGenderF,1-q95, na.rm=TRUE)

quantile(paramBoot$delayAge,q95, na.rm=TRUE)
quantile(paramBoot$delayAge,1-q95, na.rm=TRUE)
