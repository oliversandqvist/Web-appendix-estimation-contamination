#Reference: "On Lewis' simulation method for point processes". For more details, "Exact simulation of the Jump Times of a Class of Piecewise Deterministic Markov Processes".

no_states <- 3

#Computes the transition rate for a given transition at a given time and duration
#x is current state #y is destination state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0 (we use age0 = age at time 0 and gender is gender where female is 0 and male is 1)
muTilde12 <- function(t,cov,p1,p2,p3){
  return(exp(p1+p2*(t+cov)+p3*sin(0.5*pi*cov)))
}

muTilde13 <- function(t,cov,p1,p2,p3){
  return(exp(p1+p2*t^2+p3*cos(0.5*pi*cov))) 
}

muTilde23 <- function(t,s,cov,p1){
  return(exp(p1*(t-s)*cov^2))
}

params <- c(log(0.15),0.1,0.4,log(0.1),0.03,-0.3,-0.3)

transition_rate <-  function(x, y, s, t, cov){
  if(x == 1 && y==2){
    as.numeric(muTilde12(t,cov,params[1],params[2],params[3]))
  }
  else if(x == 1 && y==3){
    as.numeric(muTilde13(t,cov,params[4],params[5],params[6])) 
  }
  else if(x==2 && y==3){
    as.numeric(muTilde23(t,s,cov,params[7]))
  }
  else{
    0
  }
}

#Computes the total intensity of out the current state
#x is current state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0
jump_rate <-  function(x, s, t, cov){
  sum(sapply(1:no_states, FUN = function(y) transition_rate(x,y,s,t,cov)))
}

#Computes the mark distribution given a jump time
#x is current state #s is time of last jump #t is time of current jump
mark_dist <- function(x, s, t, cov){
  rates <- sapply(1:no_states, FUN = function(y) transition_rate(x,y,s,t,cov))
  rates/sum(rates)
}

#Simulates time from the last jump to the next jump 
#x is current state #s is time of last jump #b is bound #cov is the baseline covariates known at time 0
jump_simulator <- function(x, s, tstop, b, cov){ 
  u <- runif(1)
  e <- rexp(1, rate = 1)
  t <- e/b
  if (t > tstop) { return(Inf) }
  while(u > jump_rate(x, s, s+t, cov)/b){ 
    u <- runif(1)
    e <- rexp(1, rate = 1)
    t <- e/b+t
    if (t > tstop) { return(Inf) }
  }
  return(t)
}

#Computes the "local upper bounds" of the intensities for each of the states.
#s is the time of the last jump at time 0 #t is the current time #tstop is the terminal time #cov is the baseline covariates known at time 0
rate_bounds <- function(s, t, tstop, cov){
  Xmax <- 4
  b1 <- exp(params[1]+params[2]*(tstop+Xmax)+params[3])+exp(params[4]+params[5]*tstop^2+params[6]*1)
  b2 <- 1
  b3 <- 0
  rate_bounds <- c(b1,b2,b3)
  return(rate_bounds)
}

#Simulates one trajectory of the semi-Markov process.
#x0 is initial state #s0 is last jump time #t0 is current time (note that by convention, we include the zero'th jump-time and mark in the event history) #tstop is terminal time of the trajectory #cov is the baseline covariates known at time 0
simulator <- function(x0, s0, t0, tstop, cov){ 
  rate_bounds <- rate_bounds(s0, t0, tstop, cov)
  
  times <- c(t0)
  marks <- c(x0)
  lastJump <- s0
  lastState <- x0
  t <- t0
  bx <- rate_bounds[lastState]
  
  if(bx == 0) { break } #break if absorbed
  
  repeat{
    t <- lastJump+jump_simulator(lastState, lastJump, tstop, bx, cov)
    
    if(t > tstop){ break } #break if jump time exceeds terminal time
    
    x <- sample(1:no_states, 1, prob = mark_dist(lastState, lastJump, t, cov))
    
    times <- c(times, t)
    marks <- c(marks, x)
    lastJump <- t
    lastState <- x
    bx <- rate_bounds[lastState]
    
    if(bx == 0){ break } #break if absorbed
    
  }
  return(list(times = times, marks = marks))
}