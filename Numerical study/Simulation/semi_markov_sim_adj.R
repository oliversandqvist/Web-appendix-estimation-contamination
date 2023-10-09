#Reference: "On Lewis' simulation method for point processes". For more details, "Exact simulation of the Jump Times of a Class of Piecewise Deterministic Markov Processes".

no_states <- 3

#Computes the transition rate for a given transition at a given time and duration
#x is current state #y is destination state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0 (we use age0 = age at time 0 and gender is gender where female is 0 and male is 1)
transition_rate_adj <-  function(x, y, s, t, cov){
  if(x == 1 && y==2){
    as.numeric(0.8*(cov/(t-s+2))^2)
  }
  else if(x == 2 && y==3){
    as.numeric(exp(-1.2*(t-s)))
  }
  else{
    0
  }
}

#Computes the total intensity of out the current state
#x is current state #s is time of last jump #t is current time #cov is the baseline covariates known at time 0
jump_rate_adj <-  function(x, s, t, cov){
  sum(sapply(1:no_states, FUN = function(y) transition_rate_adj(x,y,s,t,cov)))
}

#Computes the mark distribution given a jump time
#x is current state #s is time of last jump #t is time of current jump
mark_dist_adj <- function(x, s, t, cov){
  rates <- sapply(1:no_states, FUN = function(y) transition_rate_adj(x,y,s,t,cov))
  rates/sum(rates)
}

#Simulates time from the last jump to the next jump 
#x is current state #s is time of last jump #b is bound #cov is the baseline covariates known at time 0
jump_simulator_adj <- function(x, s, tstop, b, cov){ 
  u <- runif(1)
  e <- rexp(1, rate = 1)
  t <- e/b
  if (t > tstop) { return(Inf) }
  while(u > jump_rate_adj(x, s, s+t, cov)/b){ 
    u <- runif(1)
    e <- rexp(1, rate = 1)
    t <- e/b+t
    if (t > tstop) { return(Inf) }
  }
  return(t)
}

#Computes the "local upper bounds" of the intensities for each of the states.
#s is the time of the last jump at time 0 #t is the current time #tstop is the terminal time #cov is the baseline covariates known at time 0
rate_bounds_adj <- function(s, t, tstop, cov){
  
  b1 <- (cov/(0+2))^2
  b2 <- exp(-0)
  b3 <- 0
  rate_bounds <- c(b1,b2,b3)
  
  return(rate_bounds)
}

#Simulates one trajectory of the semi-Markov process.
#x0 is initial state #s0 is last jump time #t0 is current time (note that by convention, we include the zero'th jump-time and mark in the MPP history) #tstop is terminal time of the trajectory #cov is the baseline covariates known at time 0
simulator_adj <- function(x0, s0, t0, tstop, cov){ 
  rate_bounds <- rate_bounds_adj(s0, t0, tstop, cov)
  
  times <- c(t0)
  marks <- c(x0)
  lastJump <- s0
  lastState <- x0
  t <- t0
  bx <- rate_bounds[lastState]
  
  if(bx == 0) { break } #break if absorbed
  
  repeat{
    t <- lastJump+jump_simulator_adj(lastState, lastJump, tstop, bx, cov)
    
    if(t > tstop){ break } #break if jump time exceeds terminal time
    
    x <- sample(1:no_states, 1, prob = mark_dist_adj(lastState, lastJump, t, cov))
    
    times <- c(times, t)
    marks <- c(marks, x)
    lastJump <- t
    lastState <- x
    bx <- rate_bounds[lastState]
    
    if(bx == 0){ break } #break if absorbed
    
  }
  return(list(times = times, marks = marks))
}