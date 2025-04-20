setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#states: 1 (active), 2 (disabled), 3 (reactivated), 4 (dead)
#t is time since start, u is duration in current state, x is baseline covariate (age0,time0,gender) with gender==1 being male
mu12 <- function(t,u,x){exp(0.023*(x[1]+t)+0.304*(x[2])-7.46*x[3]-8.66*(1-x[3]) )}
mu23 <- function(t,u,x){exp(-0.012*(x[1]+t)-0.125*(x[2])+0.334*x[3]+0.756*(1-x[3])-1.04*u)}
mu14 <- function(t,u,x){exp(0.07*(x[1]+t)-9.50*x[3]-9.80*(1-x[3]))}
mu24 <- function(t,u,x){exp(0.07*(x[1]+t)-6.40*x[3]-6.80*(1-x[3])-0.25*pmin(u,5))}
mu34 <- function(t,u,x){mu14(t,u,x)}


l <- 10 #steps per year
h <- 1/l #step size
x0 <- c(30,5,1) 
u0 <- 0
t0 <- 0 #initial time
t <- 90 #end time
N <- (t-t0)*l #number of steps required


#calculate p_1k(0,s,0,u; x0); first entry in matrix is s and second entry is u


kolmogorov_forward <- function(x0){
  
  #initialize time zero and first step
  p11 <- matrix(nrow = N+1, ncol = N+1)
  p12 <- matrix(nrow = N+1, ncol = N+1)
  p13 <- matrix(nrow = N+1, ncol = N+1)
  p14 <- matrix(nrow = N+1, ncol = N+1)
  
  p11[1, 1] <- 1; p11[2:(N+1), 1] <- 0 
  p12[1, 1] <- 0; p12[2:(N+1), 1] <- 0
  p13[1, 1] <- 0; p13[2:(N+1), 1] <- 0
  p14[1, 1] <- 0; p14[2:(N+1), 1] <- 0
  
  #solve first step
  p11[2, 2] <- p11[1, 1] - h * (mu12(t0,u0,x0) + mu14(t0,u0,x0))
  p12[2, 2] <- p12[1, 1] + h * mu12(t0,u0,x0)
  p13[2, 2] <- p13[1, 1] + h * 0
  p14[2, 2] <- p14[1, 1] + h * mu14(t0,u0,x0)
  
  #solve after first step
  for (n in 1:(N-1)) {
    print(t0 + h * n)
    for (d in (-n):0) {
      # Midpoints for better integral approximations
      u_vals <- u0 + h * ((1:(n+1)) - 0.5)
      u_vals_d <- u0 + h * ((1:(d+n+1)) - 0.5)
      t_now <- t0 + h * n
      
      # transition probability measures as a function of duration
      dp11 <- diff(c(0, p11[n+1, 1:(n+1)]))
      dp12 <- diff(c(0, p12[n+1, 1:(n+1)]))
      dp11_d <- diff(c(0, p11[n+1, 1:(d+n+1)]))
      dp12_d <- diff(c(0, p12[n+1, 1:(d+n+1)]))
      
      # Update p11
      p11[n+2, d+n+2] <- p11[n+1, d+n+1] - 
        h * sum(dp11_d * (mu12(t_now, u_vals_d, x0) + mu14(t_now, u_vals_d, x0)))
      
      # Update p12
      p12[n+2, d+n+2] <- p12[n+1, d+n+1] + 
        h * p11[n+1, n+1] * mu12(t_now, 0, x0) - 
        h * sum(dp12_d * (mu23(t_now, u_vals_d, x0) + mu24(t_now, u_vals_d, x0)))
      
      # Update p13
      p13[n+2, d+n+2] <- p13[n+1, d+n+1] + 
        h * sum(dp12 * mu23(t_now, u_vals, x0)) - 
        h * p13[n+1, d+n+1] * mu34(t_now, 0, x0)
      
      # Update p14
      p14[n+2, d+n+2] <- p14[n+1, d+n+1] + 
        h * p11[n+1, n+1] * mu14(t_now, 0, x0) + 
        h * sum(dp12 * mu24(t_now, u_vals, x0)) + 
        h * p13[n+1, n+1] * mu34(t_now, 0, x0)
    }
  }
  
  return(list(p11,p12,p13,p14))
}

p1kM = kolmogorov_forward(x0)

p11M <- p1kM[[1]]
p12M <- p1kM[[2]]
p13M <- p1kM[[3]]
p14M <- p1kM[[4]]

##compute probabilities (sanity checks):
#p14M[1+2*l, 1+l] #probability of being in the dead state at time 2 with duration less than 1.
#p14M[1+2*l, 1+2*l] - p14M[1+l, 1+l] #Probability of being in the dead state at time 2 minus probability of being in the dead state at time 1. 

#diag(p11M) + diag(p12M) + diag(p13M) + diag(p14M) #=1 for all times s

x0 <- c(30,5,0) 

p1kF = kolmogorov_forward(x0)

p11F <- p1kF[[1]]
p12F <- p1kF[[2]]
p13F <- p1kF[[3]]
p14F <- p1kF[[4]]

#Plot

png("../Figures/Probabilities.png", width = 10, height = 6, units = 'in', res = 300)
par(mfrow=c(1,2))
par(mar = c(4, 4, 2, 2))

plot(x0[1]+seq(t0, t, h), diag(p11M), type = "l", 
     xlab = "Age", ylab = "Probability", ylim = c(0,1), main="Transition probabilities, 30 year old male")
lines(x0[1]+seq(t0, t, h), diag(p12M), col="gray50")
lines(x0[1]+seq(t0, t, h), diag(p13M), lty = 2)
lines(x0[1]+seq(t0, t, h), diag(p14M), col="gray50", lty = 2)
legend(65,1.05,  
       legend=c("Active", "Disabled", "Reactivated", "Dead"), 
       lty = c(1,1,2,2), 
       col=c("black","gray50","black","gray50"),
       box.lty=0,
       inset=.001,
       bg="transparent")

plot(x0[1]+seq(t0, t, h), diag(p11F), type = "l", 
     xlab = "Age", ylab = "", ylim = c(0,1), main="Transition probabilities, 30 year old female")
lines(x0[1]+seq(t0, t, h), diag(p12F), col="gray50")
lines(x0[1]+seq(t0, t, h), diag(p13F), lty = 2)
lines(x0[1]+seq(t0, t, h), diag(p14F), col="gray50", lty = 2)
legend(65,1.05,  
       legend=c("Active", "Disabled", "Reactivated", "Dead"), 
       lty = c(1,1,2,2), 
       col=c("black","gray50","black","gray50"),
       box.lty=0,
       inset=.001,
       bg="transparent")

dev.off()