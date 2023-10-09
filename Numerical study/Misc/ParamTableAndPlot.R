library(tidyverse)

load(file="Results/dfParamEst.Rda")

rmse <- function(true,obs){
 obsNotNA <- na.omit(obs)
 return(sqrt(1/length(obsNotNA)*sum((true-obsNotNA)^2)))
}


##TABLE

#theta, f, g
mean(dfParam$g1, na.rm=TRUE)-0.8
sd(dfParam$g1, na.rm=TRUE)
rmse(0.8,dfParam$g1)

mean(dfParam$g2, na.rm=TRUE)-(-1.2)
sd(dfParam$g2, na.rm=TRUE)
rmse(-1.2,dfParam$g2)

mean(dfParam$f1, na.rm=TRUE)-2
sd(dfParam$f1, na.rm=TRUE)
rmse(2,dfParam$f1)

mean(dfParam$f2, na.rm=TRUE)-0.5
sd(dfParam$f2, na.rm=TRUE)
rmse(0.5,dfParam$f2)

mean(dfParam$f3, na.rm=TRUE)-0.1
sd(dfParam$f3, na.rm=TRUE)
rmse(0.1,dfParam$f3)

mean(dfParam$f4, na.rm=TRUE)-1
sd(dfParam$f4, na.rm=TRUE)
rmse(1,dfParam$f4)

mean(dfParam$f5, na.rm=TRUE)-1.5
sd(dfParam$f5, na.rm=TRUE)
rmse(1.5,dfParam$f5)

mean(dfParam$f6, na.rm=TRUE)-0.2
sd(dfParam$f6, na.rm=TRUE)
rmse(0.2,dfParam$f6)


mean(dfParam$theta1, na.rm=TRUE)-log(0.15)
sd(dfParam$theta1, na.rm=TRUE)
rmse(log(0.15),dfParam$theta1)

mean(dfParam$theta2, na.rm=TRUE)-0.1
sd(dfParam$theta2, na.rm=TRUE)
rmse(0.1,dfParam$theta2)

mean(dfParam$theta3, na.rm=TRUE)-0.4
sd(dfParam$theta3, na.rm=TRUE)
rmse(0.4,dfParam$theta3)

mean(dfParam$theta4, na.rm=TRUE)-log(0.1)
sd(dfParam$theta4, na.rm=TRUE)
rmse(log(0.1),dfParam$theta4)

mean(dfParam$theta5, na.rm=TRUE)-0.03
sd(dfParam$theta5, na.rm=TRUE)
rmse(0.03,dfParam$theta5)

mean(dfParam$theta6, na.rm=TRUE)-(-0.3)
sd(dfParam$theta6, na.rm=TRUE)
rmse(-0.3,dfParam$theta6)

mean(dfParam$theta7, na.rm=TRUE)-(-0.3)
sd(dfParam$theta7, na.rm=TRUE)
rmse(-0.3,dfParam$theta7)


#theta oracle, theta approx (remove runs that did not converge)
dfParamOracle12NoDiverge <- dfParam
dfParamOracle12NoDiverge <- dfParamOracle12NoDiverge %>% filter(!(theta1oracle > 0))

mean(dfParamOracle12NoDiverge$theta1oracle, na.rm=TRUE)-(log(0.15))
sd(dfParamOracle12NoDiverge$theta1oracle, na.rm=TRUE)
rmse(log(0.15),dfParamOracle12NoDiverge$theta1oracle)

mean(dfParamOracle12NoDiverge$theta2oracle, na.rm=TRUE)-(0.1)
sd(dfParamOracle12NoDiverge$theta2oracle, na.rm=TRUE)
rmse(0.1,dfParamOracle12NoDiverge$theta2oracle)

mean(dfParamOracle12NoDiverge$theta3oracle, na.rm=TRUE)-(0.4)
sd(dfParamOracle12NoDiverge$theta3oracle, na.rm=TRUE)
rmse(0.4,dfParamOracle12NoDiverge$theta3oracle)

mean(dfParam$theta4oracle, na.rm=TRUE)-(log(0.1))
sd(dfParam$theta4oracle, na.rm=TRUE)
rmse(log(0.1),dfParam$theta4oracle)

mean(dfParam$theta5oracle, na.rm=TRUE)-(0.03)
sd(dfParam$theta5oracle, na.rm=TRUE)
rmse(0.03,dfParam$theta5oracle)

mean(dfParam$theta6oracle, na.rm=TRUE)-(-0.3)
sd(dfParam$theta6oracle, na.rm=TRUE)
rmse(-0.3,dfParam$theta6oracle)

mean(dfParam$theta7oracle, na.rm=TRUE)-(-0.3)
sd(dfParam$theta7oracle, na.rm=TRUE)
rmse(-0.3,dfParam$theta7oracle)


mean(dfParam$theta1app, na.rm=TRUE)-(log(0.15))
sd(dfParam$theta1app, na.rm=TRUE)
rmse(log(0.15),dfParam$theta1app)

mean(dfParam$theta2app, na.rm=TRUE)-(0.1)
sd(dfParam$theta2app, na.rm=TRUE)
rmse(0.1,dfParam$theta2app)

mean(dfParam$theta3app, na.rm=TRUE)-(0.4)
sd(dfParam$theta3app, na.rm=TRUE)
rmse(0.4,dfParam$theta3app)

mean(dfParam$theta4app, na.rm=TRUE)-(log(0.1))
sd(dfParam$theta4app, na.rm=TRUE)
rmse(log(0.1),dfParam$theta4app)

mean(dfParam$theta5app, na.rm=TRUE)-(0.03)
sd(dfParam$theta5app, na.rm=TRUE)
rmse(0.03,dfParam$theta5app)

mean(dfParam$theta6app, na.rm=TRUE)-(-0.3)
sd(dfParam$theta6app, na.rm=TRUE)
rmse(-0.3,dfParam$theta6app)

mean(dfParam$theta7app, na.rm=TRUE)-(-0.3)
sd(dfParam$theta7app, na.rm=TRUE)
rmse(-0.3,dfParam$theta7app)


#theta Naive 1 (remove runs that did not converge)
dfParamCens12NoDiverge <- dfParam
dfParamCens12NoDiverge <- dfParamCens12NoDiverge %>% filter(!(theta1cens0 > 0))

mean(dfParamCens12NoDiverge$theta1cens0, na.rm=TRUE)-(log(0.15))
sd(dfParamCens12NoDiverge$theta1cens0, na.rm=TRUE)
rmse(log(0.15),dfParamCens12NoDiverge$theta1cens0)

mean(dfParamCens12NoDiverge$theta2cens0, na.rm=TRUE)-(0.1)
sd(dfParamCens12NoDiverge$theta2cens0, na.rm=TRUE)
rmse(0.1,dfParamCens12NoDiverge$theta2cens0)

mean(dfParamCens12NoDiverge$theta3cens0, na.rm=TRUE)-(0.4)
sd(dfParamCens12NoDiverge$theta3cens0, na.rm=TRUE)
rmse(0.4,dfParamCens12NoDiverge$theta3cens0)

mean(dfParam$theta4cens0, na.rm=TRUE)-(log(0.1))
sd(dfParam$theta4cens0, na.rm=TRUE)
rmse(log(0.1),dfParam$theta4cens0)

mean(dfParam$theta5cens0, na.rm=TRUE)-(0.03)
sd(dfParam$theta5cens0, na.rm=TRUE)
rmse(0.03,dfParam$theta5cens0)

mean(dfParam$theta6cens0, na.rm=TRUE)-(-0.3)
sd(dfParam$theta6cens0, na.rm=TRUE)
rmse(-0.3,dfParam$theta6cens0)

mean(dfParam$theta7cens0, na.rm=TRUE)-(-0.3)
sd(dfParam$theta7cens0, na.rm=TRUE)
rmse(-0.3,dfParam$theta7cens0)


#theta Naive 2 (remove runs that did not converge)
dfParamCens112NoDiverge <- dfParam
dfParamCens112NoDiverge <- dfParamCens112NoDiverge %>% filter(!(theta1cens1 > 0))


mean(dfParamCens112NoDiverge$theta1cens1, na.rm=TRUE)-(log(0.15))
sd(dfParamCens112NoDiverge$theta1cens1, na.rm=TRUE)
rmse(log(0.15),dfParamCens112NoDiverge$theta1cens1)

mean(dfParamCens112NoDiverge$theta2cens1, na.rm=TRUE)-(0.1)
sd(dfParamCens112NoDiverge$theta2cens1, na.rm=TRUE)
rmse(0.1,dfParamCens112NoDiverge$theta2cens1)

mean(dfParamCens112NoDiverge$theta3cens1, na.rm=TRUE)-(0.4)
sd(dfParamCens112NoDiverge$theta3cens1, na.rm=TRUE)
rmse(0.4,dfParamCens112NoDiverge$theta3cens1)

mean(dfParam$theta4cens1, na.rm=TRUE)-(log(0.1))
sd(dfParam$theta4cens1, na.rm=TRUE)
rmse(log(0.1),dfParam$theta4cens1)

mean(dfParam$theta5cens1, na.rm=TRUE)-(0.03)
sd(dfParam$theta5cens1, na.rm=TRUE)
rmse(0.03,dfParam$theta5cens1)

mean(dfParam$theta6cens1, na.rm=TRUE)-(-0.3)
sd(dfParam$theta6cens1, na.rm=TRUE)
rmse(-0.3,dfParam$theta6cens1)

mean(dfParam$theta7cens1, na.rm=TRUE)-(-0.3)
sd(dfParam$theta7cens1, na.rm=TRUE)
rmse(-0.3,dfParam$theta7cens1)


##PLOTS

#histogram
png("Figures/Histogram.png", width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(5,3))
par(mar = c(2, 2, 0, 2))


hist(dfParam$g1, main="")
abline(v=0.8, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('g'[1])), bty = "n")

hist(dfParam$g2, main="")
abline(v=-1.2, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('g'[2])), bty = "n")

hist(dfParam$f1, main="")
abline(v=2, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[1])), bty = "n")

hist(dfParam$f2, main="")
abline(v=0.5, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[2])), bty = "n")

hist(dfParam$f3, main="")
abline(v=0.1, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[3])), bty = "n")

hist(dfParam$f4, main="")
abline(v=1, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[4])), bty = "n")

hist(dfParam$f5, main="")
abline(v=1.5, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[5])), bty = "n")

hist(dfParam$f6, main="")
abline(v=0.2, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('f'[6])), bty = "n")

hist(dfParam$theta1, main="")
abline(v=log(0.15), col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[1])), bty = "n")

hist(dfParam$theta2, main="")
abline(v=0.1, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[2])), bty = "n")

hist(dfParam$theta3, main="")
abline(v=0.4, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[3])), bty = "n")

hist(dfParam$theta4, main="")
abline(v=log(0.1), col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[4])), bty = "n")

hist(dfParam$theta5, main="")
abline(v=0.03, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[5])), bty = "n")

hist(dfParam$theta6, main="")
abline(v=-0.3, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[6])), bty = "n")

hist(dfParam$theta7, main="")
abline(v=-0.3, col="black", lty=c(2), lwd=c(3))
legend("topleft", legend=expression(paste('θ'[7])), bty = "n")
dev.off()