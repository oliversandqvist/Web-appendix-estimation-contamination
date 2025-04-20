library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file="../Results/dfParamEstC&K.Rda")

##TABLE

#theta, f, g
mean(dataParam$f4, na.rm=TRUE)-1
sd(dataParam$f4, na.rm=TRUE)

mean(dataParam$f5, na.rm=TRUE)-1.5
sd(dataParam$f5, na.rm=TRUE)

mean(dataParam$f6, na.rm=TRUE)-0.2
sd(dataParam$f6, na.rm=TRUE)

mean(dataParam$theta7app, na.rm=TRUE)-(-0.3)
sd(dataParam$theta7app, na.rm=TRUE)
