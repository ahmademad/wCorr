require(wCorr)
setwd("Q:/Paul/wCorr/wCorr/vignettes")
source("sim.R")

grid <- expand.grid(fast=c(TRUE,FALSE),
                    iter=1:10,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    ML=FALSE)
grid$reset <- grid$fast
fast <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)

save(fast, file="fast.RData")
