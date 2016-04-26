require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid <- expand.grid(ML=c(TRUE,FALSE),
                    iter=1:500,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    fast=TRUE)
grid$reset <- grid$ML
ML <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)

save(ML, file="ML.RData")
