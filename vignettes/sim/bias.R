setwd("Q:/Paul/wCorr/wCorr/vignettes")
source("sim.R")

grid <- expand.grid(ML=FALSE,
                    iter=1:50,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    fast=TRUE)
grid$reset <- TRUE
bias <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)

save(bias, file="bias.RData")
