require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid <- expand.grid(ML=FALSE,
                    iter=1:5,
                    n = round(10^seq(1,6,by=0.25)),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    fast=TRUE)
ntime <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=TRUE, usew=FALSE, outstr="ntime")

save(ntime, file="ntime.RData")
