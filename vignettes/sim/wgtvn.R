require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid <- expand.grid(usew=c(FALSE,TRUE),
                    iter=1:20,
                    n = c(10,100,1000,10000),
                    rho = c(-0.99,seq(-0.95,0,by=0.05)))
grid$reset <- !grid$usew

wgtvn <- wCorrSim(n=grid$n, rho=grid$rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=grid$usew)
save(wgtvn, file="wgtvn.RData")
