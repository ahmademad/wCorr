require(wCorr)
require(mnormt)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim.R")

grid <- expand.grid(usew=c(FALSE,TRUE),
                    iter=1:100,
                    n = c(100),
                    n = c(100),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid$reset <- !grid$usew
wgtvrho <- wCorrSim(n=grid$n, rho=grid$rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=grid$usew)

save(wgtvrho, file="wgtvrho.RData")
