require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim.R")

grid1 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=80,
                     n = round(10^seq(1,4.75,by=0.25)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(5,7,by=0.25)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid1, grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)

save(speed, file="speed.RData")
