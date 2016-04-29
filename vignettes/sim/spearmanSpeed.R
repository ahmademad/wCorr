###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("Spearman_consistent_sim.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(1,4.75,by=0.25)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed1 <- spearmanSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed2")
save(speed1, file="speedSpear1.RData")


###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("Spearman_consistent_sim.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(5,6,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed2 <- spearmanSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed2")
save(speed2, file="speedSpear2.RData")

###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(6.5,6.5,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed3 <- spearmanSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed3")
save(speed3, file="speedSpear3.RData")


###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=10,
                     n = round(10^seq(7,7,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed4 <- spearmanSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed4")
save(speed4, file="speedSpear4.RData")