require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid1 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=80,
                     n = round(10^seq(1,4.75,by=0.25)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(5,7,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid1, grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed1 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed")

save(speed, file="speed.RData")


##############
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid1 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=80,
                     n = round(10^seq(1,4.75,by=0.25)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid1)

grid$reset <- (grid$ML) & (grid$fast)
speed1 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed1")

save(speed1, file="speed1.RData")

###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(5,6,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed2 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed2")
save(speed2, file="speed2.RData")


###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid2 <- expand.grid(fast=c(TRUE,FALSE),
                     ML=c(TRUE,FALSE),
                     iter=20,
                     n = round(10^seq(6.5,7,by=0.5)),
                     rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
grid <- rbind(grid2)

grid$reset <- (grid$ML) & (grid$fast)
speed3 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed3")
save(speed3, file="speed3.RData")

####
load("speed1.RData")
load("speed2.RData")
load("speed3.RData")
speed <- rbind(speed1, speed2, speed3)
save(speed, file="speed.RData")


