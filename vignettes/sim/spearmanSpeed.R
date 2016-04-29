###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("Spearman_consistent_sim.R")

grid <- expand.grid(usew=c(FALSE,TRUE),
                    iter=1:2,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0,by=0.05)))

grid$reset <- !grid$usew
spear1 <- spearmanSim(n=grid$n, rho=grid$rho, usew=grid$usew, outstr="spear1")
save(spear1, file="spear1.RData")

###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("Spearman_consistent_sim.R")

grid <- expand.grid(usew=c(FALSE,TRUE),
                    iter=1:2,
                    n = c(10000),
                    rho = c(-0.99,seq(-0.95,0,by=0.05)))

grid$reset <- !grid$usew
spear2 <- spearmanSim(n=grid$n, rho=grid$rho, outstr="spear2")
save(spear2, file="spear2.RData")

#######
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
load("spear1.RData")
spear <- spear1
load("spear2.RData")
spear <- rbind(spear, spear2)
save(spear, file="spear.RData")


