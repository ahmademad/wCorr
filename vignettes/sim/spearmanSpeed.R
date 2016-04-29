###################
require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("Spearman_consistent_sim.R")

grid <- expand.grid(usew=c(FALSE,TRUE),
                    iter=1:2,
                    n = c(10,100,1000, 10000),
                    rho = c(-0.99,seq(-0.95,0,by=0.05)))

grid$reset <- !grid$usew
spear <- spearmanSim(n=grid$n, rho=grid$rho, usew=grid$usew, outstr="spear")
save(spear, file="spear.RData")
