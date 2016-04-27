require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid <- expand.grid(fast=c(TRUE,FALSE),
                    iter=1:10,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    ML=FALSE)
grid$reset <- grid$fast
fast <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE, outstr="fast")

save(fast, file="fast.RData")

if(FALSE) {
  x <- rnorm(10)
  y <- rnorm(10)+x/2
  afast <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=TRUE, ML=FALSE)
  aslow <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=FALSE, ML=FALSE)
  c(afast-aslow,afast,aslow)
  afast2 <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=TRUE, ML=FALSE)
  c(afast-afast2,afast,afast2)
  afast <- weightedCorr(x,y, method="Spearman", weights=rep(1,10), fast=TRUE, ML=FALSE)
  aslow <- weightedCorr(x,y, method="Spearman", weights=rep(1,10), fast=FALSE, ML=FALSE)
  c( (afast-aslow)*10^16,afast,aslow)
  M <- 0 + (x > 0)
  P <- 9 + (y > 0)
  afast <- weightedCorr(x,y, method="polychoric", weights=rep(1,10), fast=TRUE, ML=FALSE)
  aslow <- weightedCorr(x,y, method="polychoric", weights=rep(1,10), fast=FALSE, ML=FALSE)
  c( (afast-aslow),afast,aslow)
}
