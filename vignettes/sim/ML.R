require(wCorr)
setwd("Q:/Paul/wCorr/vignettes/sim")
source("sim2.R")

grid <- expand.grid(ML=c(TRUE,FALSE),
                    iter=1:500,
                    n = c(10,100,1000),
                    rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                    fast=TRUE)
grid$reset <- grid$ML
ML <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)

save(ML, file="ML.RData")

require(doBy)

ml <- subset(ML, type %in% c("Polychoric", "Polyserial"))
ml$rmse <- (ml$est - ml$rho)^2

aggml <- summaryBy(rmse ~ n + rho + type + ML, data=ml, FUN=mean, na.rm=TRUE)
aggml$rmse.mean <- sqrt(aggml$rmse.mean)
aggml$ml <- ifelse(aggml$ML==TRUE, "ML=TRUE", "ML=FALSE")
aggml$nt <- factor(paste("n=",aggml$n))


ml$i <- rep(1:(nrow(ml)/2),each=2)
mml <- merge(subset(ml,ML),
             subset(ml,!ML, c("i", "est")),
             by="i",
             suffixes=c(".ml",".nonml"))
mml$absd <- abs(mml$est.ml - mml$est.nonml)
aggt1_0 <- summaryBy(absd ~ type + n + ML, data=subset(mml, type=="Polychoric"), FUN=mean, na.rm=TRUE)
aggt1_0$ML <- NULL

aggt1 <- summaryBy(rmse ~ type + n + ML, data=subset(ml, type=="Polychoric"), FUN=mean, na.rm=TRUE)



aggt2_0 <- summaryBy(absd ~ type + n + ML, data=subset(mml, type=="Polyserial"), FUN=mean, na.rm=TRUE)
aggt2_0$ML <- NULL

aggt2 <- summaryBy(rmse ~ type + n + ML, data=subset(ml, type=="Polyserial"), FUN=mean, na.rm=TRUE)
aggt2$rmse.mean <- sqrt(aggt2$rmse.mean)

save(aggml, aggt1_0, aggt1, aggt2_0, aggt2, file="aggML.RData")
