## ----packages and data, echo=FALSE, results="hide", message=FALSE,warning=FALSE----
require(knitr)
require(wCorr)
require(lattice)
require(doBy)
load("../vignettes/sim/ML.RData")
load("../vignettes/sim/fast.RData")
load("../vignettes/sim/speed.RData")

## ----ML MAD plot, echo=FALSE,fig.width=6, fig.height=2.5-----------------
ml <- subset(ML, type %in% c("Polychoric", "Polyserial"))
ml$i <- rep(1:(nrow(ml)/2),each=2)
mml <- merge(subset(ml,ML),
               subset(ml,!ML, c("i", "est")),
               by="i",
               suffixes=c(".ml",".nonml"))
mml$ML <- NULL
mml$rmsedrho <- (mml$est.ml - mml$est.nonml)^2
mml$drho <- mml$est.ml - mml$est.nonml
mml$absdrho <- abs(mml$est.ml - mml$est.nonml)

agg <- summaryBy(absdrho ~ n + rho + type, data=mml, FUN=mean, na.rm=TRUE)
xyplot(absdrho.mean ~ rho|type,
       data=agg,
       groups=n,
       scales=list(y=list(log=10, cex=0.7), x = list(cex=0.7)),
       type="l",
       ylab="MAD",
       xlab=expression(rho),
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----ML MAD table, echo=FALSE--------------------------------------------
agg <- summaryBy(absdrho ~ type+n, data=mml, FUN=mean, na.rm=TRUE)
colnames(agg) <- c("Correlation type", "n", "MAD")
kable(agg)

## ----fast MAD plot, echo=FALSE,fig.width=6, fig.height=3.5---------------
fast$i <- rep(1:(nrow(fast)/2),each=2)
mfast <- merge(subset(fast,fast),
               subset(fast,!fast, c("i", "est")),
               by="i",
               suffixes=c(".fast",".slow"))
mfast$fast <- NULL
mfast$drho <- mfast$est.fast - mfast$est.slow
mfast$absdrho <- abs(mfast$est.fast - mfast$est.slow) + 1E-17

agg <- summaryBy(absdrho ~ n + rho + type, data=mfast, FUN=mean, na.rm=TRUE)
xyplot(absdrho.mean ~ rho|type,
       data=agg,
       groups=n,
       type="l",
       ylab="MAD",
       scales=list(y=list(log=10, cex=0.7), x=list(cex=0.7)),
       xlab=expression(rho),
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2))
       )

## ----plot speed, echo=FALSE,fig.width=6, fig.height=3.5------------------
speed$class <- ifelse(speed$ML, "ML=T,", "ML=F,")
speed$class <- paste0(speed$class, ifelse(speed$fast, "fast=T", "fast=F"))
speed$t <- pmax(speed$t, 0.001)
agg <- summaryBy(t ~ n + type + class, data=speed, FUN=mean, na.rm=TRUE)
xyplot(t.mean ~ n|type,
       data=agg,
       type="l",
       ylab="Computing Time",
       scales=list(y=list(log=10, cex=0.7), x=list(log=10, cex=0.7)),
       xlab="n",
       groups=class,
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2))
       )

