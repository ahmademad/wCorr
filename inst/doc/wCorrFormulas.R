## ----packages and data, echo=FALSE, results="hide", message=FALSE, warning=FALSE----
require(wCorr)
require(lattice)
require(doBy)
load("..//vignettes/sim/ntime.RData")
load("..//vignettes/sim/bias.RData")
load("..//vignettes/sim/wgtvrho.RData")
load("..//vignettes/sim/wgtvn.RData")

## ----theta,echo=FALSE,results="hide",fig.show="hold",fig.width=6, fig.height=2.5----
x <- seq(-3,3,by=0.01) 
y <- dnorm(x)
par0 <- par(no.readonly=TRUE)
par(ann=FALSE)
par(mar=c(5,2,1,1)+0.1)
plot(x,y,type="l",xlab="y",ylab="Density", xaxt="n", yaxt="n")
axis(1,at=c(-2,-0.5,1.6), labels=expression(theta[3],theta[4],theta[5]))
text(x=c(-2.5,-1.25,0.55,2.3),y=0.05, labels=paste0("m=",1:4))
theta <- c(-2,-0.5,1.6)
for(i in 1:3) {
  lines(rep(theta[i],2), c(-1,dnorm(theta[i])))
}
par(ann=TRUE)
par(mgp=c(0.5,0,0))
title(ylab="density")
par(mgp=c(3,1,0))
title(xlab="Y")
par(par0)

## ----bias vs rho, echo=FALSE,fig.width=6, fig.height=2.5-----------------
bias$rmse <- sqrt( (bias$est - bias$rho)^2 )
bias$bias <- bias$est - bias$rho
agg <- summaryBy(bias + rmse  ~ n + rho + type, data=bias, FUN=mean, na.rm=TRUE)

xyplot(bias.mean ~ rho|type,
      data=agg,
      groups=n,
      type="l",
      ylab="Bias",
      xlab=expression(rho),
      scales=list(x=list(cex=0.7), y=list(cex=0.7)),
      auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
      par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----rmse vs rho, echo=FALSE,fig.width=6, fig.height=2.5-----------------
xyplot(rmse.mean ~ rho|type,
      data=agg,
      groups=n,
      scales=list(y=list(log=10, cex=0.7), x=list(cex=0.7)),
      ylab="RMSE",
      xlab=expression(rho),
      type="l",
      auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
      par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----rmse vs n, echo=FALSE,fig.width=6, fig.height=2.5-------------------
agg <- summaryBy(rmse  ~ n+type, data=bias, FUN=mean, na.rm=TRUE)
xyplot(rmse.mean ~ n,
     groups=type,
      data=agg,
      ylab="RMSE",
      xlab="n",
      scales=list(y=list(log=10, cex=0.7), x=list(log=10, cex=0.7)),
      type="l",
      auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
      par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----time vs n, echo=FALSE,fig.width=6, fig.height=2.5-------------------
agg <- summaryBy(t ~ n + type, data=ntime, FUN=mean, na.rm=TRUE)
agg$t.mean <- ifelse(agg$t.mean==0, 0.001,agg$t.mean)

xyplot(t.mean ~ n,
       data=agg,
       scales=list(y=list(log=10, cex=0.7), x=list(log=10, cex=0.7)),
       groups=type,
       type="l",
       ylab="Computing time (s)",
       xlab="n",
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----wgt vs rho plot, echo=FALSE,fig.width=6, fig.height=2.5-------------
wgt <- wgtvrho
wgt$absdrho <- abs(wgt$est - wgt$rho)

agg <- summaryBy(absdrho ~ rho + usew + type, data=wgt, FUN=mean, na.rm=TRUE)
agg$weight <- ifelse(agg$usew, "Weighted", "Unweighted")

xyplot(absdrho.mean ~ rho|type,
       data=agg,
       groups=weight,
       scales=list(y=list(log=10, cex=0.7), x=list(cex=0.7)),
       type="l", 
       ylab="MAD",
       xlab=expression(rho),
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

## ----wgt v n plot, echo=FALSE,fig.width=6, fig.height=2.5----------------
wgt <- wgtvn
wgt$mserho <- (wgt$est - wgt$rho)^2

agg <- summaryBy(mserho ~ n + usew + type, data=wgt, FUN=mean, na.rm=TRUE)
agg$rmserho <- sqrt(agg$mserho)
agg$weight <- ifelse(agg$usew, "Weighted", "Unweighted")

xyplot(rmserho ~ n|type,
       data=agg,
       groups=weight,
       scales=list(y=list(log=10, cex=0.7), x=list(log=10, cex=0.7)),
       type="l",
       ylab="RMSE",
       xlab="n",
       auto.key=list(lines=TRUE, points=FALSE, space="right", cex=0.7),
       par.settings=list(superpose.line=list(lwd=2), plot.line=list(lwd=2)))

