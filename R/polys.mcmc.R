polys.mcmc <- function(x,y,w=rep(1,length(x)), verbose=FALSE, nmax=1e3) {
  t1 <- list(xbar=Inf,xsd=Inf,wlpx=Inf)

  lnl <- function(x_, y_, w_, y_c, x_bar, x_sd, corr_) {
    # following notation in Olsson, Drasgow, and Dorans, 1982
    y_c <- c(-Inf,y_c,Inf)
    ti <- get("t1")
    z <- (x_ - x_bar)/x_sd
    if(ti$xbar==x_bar & ti$xsd==x_sd) {
      # we're good
    } else {
      lpx <- dnorm(x_, mean=x_bar, sd=x_sd, log=TRUE)
      wlpx <- sum(w_ * lpx)
      ti$xbar <- x_bar
      ti$xsd <- x_sd
      ti$wlpx <- wlpx
      t1 <<- ti
    }
    tj <- y_c[y_+1]
    tjm1 <- y_c[y_]
    tjstar <- (tj - corr_*z) / ((1-corr_^2)^0.5)
    tjm1star <- (tjm1 - corr_*z) / ((1-corr_^2)^0.5)
    lpy <- log( pnorm(tjstar) - pnorm(tjm1star))
    lpy[lpy==-Inf] <- min(lpy[lpy>-Inf])
    ti$wlpx + sum(w_*lpy)
  }

  yi <- as.integer(y)
  yc <- 1:(length(unique(yi))-1)
  yc <- yc - mean(yc)
  mux <- mean(x)
  sdx <- sd(x)
  corr <- cor(yi,x)
  lnl0 <- lnl(x, y, w, yc, mux, sdx, corr)
  hist <- data.frame(mux=rep(NA,nmax),
                     sdx=rep(NA,nmax),
                     corr=rep(NA,nmax),
                     lnl=rep(NA,nmax))
  factx <- 1
  factsd <- 1
  factcorr <- 1
  facty <- rep(1,length(yc))

  for(i in 1:nmax) {
    muxp <- mux + rnorm(1,0,1) * sdx / sqrt(length(y)) * factx
    lnli <- lnl(x, y, w, yc, muxp, sdx, corr)
    a <- exp(lnli - lnl0)
    if(runif(1) < a) {
      mux <- muxp
      lnl0 <- lnli
      if(runif(1) < 0.2) {
        factx <- factx * 2
      }
    } else {
      if(runif(1) < 0.2) {
        factx <- factx * 0.5
      }
    }

    sdxp <- abs(sdx + rnorm(1,0,1/100) * factsd)
    lnli <- lnl(x, y, w, yc, mux, sdxp, corr)
    a <- exp(lnli - lnl0)
    if(runif(1) < a) {
      sdx <- sdxp
      lnl0 <- lnli
      if(runif(1) < 0.2) {
        factsd <- factsd * 2
      }
    } else {
      if(runif(1) < 0.2) {
        factsd <- factsd * 0.5
      }
    }

    corrp <- tanh( atanh(corr) + rnorm(1,0,0.01) * factcorr)
    lnli <- lnl(x, y, w, yc, mux, sdx, corrp)
    a <- exp(lnli - lnl0)
    if(runif(1) < a) {
      corr <- corrp
      lnl0 <- lnli
      if(runif(1) < 0.2) {
        factcorr <- factcorr * 2
      }
    } else {
      if(runif(1) < 0.2) {
        factcorr <- factcorr * 0.5
      }
    }

    for(yci in 1:length(yc)) {
      ycp <- yc
      ycp[yci] <- ycp[yci] + rnorm(1,0,sd=0.01) * facty[yci]
      lnli <- lnl(x, y, w, ycp, mux, sdx, corr)
      a <- exp(lnli - lnl0)
      if(runif(1) < a) {
        yc <- ycp
        lnl0 <- lnli
        if(runif(1) < 0.2) {
          facty[yci] <- facty[yci] * 2
        }
      } else {
        if(runif(1) < 0.2) {
          facty[yci] <- facty[yci] * 0.5
        }
      }
    }
    hist$mux[i] <- mux
    hist$sdx[i] <- sdx
    hist$corr[i] <- corr
    hist$lnl[i] <- lnl0
  }

  mean(hist$corr)
}


polys <- function(x,y,w=rep(1,length(x)), verbose=FALSE) {
  t1 <- list(xbar=Inf,xsd=Inf,wlpx=Inf)

  lnl <- function(x_, y_, w_, y_c, x_bar, x_sd, corr_) {
    # following notation in Olsson, Drasgow, and Dorans, 1982
    y_c <- c(-Inf,y_c,Inf)
    ti <- get("t1")
    z <- (x_ - x_bar)/x_sd
    if(ti$xbar==x_bar & ti$xsd==x_sd) {
      # we're good
    } else {
      lpx <- dnorm(x_, mean=x_bar, sd=x_sd, log=TRUE)
      wlpx <- sum(w_ * lpx)
      ti$xbar <- x_bar
      ti$xsd <- x_sd
      ti$wlpx <- wlpx
      t1 <<- ti
    }
    tj <- y_c[y_+1]
    tjm1 <- y_c[y_]
    tjstar <- (tj - corr_*z) / ((1-corr_^2)^0.5)
    tjm1star <- (tjm1 - corr_*z) / ((1-corr_^2)^0.5)
    lpy <- log( pnorm(tjstar) - pnorm(tjm1star))
    lpy[lpy==-Inf] <- min(lpy[lpy>-Inf])
    ti$wlpx + sum(w_*lpy)
  }

  optf_all <- function(par, x, y, w) {
    print(par)
    c1 <- length(unique(y)) - 1
    lnl(x_=x, y_=y, w_=w, y_c=fscale_cuts(par[1:c1]), x_bar=par[c1+1], x_sd=fscale_sd(par[c1+2]), corr_=fscale_corr(par[c1+3]))
  }

  optf_corr <- function(par, parm1, x, y, w) {
    c1 <- length(unique(y)) - 1
    lnl(x_=x, y_=y, w_=w, y_c=fscale_cuts(parm1[1:c1]), x_bar=parm1[c1+1], x_sd=fscale_sd(parm1[c1+2]), corr_=fscale_corr(par))
  }

  optf_cut <- function(par, parm1, x, y, w) {
    lnl(x_=x, y_=y, w_=w, y_c=fscale_cuts(par), x_bar=parm1[1], x_sd=fscale_sd(parm1[2]), corr_=fscale_corr(parm1[3]))
  }

  optf_xpar <- function(par, parm1, x, y, w) {
    c1 <- length(unique(y)) - 1
    lnl(x_=x, y_=y, w_=w, y_c=fscale_cuts(parm1[1:c1]), x_bar=par[1], x_sd=fscale_sd(par[2]), corr_=fscale_corr(parm1[c1+1]))
  }

  fscale_cuts <- function(par) {
    cumsum(c(par[1],exp(par[-1])))
  }

  fscale_corr <- function(par) {
    tanh(par)
  }

  fscale_sd <- function(par) {
    exp(par)
  }
  yi <- as.integer(y)
  yc0 <- 1:(length(unique(yi))-1) - 2
  dx <- Inf
  mux <- mean(x)
  sdx <- sd(x)
  cor0 <- cor(yi,x)
  cor0m1 <- -1*sign(cor0)*Inf # a long way for cor0
  niter <- 1
  maxit <- 20 * length(yc0)
  while(dx > 1E-6 & niter <100) {
    # optimize y cuts
    if(length(yc0) > 1) {
      op <- optim(par=yc0, parm1=c(mux, sdx, cor0), optf_cut, x=x, y=yi, w=w, control=list(fnscale=-1, maxit=maxit))
      yc0 <- op$par
    } else {
      op <- optimize(optf_cut, interval=c(yc0-2, yc0+2), parm1=c(mux, sdx, cor0), x=x, y=yi, w=w, maximum=TRUE)
      yc0 <- op$maximum
    }
    # optimize x mean and se
    op <- optim(par=c(mux, sdx), parm1=c(yc0, cor0), optf_xpar, x=x, y=yi, w=w, control=list(fnscale=-1, maxit=maxit))
    mux <- op$par[1]
    sdx <- op$par[2]
    # optimize correlation
    op <- optimize(optf_corr, interval=c(cor0-2, cor0+2), parm1=c(yc0, mux, sdx), x=x, y=yi, w=w, maximum=TRUE)
    dx2 <- abs(fscale_corr(cor0m1) - fscale_corr(op$maximum))
    cor0m1 <- cor0
    dx <- abs(fscale_corr(cor0) - fscale_corr(op$maximum))
    if(abs(dx2) < abs(dx)) { # looks like a loop, skip some itterations
      niter <- niter * 2
      cor0 <- mean(c(cor0, op$maximum)) # split the difference
      maxit <- 20*length(yc0)
    } else {
      niter <- niter + 1
      maxit <- 10
    }
    cor0 <- op$maximum
    if(verbose) {
      cat("cor0=",cor0, " (",fscale_corr(cor0),") yc=",fscale_cuts(yc0)," mux=",mux,"sdx=",fscale_sd(sdx)," niter=",niter,"dx=",dx," dx2=",dx2,"\n")
    }
  }
  fscale_corr(cor0)
}

if(FALSE) {
  nr <- 2100
  df <- data.frame(n=rep(0,nr),
                   cor=rep(NA,nr),
                   cuts1=rep(NA,nr),
                   cuts2=rep(NA,nr),
                   ps1c=rep(NA,nr),
                   psm=rep(NA,nr),
                   Pe=rep(NA,nr))
  corv <- seq(-0.95,0.95,len=21)

  for(i in 1:nr) {
    cat("i=",i,"\n")
    n <- 1e4
    cuts <- c(-0.2,0.2)
    if(runif(1) < 0.2) { cuts <- c(-1,1) }
    if(runif(1) < 0.2) { cuts <- c(-1,0,1) }
    if(runif(1) < 0.2) { cuts <- c(-0.5,0,0.5) }

    cor <- sample(corv,1)
    df$cuts1[i] <- cuts[1]
    df$cuts2[i] <- cuts[2]
    df$cor[i] <- cor
    df$n[i] <- n

    S <- matrix(c(1,cor,cor,1), nrow=2)

    xy <- mvrnorm(n=n, mu=c(0,0), Sigma=S)
    x <- xy[,1]
    y <- length(cuts)+1
    for(j in rev(1:length(cuts))) {
      y <- ifelse(xy[,2] <= cuts[j], j, y)
    }

    df$ps1c[i] <- polys(x,y)
    #  df$psm[i]<- polys.mcmc(x,y)
    df$Pe[i] <- cor(xy)[1,2]
    print(df[i,])

    dfi <- df[1:i,]
    dfi$dps1c <- dfi$ps1c - dfi$cor
    dfi$dPe <- dfi$Pe - dfi$cor
    dfi$dpsm <- dfi$psm - dfi$cor

    df2 <- sqldf("SELECT avg(dps1c) AS dps1c, avg(dPe) as dPe, avg(dpsm) as dpsm, cor FROM dfi GROUP BY cor")
    df2 <- df2[order(df2$cor),]

    plot(c(-1,1), range(c(df2$dps1c, df2$dPe, df2$dpsm)), , type="n", lwd=3, lty=1, col="black")
    abline(h = 0)
    lines(df2$cor, df2$dps1c , lwd=1, lty=2, col="blue")
    lines(df2$cor, df2$dPe   , lwd=1, lty=1, col="orange")
    lines(df2$cor, df2$dpsm  , lwd=1, lty=3, col="green")

  }

  require(lattice)
  df <- df[order(df$cor),]

  plot(df$cor, df$ps1c, type="l", lwd=1, lty=2, col="blue")
  abline(0,1)
  lines(df$cor, df$Pe, lwd=3, lty=1)
}
