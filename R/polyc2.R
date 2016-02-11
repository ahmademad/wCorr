# based losely on Olsson, Ulf (1979), "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient", Psychometrica, 44(4), 443-460.
polyc2 <- function(x,y,w) {
  lnl <- function(xytab, cc, rc, corr) {
    cc <- c(-Inf, cc, Inf)
    rc <- c(-Inf, rc, Inf)
    pm <- sapply(1:(length(cc)-1), function(c) {
            sapply(1:(length(rc)-1), function(r) {
              pmvnorm(lower=c(cc[c], rc[r]),
                      upper=c(cc[c+1], rc[r+1]),
                      mean=c(0,0),
                      corr=matrix(c(1,corr,corr,1), nrow=2, ncol=2, byrow=TRUE))
          })
    })
    suppressWarnings(lpm <- log(pm))
    lpm[is.nan(lpm)] <- 0
    lpm[!is.finite(lpm)] <- log(.Machine$double.xmin)
    sum(xytab * lpm)
  }

  optf_all <- function(par, xytab) {
    c1 <- ncol(xytab)-1
    c2 <- c1 + nrow(xytab)-1
    lnl(xytab, cc=fscale_cuts(par[1:c1]), rc=fscale_cuts(par[(c1+1):c2]), corr=fscale_corr(par[length(par)] ))
  }

  optf_cut1 <- function(par, parm1, xytab) {
    lnl(xytab, cc=fscale_cuts(par), rc=fscale_cuts(parm1[1:(length(parm1)-1)]), corr=fscale_corr(parm1[length(parm1)] ))
  }

  optf_cut2 <- function(par, parm1, xytab) {
    lnl(xytab, cc=fscale_cuts(parm1[1:(length(parm1)-1)]), rc=fscale_cuts(par), corr=fscale_corr(parm1[length(parm1)] ))
  }

  optf_corr <- function(par, parm1, xytab) {
    c1 <- ncol(xytab)-1
    lnl(xytab, cc=fscale_cuts(parm1[1:c1]), rc=fscale_cuts(parm1[(c1+1):length(parm1)]), corr=fscale_corr(par))
  }

  fscale_cuts <- function(par) {
    cumsum(c(par[1],exp(par[-1])))
  }

  fscale_corr <- function(par) {
    tanh(par)
  }

  weightedTable <- function(x,y,w=rep(1,length(x))) {
    tab <- table(x,y)
    for(i in 1:nrow(tab)) {
      for(j in 1:ncol(tab)) {
        tab[i,j] <- sum(w[ x==dimnames(tab)[[1]][i] & y == dimnames(tab)[[2]][j] ])
      }
    }
    tab
  }
  xytab <- weightedTable(x,y,w)
  #op <- optim(par=c(log(1:(ncol(xytab)-1)), log(1:(nrow(xytab)-1)),cor(x,y)), optf_all, xytab=xytab, control=list(fnscale=-1), method="BFGS")
  #fscale_corr(op$par[length(op$par)])
  cut1 <- log(1:(ncol(xytab)-1))
  cut2 <- log(1:(nrow(xytab)-1))
  cor0 <- cor(as.numeric(x),as.numeric(y))
  dx <- Inf
  cor0m1 <- -1*sign(cor0)*Inf # a long way for cor0
  niter <- 1
  maxit <- 20*max(c(length(cut1), length(cut2)))
  while(dx>1E-6 & niter < 100) {
  	# optimize cut1 alone
  	if(length(cut1) > 1) { # 2 or more values in cut1, use optim
      op <- optim(par=cut1, optf_cut1, parm1=c(cut2, cor0), xytab=xytab, control=list(fnscale=-1, maxit=maxit))
      cut1 <- op$par
    } else { # 1 value in cut1, use optimize
      op <- optimize(optf_cut1, interval=c(cut1-2, cut1+2), parm1=c(cut2, cor0), xytab=xytab, maximum=TRUE)
      cut1 <- op$maximum
    }
    # optimize cut2 alone
  	if(length(cut2) > 1) { # 2 or more values in cut2, use optim
      op <- optim(par=cut2, optf_cut2, parm1=c(cut1, cor0), xytab=xytab, control=list(fnscale=-1, maxit=maxit))
      cut2 <- op$par
    } else { # 1 value in cut2, use optimize
      op <- optimize(optf_cut2, interval=c(cut2-10, cut2+10), parm1=c(cut1, cor0), xytab=xytab, maximum=TRUE)
      cut2 <- op$maximum
    }
    # optimize cor alone
    # cor is always a single parameter
    op <- optimize(optf_corr, interval=c(cor0-2, cor0+2), parm1=c(cut1, cut2), xytab=xytab, maximum=TRUE)
    dx2 <- abs(fscale_corr(cor0m1) - fscale_corr(op$maximum))
    cor0m1 <- cor0
    dx <- abs(fscale_corr(cor0) - fscale_corr(op$maximum))
    if(abs(dx2) < abs(dx)) { # looks like a loop, skip some itterations
      niter <- niter * 2
      cor0 <- mean(c(cor0, op$maximum)) # split the difference
      maxit <- 20*max(c(length(cut1), length(cut2)))
    } else {
      niter <- niter + 1
      maxit <- 10
    }
    cor0 <- op$maximum
    if(verbose) {
      cat("cor0=",cor0, " (",fscale_corr(cor0),") cut1=",cut1," cut2=",cut2," niter=",niter," dx=",dx," dx2=", dx2 ," \n")
    }
    niter <- niter + 1
  }
  fscale_corr(cor0)
}
