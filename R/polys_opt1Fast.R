
polys_opt1Fast <- function(x, M, w, ML=FALSE) {
    polysLnL <- function(x,M,rho,theta,w) {
      R <- (1-rho^2)^0.5
      Qp2 <- (theta[M+2] - rho*x) / R
      Qp1 <- (theta[M+1] - rho*x) / R
      #-log(R) + sum(w*dnorm(x,log=TRUE)) + sum(w*log(pnorm(Qp2) - pnorm(Qp1)))
      #-log(R) + sum(w*dnorm(x,log=TRUE)) + sum(w*log(Phi(x,M+2,rho,theta) - Phi(x,M+1,rho,theta)))
      sum(w*dnorm(x,log=TRUE)) + sum(w* log(pnorm(Qp2) - pnorm(Qp1)))
    }

    mapTheta <- function(v) {
         vv <- cumsum(c(v[1],exp(v[-1])))
         c(NA,-Inf,vv,Inf)
       }
    
  optF <- function(x,M,w) {
    w <- w/sum(w)
    fx <- fixxFast(x,w)
    function(par) {
      res <- polysLnL(fx,M,tanh(par[1]),mapThetaFast(par[-1]),w)
      ifelse(res==-Inf,.Machine$double.xmax,-1*res)
    }
  }
  
  
  M <- as.numeric(as.factor(M))
  uM <- sort(unique(M))
  theta0 <- sapply(uM[-length(uM)],function(z) qnorm(mean(M<=z)) )
  imapTheta <- function(theta0) {
    c(theta0[1], log(theta0[-1]-theta0[-length(theta0)]))
  }
  #bob <- bobyqa(par=c(imapCor(cor(x,M)),imapTheta(theta0)), fn=optF(x,M,w))
  if(ML) {
    bob <- bobyqa(par=c(imapCor(cor(x,M)),imapTheta(theta0)), fn=optF(x,M,w))
    return(  tanh(bob$par[1]))
  } else {
    opt <- optimize(optFcFast(x,M,w,imapTheta(theta0), c(0)), imapCorFast(cor(x,M)) + c(-3,3))
    return( tanh(opt$minimum) )
  }
}
