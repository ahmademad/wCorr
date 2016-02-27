
polys_opt1Fast <- function(x, M, w, ML=FALSE) {
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
    bob <- bobyqa(par=c(imapCorFast(cor(x,M)),imapTheta(theta0)), fn=optF(x,M,w))
    return(  tanh(bob$par[1]))
  } else {
    w <- w/sum(w)
    fx = fixxFast(x, w);
    ftheta0 = mapThetaFast(imapTheta(theta0));
    temp1 <- ftheta0[M+2]
    temp2 <- ftheta0[M+1]
    temp3 <- w*dnorm(x)
    
    opt <- optimize(optFcFast, interval = imapCorFast(cor(x,M)) + c(-3,3), x,w,temp1, temp2, temp3)
    return( tanh(opt$minimum) )
  }
}


