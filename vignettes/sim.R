wCorrSim <- function(n, rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=FALSE) {
  len <- max(c(length(n), length(rho), length(ML), length(fast), length(reset), length(usew)))
  vec <- c("n", "rho", "ML", "fast", "reset", "usew")
  for(i in 1:length(vec)) {
    var <- get(vec[i])
    if(length(var) != len) {
      if(length(var) != 1) {
        stop("length of ", sQuote(vec[i]), " must be 1 or the same as the longest vector passed to sim")
      } else {
        var <- rep(var,len)
      }
    }
    assign(vec[i],var)
  }
  
  if(sum(usew)>1) {
    stop("usew not implemented")
  }
  
  ns <- n
  cor0 <- rho
  df <- data.frame(n=n,rho=rho, ML=ML, fast=fast, reset=reset, usew=usew)
  

  df$spear <- df$speart <- NA
  df$pear <- df$peart <- NA
  df$pc <- df$pct <- NA
  df$ps <- df$pst <- NA

  ii <- 1
  for(n in ns) {  
    for(cori in cor0) {
      if(interactive()) {
        cat("n=",n,"cori=",cori,"\n")
      }
      cori <- df$rho[ii]
      n <- df$n[ii]
      ML <- df$ML[ii]
      fast <- df$fast[ii]
      reset <- df$reset[ii]
      usew <- df$usew[ii]
      
      if(reset) {
        M <- 1  
        Q <- 1
        nm <- sample(2:5,1)
        tm <- sort(rnorm(nm))
        nq <- sample(2:5,1)
        tq <- sort(rnorm(nq))
        while( (length(unique(M)) < 2) | (length(unique(Q)) < 2) ) {
          theta <- c(NA,-Inf,tq,Inf)
          x <- rnorm(n)
          cr <- cori
          Q <- rep(NA,n)
          for(i in 2:length(theta)) {
            Q <- ifelse(x>theta[i], i, Q)
          }
          Q <- Q - 1
          
          y <- sqrt(1-cr^2)*rnorm(n) + cr*x
          theta <- c(NA,-Inf,tm,Inf)
          M <- rep(NA,n)
          for(i in 2:length(theta)) {
            M <- ifelse(y>theta[i], i, M)
          }
          M <- M - 1
          theta <- c(NA,-Inf,tq,Inf)
          
          M <- as.numeric(as.factor(M))
        }
        df$cor[ii] <- cor(x,y)
        df$M[ii] <- length(unique(M))
        df$Q[ii] <- length(unique(Q))
      } else {
        df$cor[ii] <- df$cor[ii-1]
        df$M[ii] <- length(unique(M))
        df$Q[ii] <- length(unique(Q))
      }
      
      st0 <- system.time(fcorp <- weightedCorr(x,y, method="Pearson", weights=rep(1,length(x)), fast=fast, ML=ML))
      df$peart[ii] <- st0[3]
      df$pear[ii] <- fcorp

      st0 <- system.time(fcorp <- weightedCorr(x,y, method="Spearman", weights=rep(1,length(x)), fast=fast, ML=ML))
      df$speart[ii] <- st0[3]
      df$spear[ii] <- fcorp
      
      st0 <- system.time(fcorp <- weightedCorr(x,M, method="Polyserial", weights=rep(1,length(x)), fast=fast, ML=ML))
      df$pst[ii] <- st0[3]
      df$ps[ii] <- fcorp
      
      st1 <- system.time(fpolyc <- weightedCorr(M, Q, method="Polychoric", weights=rep(1,length(x)), fast=fast, ML=ML))
      df$pct[ii] <- st1[3]
      df$pc[ii] <- fpolyc
      
      ii <- ii + 1
    }
  }
  dfout <- data.frame(n=rep(df$n,4),
                      rho=rep(df$rho,4),
                      ML=rep(df$ML,4),
                      usew=rep(df$usew,4),
                      fast=rep(df$fast,4),
                      est=c(df$pear, df$spear, df$ps, df$pc),
                      t=c(df$tpear, df$tspear, df$tps, df$tpc),
                      type=rep(c("Pearson", "Spearman", "Polyserial", "Polychoric")))
  dfout
}


if(FALSE) {
  wCorrSim(n=c(10,100,1000), rho=c(0,0,0))
}

