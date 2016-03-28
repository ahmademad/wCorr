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
  everusew <- sum(usew)>0
  
  ns <- n
  cor0 <- rho
  df <- data.frame(n=n,rho=rho, ML=ML, fast=fast, reset=reset, usew=usew)
  

  df$spear <- df$speart <- NA
  df$Q <- df$M <- NA
  df$pear <- df$peart <- NA
  df$pc <- df$pct <- NA
  df$ps <- df$pst <- NA
  
  ii <- 1
  while(ii <= nrow(df)) {
    cori <- df$rho[ii]
    n <- df$n[ii]
    ML <- df$ML[ii]
    fast <- df$fast[ii]
    reset <- df$reset[ii]
    usew <- df$usew[ii]

    if(interactive()) {
      cat("n=",n,"cori=",cori,"\n")
    }
    if(reset) {
      M <- 1  
      Q <- 1
      nm <- sample(2:5,1)
      tm <- sort(rnorm(nm))
      nq <- sample(2:5,1)
      tq <- sort(rnorm(nq))
      n <- ifelse(everusew, 5*df$n[ii], df$n[ii])
      while( (length(unique(M)) < 2) | (length(unique(Q)) < 2) ) {
        theta1 <- c(NA,-Inf,tq,Inf)
        theta2 <- c(NA,-Inf,tm,Inf)
        cr <- cori
        x <- y <- w <- M <- Q <- c()
        while(length(w) < df$n[ii]) {
          xp <- rnorm(n)
          Qp <- rep(NA,n)
          for(i in 2:length(theta1)) {
            Qp <- ifelse(xp>theta1[i], i, Qp)
          }
          Qp <- Qp - 1
          
          yp <- sqrt(1-cr^2)*rnorm(n) + cr*xp
          Mp <- rep(NA,n)
          for(i in 2:length(theta2)) {
            Mp <- ifelse(yp>theta2[i], i, Mp)
          }
          Mp <- Mp - 1
  
          Mp <- as.numeric(as.factor(Mp))
          if(everusew) {
            wp <- (xp-yp)^2+1
  
            pr <- 1/wp
            pr <- pr/(sum(pr) * 3)
            wp <- 1/pr
            #samp <- sample(1:n, size=df$n[ii], replace=FALSE, prob=pr)
            samp <- (1:n)[runif(n)<pr]
            M <- c(M,Mp[samp])
            x <- c(x,xp[samp])
            Q <- c(Q,Qp[samp])
            y <- c(y,yp[samp])
            w <- c(w,wp[samp])
          } else {
            M <- Mp
            x <- xp
            Q <- Qp
            y <- yp
            w <- rep(1/n, n)
          }
        }
        M <- M[1:df$n[ii]]
        x <- x[1:df$n[ii]]
        Q <- Q[1:df$n[ii]]
        y <- y[1:df$n[ii]]
        w <- w[1:df$n[ii]]
      }
      df$M[ii] <- length(unique(M))
      df$Q[ii] <- length(unique(Q))
    } else {
      df$cor[ii] <- df$cor[ii-1]
      df$M[ii] <- length(unique(M))
      df$Q[ii] <- length(unique(Q))
    }
    
    if(usew) {
      wu <- w
    } else {
      wu <- rep(1,length(x))
    }
    #cat("usew=",usew,"ii=",ii,"\n")
    #save(x,wu,y,Q,M,file="tmp.RData")
    #print(wu)
    
    
    st0 <- system.time(fcorp <- weightedCorr(x,y, method="Pearson", weights=wu, fast=fast, ML=ML))
    df$peart[ii] <- st0[3]
    df$pear[ii] <- fcorp

    st0 <- system.time(fcorp <- weightedCorr(x,y, method="Spearman", weights=wu, fast=fast, ML=ML))
    df$speart[ii] <- st0[3]
    df$spear[ii] <- fcorp
    
    st0 <- system.time(fcorp <- weightedCorr(x,M, method="Polyserial", weights=wu, fast=fast, ML=ML))
    df$pst[ii] <- st0[3]
    df$ps[ii] <- fcorp
    
    st1 <- system.time(fpolyc <- weightedCorr(M, Q, method="Polychoric", weights=wu, fast=fast, ML=ML))
    df$pct[ii] <- st1[3]
    df$pc[ii] <- fpolyc
    
    ii <- ii + 1
  }
  dfout <- data.frame(n=rep(df$n,4),
                      rho=rep(df$rho,4),
                      ML=rep(df$ML,4),
                      usew=rep(df$usew,4),
                      fast=rep(df$fast,4),
                      est=c(df$pear, df$spear, df$ps, df$pc),
                      t=c(df$peart, df$speart, df$pst, df$pct),
                      type=rep(c("Pearson", "Spearman", "Polyserial", "Polychoric"),each=nrow(df)))
  dfout
}



