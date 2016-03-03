if(FALSE) {
  install.packages("RcppArmadillo")
  install.packages("Q:/Paul/wCorr/wCorr_1.3.tar.gz", repos=NULL, type="source")
}

require(wCorr)
require(vcdExtra)
require(doBy)
require(lattice)
setwd("Q:/Paul/wCorr/wCorr/vignettes")
ns <- 10^c(3,4)#,7)
cor0 <- c(-0.99, seq(-0.95, 0.95, by=0.05), 0.99)
reps <- 5

df <- data.frame(i=1:(length(ns)*reps*length(cor0)))

df$cor <- df$cor0 <- df$n <- NA
df$Q <- df$M <- NA

df$pc <- df$pct <- NA
df$ps <- df$pst <- NA

dim(df)

ii <- 1
for(n in ns) {
  for(cori in cor0) {
    if(interactive()) {
      cat("n=",n,"cori=",cori,"\n")
    }
    for(j in 1:reps) {
      df$cor0[ii] <- cori
      df$n[ii] <- n
      
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
      
      st0 <- system.time(fcorp <- weightedCorr(x,M, method = c("Polyserial"), weights=rep(1,length(x)), fast=TRUE, ML=FALSE))
      df$pst[ii] <- st0[3]
      df$ps[ii] <- fcorp
      
      st1 <- system.time(fpolyc <- weightedCorr(M, Q, method = c("Polychoric"), weights=rep(1,length(x)), fast=TRUE, ML=FALSE))
      df$pct[ii] <- st1[3]
      df$pc[ii] <- fpolyc
      
      ii <- ii + 1
    }
  }
}
df$pcb <- df$pc - df$cor0
df$psb <- df$ps - df$cor0

biasdf <- df

save(biasdf, file="bias.RData")
