if(FALSE) {
  install.packages("RcppArmadillo")
  install.packages("Q:/Paul/wCorr/wCorr_1.3.tar.gz", repos=NULL, type="source")
}

require(vcdExtra)
require(doBy)
require(lattice)
setwd("Q:/Paul/wCorr/wCorr/vignettes")
require(wCorr)
ns <- 10^c(2)
#ns <- 20
#cor0 <- rev(c(0.99, 0.98, 0.95, 0.9, 0.8, 0.5, 0.2, 0, -0.99))
cor0 =  -0.99
reps <- 10

df <- data.frame(i=1:(length(ns)*reps*length(cor0)))

df$cor <- df$cor0 <- df$n <- NA
df$Q <- df$M <- NA

df$pst <- df$fpst <- df$dps <- NA
df$pct <- df$fpct <- df$dpc <- NA
df$pc <- NA
df$ps <- NA
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

      st0 <- system.time(corp <- weightedCorr(x,M, method = c("Polyserial"), weights=rep(1,length(x)), fast=FALSE, ML=TRUE))
      df$pst[ii] <- st0[3]
      df$ps[ii] <- corp

      st0 <- system.time(fcorp <- weightedCorr(x,M, method = c("Polyserial"), weights=rep(1,length(x)), fast=TRUE, ML=TRUE))
      df$fpst[ii] <- st0[3]
      df$dps[ii] <- corp - fcorp

      st1 <- system.time(polyc <- weightedCorr(M, Q, method = c("Polychoric"), weights=rep(1,length(x)), fast=FALSE, ML=TRUE))
      df$pct[ii] <- st1[3]
      df$pc[ii] <- polyc

      st1 <- system.time(fpolyc <- weightedCorr(M, Q, method = c("Polychoric"), weights=rep(1,length(x)), fast=TRUE, ML=TRUE))
      df$fpct[ii] <- st1[3]
      df$dpc[ii] <- polyc - fpolyc


      #if( interactive() &  (df$dps[ii]) > 1E-4 ) {
      #  stop("non zero difference in PS")
      #}

      if( interactive() & (df$dpc[ii]) > 0.001 ) {
        stop("non zero difference in PC")
      }
      ii <- ii + 1
    }
  }
}
df$pcb <- df$pc - df$cor0
df$psb <- df$ps - df$cor0

weightedCorr(x,M, method = c("Polyserial"),
             weights=rep(1,length(x)), fast=FALSE, ML=FALSE)

save(df, file="fast.RData")
