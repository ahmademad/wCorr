setwd("/Users/ahmademad/Documents/wCorr/R")
require(Rcpp)
require(rbenchmark)
U <- c(72, 88, 112, 93, 86, 87, 78, 69, 101, 108, 104, 
       92, 77, 80, 99, 92, 72, 81, 88, 104, 103, 85, 96, 96, 104)
M <- c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2)
w <- rep(1, length(U))
x <- U
M <- as.numeric(as.factor(M))
uM <- sort(unique(M))
theta0 <- sapply(uM[-length(uM)],function(z) qnorm(mean(M<=z)) )
sourceCpp("helpers.cpp")

benchmark(polys_opt1(rep(U, 1000), rep(M, 1000), rep(w, 1000), FALSE),polys_opt1(rep(U, 1000), rep(M, 1000), rep(w, 1000), TRUE),
          polys_opt1Fast(rep(U, 1000), rep(M, 1000), rep(w, 1000), FALSE), polys_opt1Fast(rep(U, 1000), rep(M, 1000), rep(w, 1000), TRUE))
benchmark(theta(rep(theta0,10000)), test2(rep(theta0,10000)))

test2 <- function(uM, M) {
  sapply(uM[-length(uM)],function(z) qnorm(mean(M<=z)) )
}