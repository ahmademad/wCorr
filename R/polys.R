polys <- function(M, U, w) {
  # M is discrete, U is continuous.
  #Find unique values of x, and sort them.
  uqx <- sort(unique(M))
  N <- length(U)
  n <- c()
  n[1] <- 0
  t_ <- length(uqx)
  for(i in 1:t_) {
    n[i+1] <- n[i] + sum(M==uqx[i])
  }
  theta0 <- mean(U)
  theta <- c(sd(U)* sqrt((N-1)/N), 0, qnorm(n[2:(length(n)-1)]/N))
  #Okay, R stupidity of not having 0 indexed. So separating out theta0
  r <- cor(M, U) 
  s <- sqrt((1/N) * sum((M-mean(M))^2))
  # Setting the value of theta 2 star using the formula
  #theta[2] <- s * r / sum(dnorm(theta[3:length(theta)]))
  thetaStar <- s * r / sum(dnorm(theta[3:length(theta)]))
  theta[t_+2] <- 0 
  R <- sqrt(1-theta[2]^2)
  x <- (U-theta0)/theta[1]
  
  
  returnJ <- function(theta) {
    R <- sqrt(1-theta[2]^2)
    I <- diag(0, nrow=t_+2, ncol=t_+2)
    gamma0 <- c(rep(0, length(theta)))
    gamma1 <- c(rep(0, length(theta)))
    gamma2 <- c(rep(0, length(theta)))
    xi0 <- c(rep(0, length(theta)))
    xi1 <- c(rep(0, length(theta)))
    xi2 <- c(rep(0, length(theta)))
    
    for(s in (3:(t_+1))) {
      gamma0[s] <- integrate(gammaToIntegrate(s,0, theta), lower=-6, upper=6)$value
      gamma1[s] <- integrate(gammaToIntegrate(s,1, theta), lower=-6, upper=6)$value
      gamma2[s] <- integrate(gammaToIntegrate(s,2, theta), lower=-6, upper=6)$value
      
      
      xi0[s] <- integrate(xiToIntegrate(s,0, theta), lower=-6, upper=6)$value
      xi1[s] <- integrate(xiToIntegrate(s,1, theta), lower=--6, upper=6)$value
      xi2[s] <- integrate(xiToIntegrate(s,2, theta), lower=-6, upper=6)$value
      I[s+1,s+1] <- (N/(R^2)) * gamma0[s]
      
      I[1,s+1] <- (N*theta[2]/(theta[1]*R^2))*(gamma0[s] + xi0[s] + xi0[s-1]) # I[0,s]
      I[s+1,1] <- I[1,s+1] #I[s,0]
      
      I[2,s+1] <- (N*theta[2]/(theta[1]*R^2))*(gamma1[s] + xi1[s] + xi1[s-1]) #I[1,s]
      I[s+1,2] <- I[2,s+1] # I[s,1]
      
      I[3,s+1] <- (N*theta[2]/R^4) * (theta[s]*gamma0[s] + theta[s+1]*xi0[s] + xi0[s-1]*theta[s-1]) - (N/R^4)*(gamma1[s] + xi1[s] + xi1[s-1]) #I[2,s]
      I[s+1, 3] <- I[3,s+1]
      
      if(s<=t_) {
        I[s+1, s+2] <- (N/(R^2)) * xi0[s]
        I[s+2, s+1] <- I[s+1, s+2] 
      }
      
    }
    I[1,1] <- (N/theta[1]^2) + (theta[2]/theta[1]) * sum(I[1,(3:(t_+2))])
    I[1,2] <- (N/theta[1]^2) + (theta[2]/theta[1]) * sum(I[2,(3:(t_+2))])
    I[2,1] <- I[1,2]
    I[1,3] <- (N/theta[1]^2) + (theta[2]/theta[1]) * sum(I[3,(3:(t_+2))])
    I[3,1] <- I[1,3]
    
    I[2,2] <- 2*N/(theta[1]^2) + ((N * theta[2]^2)/(theta[1]^2 *R^2)) *(sum(gamma2[3:(t_+1)]) + 2*sum(xi2[3:(t_+1)]))
    
    temp <- (1/R^2) * sum(theta[3:(t_+1)] * I[2,(4:(t_+2))]) - (N)/(theta[1]*R^2)*(sum(gamma2[3:(t_+1)]) + 2*sum(xi2[3:(t_+1)]))
    I[2,3] <- theta[2]* temp
    I[3,2] <- I[2,3]
    
    I[3,3] <- (theta[2]/R^2) * sum(theta[3:(t_+1)] * I[3,(4:(t_+2))]) - (theta[1]/(R^2))*temp
    
    I
  }
  
  Phi <- function(x, s, theta) {
    R <- sqrt(1-theta[2]^2)
    if (s==1) {
      return(rep(NA,length(x))) 
    }
    if(s==2) {
      return(rep(0,length(x))) 
    }
    if(s==t_+2) {
      return(rep(1, length(x)))
    }
    pnorm((theta[s]-theta[2]*x)/R)
  }
  
  phi <- function(x, s, theta) {
    R <- sqrt(1-theta[2]^2)
    if (s==1) {
      return(rep(NA,length(x))) 
    }
    if(s==2 | s==t_+2) {
      return(rep(0,length(x))) 
    }
    
    dnorm((theta[s]-theta[2]*x)/R)
  }
  
  gammas <- function(x, s, theta) {
    denom1 <- Phi(x,s+1, theta)- Phi(x,s, theta) 
    denom2 <- Phi(x,s, theta) - Phi(x,s-1, theta)
    denom1[denom1==0] <- 1 # this is a hack. When the denominator is 0 the numerator will be too (I hope)
    denom2[denom2==0] <- 1
    phi(x,s, theta)^2*(1/denom1 + 1/denom2)
    
  }
  
  xis <- function(x, s, theta) {
    if (s < 2 | s > t_ + 1) {stop("ERROR ERROR")}
    if(s==2 | s==t_+1) {
      return(rep(0,length(x))) 
    }
    denom <- (Phi(x, s+1, theta)-Phi(x, s, theta))
    denom[denom==0] <- 1 # this is a hack. When the denominator is 0 the numerator will be too (I hope)
    -1*phi(x, s, theta)*phi(x, s+1, theta)/denom
  }
  
  
  gammaToIntegrate <- function(s,i, theta) {
    if(i==0) {
      function(x) {
        gammas(x, s, theta) * dnorm(x)
      }
    } else {
      function(x) {
        x^i * gammas(x, s, theta) * dnorm(x)
      }
    }
  }
  
  xiToIntegrate <- function(s,i, theta) {
    function(x) {
      res <- x^i * xis(x, s, theta) * dnorm(x)
      res
    }
  }
  
  alpha <- function(y, theta) {
    function(k, i){
      denom <- Phi(y[i], k+2, theta) - Phi(y[i], k+1, theta)
      #if(denom == 0) {
      #  denom <- 1
      #}
      (phi(y[i], k+2, theta) - phi(y[i], k+1, theta)) / denom
    }
    
  }
  
  beta <- function(y, theta) {
    function(k,i) {
      denom <- (Phi(y[i], k+2, theta) - Phi(y[i], k+1, theta))
      if(denom == 0) {
        denom=1
      }
      ((theta[k+2] * phi(y[i], k+2, theta)) - (theta[k+1] * phi(y[i], k+1, theta))) / denom
    }
  }
  
  deltas <- function(x, theta) {
    R <- sqrt(1-theta[2]^2)
    delta <- c()
    delta0 <- (theta[2]/(theta[1]*R)) * doubleSum(alpha(x, theta)) + (1/theta[1]) * sum(x)
    delta1 <- (theta[2]/(theta[1]*R)) *doubleSum(alpha(x, theta), x) + (1/theta[1]) * sum(x*x - 1)
    
    delta2 <- (theta[2]/(R^3)) * doubleSum(beta(x, theta)) - (1/R^3) * doubleSum(alpha(x, theta), x)
    
    for (s in seq(3,(t_+1))) {
      temp <- (1/R) * (sum(phi(x[(n[s-2]+1):(n[s-1])], s, theta)/ (Phi(x[(n[s-2]+1):(n[s-1])], s, theta) - Phi(x[(n[s-2]+1):(n[s-1])], s-1, theta))) - 
                         sum(phi(x[(n[s-1]+1):(n[s])], s, theta)/ (Phi(x[(n[s-1]+1):(n[s])], s+1, theta) - Phi(x[(n[s-1]+1):(n[s])], s, theta))))
      delta <- c(delta, temp)
    }
    c(delta2, delta)  
  }
  
  doubleSum <- function(func, y=rep(1, length(x))) {
    tot <- 0
    for(k in seq(1,t_)) {
      for(i in seq((n[k] + 1),n[k+1] )) {
        tot <- tot + y[i] * func(k, i)
      }
    }
    tot
  }


g <- c(theta[2:(t_+1)])
d <- deltas(x, theta)
I <- returnJ(theta)
J <- I[-c(1,2), -c(1,2)]
temp <- qr.solve(qr(J), d)
i <- 0
while(TRUE) {
  tempnew <- qr.solve(qr(J), d)
  if(sum(abs(temp - tempnew)) < 0.0001 & i != 0) {
    break
  }
  temp <- tempnew  
  i <- i+1
  g <- c(theta[2:(t_+1)]) + temp
    theta <- c(theta[1], g, 0)
  I <- returnJ(theta)
  J <- I[-c(1,2), -c(1,2)]
  d <- deltas(x, theta)
}

#temp <- sqrt(solve(tail(I, n=1)[[1]])[3,3])

theta[2]
}


