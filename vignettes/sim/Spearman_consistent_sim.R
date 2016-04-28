n <- 1000
i <- 1
biast <- c()
rmset <- c()
nt <- c()
for(n in c(100, 1000, 10000, 100000)) {
for(ind in 1:200) {  
x <- sort(rnorm(n))
y <- rnorm(n)
y <- y[ order(y + rnorm(n)) ]

r <- rank(x)
s <- rank(y)

wi <- rep(5+1:10, each=n/10)
pi <- 1/wi

samp <- (1:n)[runif(n) < pi]

xs <- x[samp]
ys <- y[samp]
wis <- wi[samp]

rhat <- wCorr:::wrank(xs, wis)
shat <- wCorr:::wrank(ys, wis)
#plot(r[samp],rhat)
#plot(r[samp], (r[samp]-rhat) / n)

rbarhat <- sum(wis * rhat) / sum(wis)
rbar <- mean(r)
c(rbar, rbarhat, rbar-rbarhat)

sbarhat <- sum(wis * shat) / sum(wis)
sbar <- mean(s)
c(sbar, sbarhat, sbar-sbarhat)

covhat <- sum(wis * (rhat - rbarhat) * (shat - sbarhat))
cov    <- sum(      (r    - rbar)    * (s    - sbar))
c(cov, covhat, cov-covhat)

#var1hat <- sum(wis* (rhat - rbarhat)^2)
#var1 <- sum( (r - rbar)^2 )
#c(var1, var1hat, var1-var1hat)

shat <- sqrt(sum(wis* (rhat - rbarhat)^2) * sum(wis* (shat - sbarhat)^2))
st <- sqrt(sum((r - rbar)^2) * sum((s - sbar)^2))
c(st, shat, st-shat)

c(cov/st, covhat/shat, cov/st - covhat/shat)

rhohat <- covhat/shat
rhohat <- wCorr::weightedCorr(x=xs,y=ys,method="Spearman",weights = wis, fast=TRUE)

biast[i] <- cov/st - rhohat
rmset[i] <- (cov/st - rhohat)^2
nt[i] <- n
if(i %%100 == 0) {
  agg <- aggregate(rmse ~ n, data.frame(n=nt[1:i],rmse=rmset[1:i]), mean)
  agg$rmse <- sqrt(agg$rmse)
  plot(agg, log="xy")
}
i <- i + 1
}
}

