
# # example:
# n <- 10 # number of units to rank
# n2 <- 4 # the first n2 units are tied, demonstrating how ties are handled
# x <- rnorm(n) # make random data
# if(n2>1) { # implement replication of first n2 units
#   for(i in 2:n2) {
#     x[i] <- x[1]
#   }
# }
# rk <- rank(x, ties.method="average") # traditional rank with average
# rk2 <- wrank(x)
# all.equal(rk, rk2) # the ranks are equal when weights are all 1
# w <- rchisq(n, df=1) # make random ranks with strictly positive values
# rk3 <- wrank(x, w)
# cor(rk, rk3, method="spearman") # ordering is preserved
# cbind(x, rk, rk2, rk3, w)



wrank <- function(x, w=rep(1,length(x))) {
  sapply(1:n, function(i) {
    t1 <- sum(w[x<x[i]]) # ranked below every unit below it
    t2 <- w[x==x[i]] #  ties
    # Note: when selecting the range to average over you have to figure out which unit is first.
    #       this method assumes that all of the units are exchangable and integrates over all
    #       units in the tie.
    # mean(t2) brings the unit up to the smallest rank
    # sum(t2) - mean(t2) /2 is then the middle of the ranks
    t1  + mean(t2) + (sum(t2) -mean(t2))/2
  })
}
