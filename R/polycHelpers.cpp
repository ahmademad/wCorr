#include <RcppArmadillo.h>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
const arma::vec imapThetaFast(const arma::vec& theta0) {
  int n = theta0.size();
  arma::vec temp(n);
  temp(0) = theta0(0);
  temp.subvec(1, n-1) = arma::log(theta0.subvec(1, n-1) - theta0.subvec(0, n-2));
  return temp;
}

// [[Rcpp::export]]
const arma::vec fscale_cutsFast(const arma::vec& par) {
  arma::vec temp(par.size());
  temp(0) = par(0);
  temp.subvec(1,par.size()-1) = arma::cumsum(arma::exp(par.subvec(1, par.size()-1)))+par(0);
  return temp;
}

// [[Rcpp::export]]
const arma::mat tableFast(const arma::vec& x, const arma::vec& y) {
  const arma::vec xUni = arma::sort(arma::unique(x)); 
  const arma::vec yUni = arma::sort(arma::unique(y)); 
  
  arma::mat tab(xUni.size(), yUni.size());
  for (int i = 0 ; i < xUni.size(); i+=1) {
    const arma::uvec x_ind = arma::find(x==xUni(i));
    for(int j = 0; j < yUni.size(); j+=1) {
      const arma::vec y_val = y.elem(arma::conv_to<arma::uvec>::from(x_ind));
      const arma::vec temp = arma::conv_to<arma::vec>::from(arma::find(y_val==yUni(j)));
      int n = temp.size();
      tab(i,j) = n;
      
    }
  }
  for (int i = 0; i < xUni.size(); i+=1) {
    for(int j = 0; j < yUni.size(); j+=1) 
      tab(i,j) = 
  }
  return tab;
}


weightedTable <- function(x,y,w=rep(1,length(x))) {
  tab <- tableFast(x,y)
  t1 <- sort(unique(x))
  t2 <- sort(unique(y))
  for(i in 1:nrow(tab)) {
    for(j in 1:ncol(tab)) {
      tab[i,j] <- sum(w[ x==t1[i] & y == t2[j] ])
    }
  }
  tab
}

