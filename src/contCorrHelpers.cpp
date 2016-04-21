#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rcpp.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double cont(arma::vec x,  arma::vec y, const arma::vec& w) {
  double sumw = arma::sum(w);
  double xb = arma::sum(w%x)/sumw;
  double yb = arma::sum(w%y)/sumw;
  const arma::vec temp1 = x-xb;
  const arma::vec temp2 = y-yb;
  double numerator = arma::sum(w%temp1%temp2);
  double denominator = pow(arma::sum(w%arma::square(temp1)) * sum(w%arma::square(temp2)), 0.5);
  return numerator/denominator;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec wrankFast(arma::vec x, const arma::vec& w) {
  int size = x.size();
  arma::vec ranks(size);
  for(int i =0; i < size; i++ ) {
    const arma::uvec ind = arma::find(x<x[i]);
    double t1 = arma::sum(w.elem(ind));
    const arma::uvec ind1 = arma::find(x==x[i]);
    arma::vec t2 = w.elem(ind1);
    ranks[i] = t1 + arma::mean(t2) + (arma::sum(t2) - arma::mean(t2))/2;
  }
  return(ranks);
}