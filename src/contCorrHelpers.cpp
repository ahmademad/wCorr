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
