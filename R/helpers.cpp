#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


double square(double x) {
  return x*x;
}
// [[Rcpp::export]]
arma::vec fixxFast(vec x, vec w) {
  omp_set_num_threads(4);
  arma::vec mux(x.size());
  double s = 0;
  for(int i=0; i < x.size(); i+=1)
    s+=x(i)*w(i);
  
  for(int i = 0 ; i < x.size(); i+=1) {
    mux(i) = s;
  }
 arma::vec temp = x-mux;
 arma::vec temp2(x.size());
 for(int i=0; i < x.size(); i+=1)
   temp2(i) = pow(temp(i)*w(i),2);

 double sdx = arma::sum(temp2);
 return (temp)/pow(sdx, 0.5);
}

// [[Rcpp::export]]
double imapCorFast(double cor) {
  return atanh(cor);
}

// [[Rcpp::export]]

NumericVector mapThetaFast(NumericVector v) {
  NumericVector temp = exp(v);
  temp[0] = v[0];
  NumericVector vv = cumsum(temp);
  NumericVector temp2(vv.size()+3);
  temp2[0] = NAN;
  temp2[1] = -std::numeric_limits<double>::infinity();;
  int i = 1;
  for(i =0; i < vv.size(); i+=1)
    temp2[i+2] = vv[i];
  temp2[vv.size()+2] = std::numeric_limits<double>::infinity();
  return wrap(temp2);
}

// [[Rcpp::export]]

double optFcFast(NumericVector par, NumericVector x, NumericVector w, NumericVector temp1, NumericVector temp2,
                 NumericVector temp3) {
 
  double rho = tanh(par[0]);
  double R = pow((1-rho*rho), 0.5);

  NumericVector Qp2 = (temp1 - (rho*x))/R;
  NumericVector Qp1 = (temp2 - (rho*x))/R;
  double res= sum(temp3) + sum(w* (pnorm(Qp2) - pnorm(Qp1)));
  return -res;
  
}



// 
// polysLnL <- function(x,M,rho,theta,w) {
//   R <- (1-rho^2)^0.5
//   Qp2 <- (theta[M+2] - rho*x) / R
//   Qp1 <- (theta[M+1] - rho*x) / R
// #-log(R) + sum(w*dnorm(x,log=TRUE)) + sum(w*log(pnorm(Qp2) - pnorm(Qp1)))
// #-log(R) + sum(w*dnorm(x,log=TRUE)) + sum(w*log(Phi(x,M+2,rho,theta) - Phi(x,M+1,rho,theta)))
//   sum(w*dnorm(x,log=TRUE)) + sum(w* log(pnorm(Qp2) - pnorm(Qp1)))
// }
