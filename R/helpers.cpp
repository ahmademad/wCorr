#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


double square(double x) {
  return x*x;
}
// [[Rcpp::export]]
NumericVector fixxFast(NumericVector x, NumericVector w) {
  NumericVector mux(x.size(), sum(x*w));
  NumericVector temp = w*(x-mux);
  double sdx = sum(sapply(temp, square));
  return wrap((x-mux)/pow(sdx, 0.5));
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
  temp2[1] = -89585;
  int i = 1;
  for(i =0; i < vv.size(); i+=1)
    temp2[i+2] = vv[i];
  temp2[vv.size()+2] = 89898;
  return wrap(temp2);
}

// [[Rcpp::export]]

double temp1(NumericVector x, IntegerVector M, NumericVector w, NumericVector theta0, NumericVector par) {
  w = w/sum(w);
  NumericVector fx = fixxFast(x, w);
  NumericVector ftheta0 = mapThetaFast(theta0);
  double rho = tanh(par[0]);
  double R = pow((1-rho*rho), 0.5);
  IntegerVector idx1 = M +IntegerVector(M.size(), 1);
  IntegerVector idx2 = M +IntegerVector(M.size(), 0);
  
  NumericVector temp1 = ftheta0[idx1];
  NumericVector Qp2 = (temp1 - (rho*x))/R;
  NumericVector temp2 = ftheta0[idx2];
  NumericVector Qp1 = (temp2 - (rho*x))/R;
  NumericVector temp3 =(w*dnorm(x));
  double res= sum(temp3) + sum(w* (pnorm(Qp2) - pnorm(Qp1)));
  cout << res;
  return -res;
  
}

double optFcFast(NumericVector pars) {
  
  temp1
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
