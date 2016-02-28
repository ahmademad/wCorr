// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
const arma::vec fixxFast(const arma::vec x, const arma::vec w) {
  const arma::vec temp = x-arma::sum(x%w);
  double sdx = arma::sum(arma::square(temp)%w);
  return (temp)/pow(sdx, 0.5);
}


// [[Rcpp::export]]

const arma::vec mapThetaFast(const arma::vec& v) {
  int s = v.size();
  arma::vec temp(s+3);
  temp(0) = R_NaReal;
  temp(1) = -std::numeric_limits<double>::infinity();
  temp(s+2) = std::numeric_limits<double>::infinity();
  temp.subvec(3,s+1) = arma::cumsum(arma::exp(v.tail(s-1)))+v(0);
  temp(2) = v(0);
  return temp;
}
// [[Rcpp::export]]
double optFcFast(const arma::vec& par, const arma::vec& x, const arma::vec& w, 
                 const arma::vec& temp1, const arma::vec& temp2, const arma::vec& temp3) {
 
  double rho = tanh(par[0]);
  double R = pow((1-rho*rho), 0.5);

  const arma::vec Qp2 = (temp1 - (rho*x))/R;
  const arma::vec Qp1 = (temp2 - (rho*x))/R;
  const arma::vec t = (Rcpp::pnorm(Rcpp::NumericVector(Qp2.begin(),Qp2.end())) -
                       Rcpp::pnorm(Rcpp::NumericVector(Qp1.begin(),Qp1.end())));
  double res= temp3(0) + arma::sum(w%arma::log(t));
  if (res == -std::numeric_limits<double>::infinity()) 
    res = -std::numeric_limits<double>::max();
  else
    res = -res;
  return res;
  
}

// [[Rcpp::export]]

const arma::vec theta(const arma::vec& uM, const arma::vec& M) {
  arma::vec theta0(uM.size()-1);
  double s = 0;
  for(int i=0; i < uM.size()-1; i+=1) {
    s = arma::mean(arma::conv_to<arma::vec>::from(M<=uM(i)));
    theta0(i) =  R::qnorm(s, 0.0, 1.0, 1, 0);
   
  }
  return theta0;
}

// [[Rcpp::export]]
const arma::vec imapThetaFast(const arma::vec& theta0) {
  int n = theta0.size();
  arma::vec temp(n);
  temp(0) = theta0(0);
  temp.subvec(1, n-1) = arma::log(theta0.subvec(1, n-1) - theta0.subvec(0, n-2));
  return temp;
}
// [[Rcpp::export]]

arma::field<arma::vec> mainF(const arma::vec x, const arma::vec M, arma::vec  w) {
  const arma::vec uM = arma::sort(arma::unique(M));
  const arma::vec theta0 = theta(uM, M);
  

    const arma::vec ftheta0 = mapThetaFast(imapThetaFast(theta0));
    const arma::vec temp1 = ftheta0.elem(arma::conv_to<arma::uvec>::from(M)+1); 
    const arma::vec temp2 = ftheta0.elem(arma::conv_to<arma::uvec>::from(M)); 
    w = w/arma::sum(w);
    const arma::vec fx= fixxFast(x, w);
    
    arma::vec temp3 = Rcpp::dnorm(Rcpp::NumericVector(fx.begin(),fx.end()), 0.0, 1.0, true);
    arma::vec temp(2);
    const arma::mat s = arma::cor(x,M);
    double c = s(0,0); 
    double s1 = arma::sum(w%temp3);
    
    arma::vec t(1);
    t.fill(s1);
    temp(0) = c-3;
    temp(1) = c+3;
    arma::field<arma::vec> F(5);
    F(0) = temp;
    F(1) = temp1;
    F(2) = temp2;
    F(3) = t;
    F(4) = fx;
    //double r  = R.parseEval("optimize(optFcFast, interval = temp, x,w,temp1, temp2, temp3)");
    return F;
  
}

