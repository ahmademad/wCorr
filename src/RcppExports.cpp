// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cont
double cont(arma::vec x, arma::vec y, const arma::vec& w);
RcppExport SEXP wCorr_cont(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(cont(x, y, w));
    return rcpp_result_gen;
END_RCPP
}
// wrankFast
arma::vec wrankFast(arma::vec x, const arma::vec& w);
RcppExport SEXP wCorr_wrankFast(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(wrankFast(x, w));
    return rcpp_result_gen;
END_RCPP
}
// fixxFast
const arma::vec fixxFast(const arma::vec x, const arma::vec w);
RcppExport SEXP wCorr_fixxFast(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(fixxFast(x, w));
    return rcpp_result_gen;
END_RCPP
}
// mapThetaFast
const arma::vec mapThetaFast(const arma::vec& v);
RcppExport SEXP wCorr_mapThetaFast(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(mapThetaFast(v));
    return rcpp_result_gen;
END_RCPP
}
// optFcFast
double optFcFast(const arma::vec& par, const arma::vec& x, arma::vec w, const arma::vec& M, double temp3, const arma::vec& theta0);
RcppExport SEXP wCorr_optFcFast(SEXP parSEXP, SEXP xSEXP, SEXP wSEXP, SEXP MSEXP, SEXP temp3SEXP, SEXP theta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type temp3(temp3SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    rcpp_result_gen = Rcpp::wrap(optFcFast(par, x, w, M, temp3, theta0));
    return rcpp_result_gen;
END_RCPP
}
// optFFast
double optFFast(const arma::vec& par, const arma::vec& x, const arma::vec w, const arma::vec& M, double temp3);
RcppExport SEXP wCorr_optFFast(SEXP parSEXP, SEXP xSEXP, SEXP wSEXP, SEXP MSEXP, SEXP temp3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type temp3(temp3SEXP);
    rcpp_result_gen = Rcpp::wrap(optFFast(par, x, w, M, temp3));
    return rcpp_result_gen;
END_RCPP
}
// theta
const arma::vec theta(const arma::vec& M);
RcppExport SEXP wCorr_theta(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(theta(M));
    return rcpp_result_gen;
END_RCPP
}
// imapThetaFast2
const arma::vec imapThetaFast2(const arma::vec& theta0);
RcppExport SEXP wCorr_imapThetaFast2(SEXP theta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    rcpp_result_gen = Rcpp::wrap(imapThetaFast2(theta0));
    return rcpp_result_gen;
END_RCPP
}
// mainF
arma::field<arma::vec> mainF(const arma::vec& x, const arma::vec& M, arma::vec w, const arma::vec& theta0);
RcppExport SEXP wCorr_mainF(SEXP xSEXP, SEXP MSEXP, SEXP wSEXP, SEXP theta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    rcpp_result_gen = Rcpp::wrap(mainF(x, M, w, theta0));
    return rcpp_result_gen;
END_RCPP
}
// imapThetaFast
const arma::vec imapThetaFast(const arma::vec& theta0);
RcppExport SEXP wCorr_imapThetaFast(SEXP theta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta0(theta0SEXP);
    rcpp_result_gen = Rcpp::wrap(imapThetaFast(theta0));
    return rcpp_result_gen;
END_RCPP
}
// fscale_cutsFast
const arma::vec fscale_cutsFast(const arma::vec& par);
RcppExport SEXP wCorr_fscale_cutsFast(SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(fscale_cutsFast(par));
    return rcpp_result_gen;
END_RCPP
}
// tableFast
const arma::mat tableFast(const arma::vec& x, const arma::vec& y, const arma::vec& w);
RcppExport SEXP wCorr_tableFast(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(tableFast(x, y, w));
    return rcpp_result_gen;
END_RCPP
}
// discord
int discord(const arma::mat& xytab);
RcppExport SEXP wCorr_discord(SEXP xytabSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xytab(xytabSEXP);
    rcpp_result_gen = Rcpp::wrap(discord(xytab));
    return rcpp_result_gen;
END_RCPP
}
// lnlFast
double lnlFast(const arma::mat& xytab, const arma::mat& pm);
RcppExport SEXP wCorr_lnlFast(SEXP xytabSEXP, SEXP pmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xytab(xytabSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pm(pmSEXP);
    rcpp_result_gen = Rcpp::wrap(lnlFast(xytab, pm));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP wCorr_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP wCorr_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP wCorr_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP wCorr_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
