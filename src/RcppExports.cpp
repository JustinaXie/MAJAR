// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logsumexp
double logsumexp(NumericVector x);
RcppExport SEXP _MAJAR_logsumexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp(x));
    return rcpp_result_gen;
END_RCPP
}
// llbR0_j
double llbR0_j(NumericVector betajk_j_G, NumericVector betajk_j_GT, NumericVector sjk2_j_G, NumericVector sjk2_j_GT, double lambda, double alpha);
RcppExport SEXP _MAJAR_llbR0_j(SEXP betajk_j_GSEXP, SEXP betajk_j_GTSEXP, SEXP sjk2_j_GSEXP, SEXP sjk2_j_GTSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j_G(betajk_j_GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j_GT(betajk_j_GTSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j_G(sjk2_j_GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j_GT(sjk2_j_GTSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(llbR0_j(betajk_j_G, betajk_j_GT, sjk2_j_G, sjk2_j_GT, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// deltis
NumericVector deltis(NumericVector betajk_j_G, NumericVector betajk_j_GT, NumericVector sjk2_j_G, NumericVector sjk2_j_GT, double lambda, double alpha);
RcppExport SEXP _MAJAR_deltis(SEXP betajk_j_GSEXP, SEXP betajk_j_GTSEXP, SEXP sjk2_j_GSEXP, SEXP sjk2_j_GTSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j_G(betajk_j_GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j_GT(betajk_j_GTSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j_G(sjk2_j_GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j_GT(sjk2_j_GTSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(deltis(betajk_j_G, betajk_j_GT, sjk2_j_G, sjk2_j_GT, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ll_R1_j
double ll_R1_j(double bs2_j_G, double bs2_j_GT, double os22_j_G, double os22_j_GT, double b2s2_j_G, double b2s2_j_GT, double lambda, double alpha, double rho, double tau2, arma::vec sjk2_j_G, arma::vec sjk2_j_GT);
RcppExport SEXP _MAJAR_ll_R1_j(SEXP bs2_j_GSEXP, SEXP bs2_j_GTSEXP, SEXP os22_j_GSEXP, SEXP os22_j_GTSEXP, SEXP b2s2_j_GSEXP, SEXP b2s2_j_GTSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP tau2SEXP, SEXP sjk2_j_GSEXP, SEXP sjk2_j_GTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bs2_j_G(bs2_j_GSEXP);
    Rcpp::traits::input_parameter< double >::type bs2_j_GT(bs2_j_GTSEXP);
    Rcpp::traits::input_parameter< double >::type os22_j_G(os22_j_GSEXP);
    Rcpp::traits::input_parameter< double >::type os22_j_GT(os22_j_GTSEXP);
    Rcpp::traits::input_parameter< double >::type b2s2_j_G(b2s2_j_GSEXP);
    Rcpp::traits::input_parameter< double >::type b2s2_j_GT(b2s2_j_GTSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sjk2_j_G(sjk2_j_GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sjk2_j_GT(sjk2_j_GTSEXP);
    rcpp_result_gen = Rcpp::wrap(ll_R1_j(bs2_j_G, bs2_j_GT, os22_j_G, os22_j_GT, b2s2_j_G, b2s2_j_GT, lambda, alpha, rho, tau2, sjk2_j_G, sjk2_j_GT));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _MAJAR_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _MAJAR_rcpparma_outerproduct(SEXP xSEXP) {
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
RcppExport SEXP _MAJAR_rcpparma_innerproduct(SEXP xSEXP) {
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
RcppExport SEXP _MAJAR_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MAJAR_logsumexp", (DL_FUNC) &_MAJAR_logsumexp, 1},
    {"_MAJAR_llbR0_j", (DL_FUNC) &_MAJAR_llbR0_j, 6},
    {"_MAJAR_deltis", (DL_FUNC) &_MAJAR_deltis, 6},
    {"_MAJAR_ll_R1_j", (DL_FUNC) &_MAJAR_ll_R1_j, 12},
    {"_MAJAR_rcpparma_hello_world", (DL_FUNC) &_MAJAR_rcpparma_hello_world, 0},
    {"_MAJAR_rcpparma_outerproduct", (DL_FUNC) &_MAJAR_rcpparma_outerproduct, 1},
    {"_MAJAR_rcpparma_innerproduct", (DL_FUNC) &_MAJAR_rcpparma_innerproduct, 1},
    {"_MAJAR_rcpparma_bothproducts", (DL_FUNC) &_MAJAR_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MAJAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
