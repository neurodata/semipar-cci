// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// random_factor_decomp
SEXP random_factor_decomp(SEXP A_r, int n, int m, SEXP g_r, int rank);
RcppExport SEXP RandomFactorTensor_random_factor_decomp(SEXP A_rSEXP, SEXP nSEXP, SEXP mSEXP, SEXP g_rSEXP, SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A_r(A_rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< SEXP >::type g_r(g_rSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    rcpp_result_gen = Rcpp::wrap(random_factor_decomp(A_r, n, m, g_r, rank));
    return rcpp_result_gen;
END_RCPP
}