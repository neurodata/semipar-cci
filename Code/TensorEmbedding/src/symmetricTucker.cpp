// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "Tensor.hpp"

using namespace Rcpp;
using namespace arma;




// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
SEXP rcpparma_bothproducts(SEXP A_r, int n, int p, int k) {

    
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n,n, p, false);    


  SymmTensor tensorA (A);
  tensorA.setK(k);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

    return List::create(
        Named("A")=A,
        Named("L")=tensorA.L,
        Named("C")=tensorA.C,
        Named("gradL")=gradL,
        Named("gradC")=gradC
        );
}
