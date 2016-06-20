// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "Tensor.hpp"

using namespace Rcpp;
using namespace arma;

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
SEXP symmetric_tensor_decomp(SEXP A_r, int n, int p, int k, int steps = 1000,
                             double delta1 = 1E-2, double delta2 = 1E-2,
                             int loss_type = 1, double tol = 1E-8,
                             bool restrictCoreToDiag = true) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, p, false);

  LossType lossTypePick;
  if (loss_type == 1) {
    lossTypePick = logistic;
  }
  if (loss_type == 2) {
    lossTypePick = sqr;
  }

  SymmTensor tensorA(A, lossTypePick, restrictCoreToDiag);
  tensorA.setK(k);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

  tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

  return List::create(Named("A") = A, Named("L") = tensorA.L,
                      Named("C") = tensorA.C, Named("gradL") = gradL,
                      Named("gradC") = gradC);
}
