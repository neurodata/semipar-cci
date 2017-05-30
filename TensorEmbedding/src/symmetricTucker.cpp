// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "Tensor.hpp"
#include "TensorEM.hpp"
#include "TensorCovariate.hpp"

using namespace Rcpp;
using namespace arma;

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
SEXP symmetric_tensor_decomp(SEXP A_r, int n, int m, int k, int steps = 1000,
                             double delta1 = 1E-2, double delta2 = 1E-2,
                             double tol = 1E-8, int loss_type = 1,
                             bool restrictCoreToDiag = true) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, m, false);

  SymmTensor tensorA(A, loss_type, restrictCoreToDiag);
  tensorA.setK(k);

  tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

  cout << "finish" << endl;
  cout << tensorA.computeLoss(tensorA.L, tensorA.C) << endl;

  return List::create(Named("A") = A, Named("L") = tensorA.L,
                      Named("C") = tensorA.C, Named("gradL") = gradL,
                      Named("gradC") = gradC);
}



// [[Rcpp::export]]
SEXP symmetric_tensor_decomp_cov(SEXP A_r, SEXP X_r, int n, int m, int k, int steps = 1000, double lam =10,
                             double delta1 = 1E-2, double delta2 = 1E-2,
                             double tol = 1E-8, int loss_type = 1,
                             bool restrictCoreToDiag = true) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, m, false);

  Rcpp::NumericVector Xr(X_r);
  const mat X(Xr.begin(), n, k, false);

  SymmTensorCov tensorA(A,X,lam, loss_type, restrictCoreToDiag);
  tensorA.setK(k);

  tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

  cout << "finish" << endl;
  cout << tensorA.computeLoss(tensorA.L, tensorA.C) << endl;

  return List::create(Named("A") = A, Named("L") = tensorA.L,
                      Named("C") = tensorA.C, Named("gradL") = gradL,
                      Named("gradC") = gradC);
}
