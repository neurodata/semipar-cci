// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include "Tensor.hpp"

#include "TensorSharedCore.hpp"

using namespace Rcpp;
using namespace arma;

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
SEXP symmetric_tensor_decomp(SEXP A_r, int n, int m, int k, int steps = 1000,
                             double delta1 = 1E-2, double delta2 = 1E-2,
                             double tol = 1E-8,
                             bool restrictCoreToDiag = true) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, m, false);

  bool logistic = true;
  SymmTensor tensorA(A, logistic, restrictCoreToDiag);
  tensorA.setK(k);

  tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

  cout << "finish" << endl;
  cout << tensorA.logisticLoss(tensorA.L, tensorA.C) << endl;

  return List::create(Named("A") = A, Named("L") = tensorA.L,
                      Named("C") = tensorA.C, Named("gradL") = gradL,
                      Named("gradC") = gradC);
}

// [[Rcpp::export]]
SEXP symm_group_tensor_decomp(SEXP A_r, SEXP group_r, int n, int m, int k,
                              int steps = 1000, double delta1 = 1E-2,
                              double delta2 = 1E-2, double tol = 1E-8,
                              bool restrictCoreToDiag = true) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, m, false);

  Rcpp::NumericVector gr(group_r);
  const vec g(gr.begin(), m, false);

  vec uniqueG = unique(g);

  int mG = uniqueG.n_elem;

  cout << " Found " << mG << " groups" << endl;

  cube groupAvgA = zeros<cube>(n, n, mG);
  vec groupCount = zeros<vec>(mG);

  for (int i = 0; i < mG; ++i) {
    int label = uniqueG(i);
    uvec pick = find(g == label);

    groupCount(i) = pick.n_elem;

    mat A_pick_sum = zeros(n,n);

    for (int j = 0; j <  groupCount(i); ++j)
    {
      A_pick_sum += A.slice(pick(j));
    }

    groupAvgA.slice(i) = A_pick_sum / (double)groupCount(i);
  }

  bool logistic = true;
  SymmTensor tensorA(groupAvgA, logistic, restrictCoreToDiag);
  tensorA.setK(k);
  tensorA.setW(groupCount);

  tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

  mat gradL = tensorA.gradL(tensorA.L, tensorA.C);
  cube gradC = tensorA.gradC(tensorA.L, tensorA.C);

  cout << "finish" << endl;
  cout << tensorA.logisticLoss(tensorA.L, tensorA.C) << endl;

  return List::create(Named("A") = A, Named("L") = tensorA.L,
                      Named("C") = tensorA.C, Named("gradL") = gradL,
                      Named("gradC") = gradC);
}

// SEXP symmetric_decomp_shared_core(SEXP A_r, int n, int p, int k,
//                                   int steps = 1000, double delta1 = 1E-2,
//                                   double delta2 = 1E-2, bool logistic = true,
//                                   double tol = 1E-8,
//                                   bool restrictCoreToDiag = true) {
//   Rcpp::NumericVector Ar(A_r);
//   const cube A(Ar.begin(), n, n, p, false);

//   SymmTensorSharedCore tensorA(A, logistic, restrictCoreToDiag);
//   tensorA.setK(k);

//   tensorA.OptConjugateGradient(steps, delta1, delta2, tol);

//   cube gradL = tensorA.gradL(tensorA.L, tensorA.C);
//   mat gradC = tensorA.gradC(tensorA.L, tensorA.C);

//   cout << "finish" << endl;
//   cout << tensorA.logisticLoss(tensorA.L, tensorA.C) << endl;

//   return List::create(Named("A") = A, Named("L") = tensorA.L,
//                       Named("C") = tensorA.C, Named("gradL") = gradL,
//                       Named("gradC") = gradC);
// }
