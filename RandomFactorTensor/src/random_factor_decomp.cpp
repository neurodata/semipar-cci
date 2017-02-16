// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <vector>
#include "RcppArmadillo.h"

#include "rfTensor.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP random_factor_decomp(SEXP A_r, int n, int m, SEXP g_r, int rank) {
  Rcpp::NumericVector Ar(A_r);
  const cube A(Ar.begin(), n, n, m, false);

  Rcpp::NumericVector gr(g_r);
  const vec g(gr.begin(), m, false);

  std::vector<cube> A_list;

  vec idx = linspace<vec>(0, g.n_elem - 1, g.n_elem);
  for (int i = 0; i <= max(g); ++i) {
    vec idx_g = idx(find(g == i));
    int idx0 = min(idx_g);
    int idx1 = max(idx_g);
    A_list.push_back(A.slices(idx0, idx1));
  }

  RFTensor(A_list, rank);

  return List::create(Named("A") = A);
}
