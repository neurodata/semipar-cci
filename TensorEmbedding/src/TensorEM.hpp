using namespace arma;

// symmetric 3D tensor
// n x n x p
// latent factoring dimension k
class SymmTensorEM {
 public:
  int n, p, k;
  cube A;
  cube W;  // expectation of the Polya-Gamma

  // loss type
  // bool logistic;
  // bool poisson;

  int loss_type;  // 1: logistic 2: poisson
  bool coreDiag;  // restrict each layer of core to be diagonal

  // A =~ L C L
  mat L;   // left matrix n,k
  cube C;  // core tensor k,k,p
  vec w;   // weight vec p
  cube theta;
  double r;

  mat mapL;
  mat mapC;
  cube mapCcube;

  // direction for L and C to drop
  mat directionL;
  cube directionC;

  // variance for elements in L and C
  // regularization to ensure convexity
  double varLC;

  SymmTensorEM(cube _A, int _loss_type = 1, bool _coreDiag = false) {
    loss_type = _loss_type;

    coreDiag = _coreDiag;

    A = _A;
    n = A.n_rows;
    p = A.n_slices;

    theta = zeros(n, n, p);
    W = ones(n, n, p);

    r = 1000;  // aux parameter in Polya-Gamma

    w = ones(p);

    varLC = 1E4;
  }

  // set latent dimension
  void setK(int _k) {
    k = _k;

    mat Amean = sum(A, 2) / p;

    vec s;
    mat V;

    if (loss_type == 1) svd(L, s, V, Amean);
    if (loss_type == 2) svd(L, s, V, log(Amean + 0.001));

    // L = trans(chol(Amean));

    L.resize(n, k);
    s.resize(k);

    C = zeros(k, k, p);
    for (int i = 0; i < p; ++i) {
      C.slice(i) = diagmat(s);  // making C symmetric
    }
  }

  // set weight for each 2-mode array

  void setW(vec _w) { w = _w; }

  // loss functions and its derivative wrt LCL:

  double expectedLoss(mat L, cube C) {
    // vec C_flat = reshape(C, k * k * p, 1, 1);

    // if (any(C_flat < 0)) return INFINITY;

    // cube theta = zeros(n, n, p);
    cube loss = zeros(n, n, p);

#pragma omp parallel for
    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      mat W_local = W.slice(i);

      theta.slice(i) = L * C_local * L.t();

      // loss.slice(i) =
      // -theta.slice(i) % A.slice(i) + log(1 + exp(theta.slice(i)));

      if (loss_type == 1) {
        loss.slice(i) = -theta.slice(i) % (A.slice(i) - 0.5) +
                        W_local % theta.slice(i) % theta.slice(i) / 2.0;
      }
      if (loss_type == 2) {
        double logr = log(r);
        loss.slice(i) =
            -theta.slice(i) % (A.slice(i) - 0.5 * r) +
            W_local % (theta.slice(i) - logr) % (theta.slice(i) - logr) / 2.0;
      }

      loss.slice(i) *= w(i);

      // loss.slice(i) = trimatu(loss.slice(i));
      // loss.slice(i).diag().fill(0);
    }

    return accu(loss) + accu(C % C) / 2 / varLC + accu(L % L) / 2 / varLC;
  }

  cube deriExpectedLoss(mat L, cube C) {
    cube theta = zeros(n, n, p);
    cube deriLoss = zeros(n, n, p);

#pragma omp parallel for
    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      mat W_local = W.slice(i);

      theta.slice(i) = L * C_local * L.t();

      // deriLoss.slice(i) = -A.slice(i) + 1.0 / (1.0 + exp(-theta.slice(i)));
      if (loss_type == 1) {
        deriLoss.slice(i) = -(A.slice(i) - 0.5) + W_local % theta.slice(i);
      }

      if (loss_type == 2) {
        double logr = log(r);
        deriLoss.slice(i) =
            -(A.slice(i) - 0.5 * r) + W_local % (theta.slice(i) - logr);
      }
      // deriLoss.slice(i) = trimatu(deriLoss.slice(i));
      // deriLoss.slice(i).diag().fill(0);

      deriLoss.slice(i) *= w(i);
    }

    return deriLoss;
  }

  // gradient for L
  mat gradL(mat L, cube C) {
    cube grad(n, k, p);

    cube grad_LCL;
    grad_LCL = deriExpectedLoss(L, C);
// else
// grad_LCL = deriSqrLoss(L, C);

#pragma omp parallel for
    for (int i = 0; i < p; ++i) {
      mat LC = L * C.slice(i);
      grad.slice(i) = 2 * grad_LCL.slice(i) * LC;
    }
    mat loglikDeri = sum(grad, 2);

    return loglikDeri + L / varLC;  // sum over last dim
  }

  // gradient for C

  cube gradC(mat L, cube C) {
    cube grad = zeros(k, k, p);

    cube grad_LCL;
    grad_LCL = deriExpectedLoss(L, C);

#pragma omp parallel for
    for (int i = 0; i < p; ++i) {
      mat LC = L * C.slice(i);

      if (coreDiag) {
        for (int a = 0; a < k; ++a) {
          mat deriC = zeros(k, k);
          deriC(a, a) = 1;
          grad(a, a, i) = accu(grad_LCL.slice(i) % (L * deriC * L.t()));
        }
      } else {
        for (int a = 0; a < k; ++a) {
          for (int b = 0; b <= a; ++b) {
            mat deriC = zeros(k, k);
            deriC(a, b) = 1;
            deriC(b, a) = 1;
            grad(a, b, i) = accu(grad_LCL.slice(i) % (L * deriC * L.t()));
            grad(b, a, i) = grad(a, b, i);
          }
        }
      }
    }

    return grad + C / varLC;
  }

  // compute loss when move delta
  double lossAtDelta(double delta, int paramIdx = 0) {
    double cur_loss = expectedLoss(L, C);
    double result = 0;
    if (paramIdx == 0) {
      mat prop = L - delta * directionL;
      result = expectedLoss(prop, C);
    }
    if (paramIdx == 1) {
      cube prop = C - delta * directionC;
      result = expectedLoss(L, prop);
    }
    return result;
  }

  double lineSearch(double l, double r, double loss_l, double loss_r,
                    int paramIdx = 0) {
    if (fabs(loss_l - loss_r) < 1E-3) return l;

    double m = (l + r) / 2;

    double loss_m = lossAtDelta(m, paramIdx);

    if (loss_l < loss_m & loss_l < loss_r) return l;

    if (loss_r < loss_m & loss_r < loss_l) return r;

    if (loss_m < loss_l & loss_m < loss_r) {
      double opt_l = lineSearch(l, m, loss_l, loss_m, paramIdx);
      double opt_r = lineSearch(m, r, loss_m, loss_r, paramIdx);

      double loss_opt_l = lossAtDelta(opt_l, paramIdx);
      double loss_opt_r = lossAtDelta(opt_r, paramIdx);
      if (loss_opt_l < loss_opt_r)
        return opt_l;
      else
        return opt_r;
    }
    return l;
  }

  void expectationW() {
    if (loss_type == 1) W = 0.5 / abs(theta) % tanh(abs(theta) / 2.0);
    if (loss_type == 2)
      W = 0.5 * r / abs(theta - log(r)) % tanh(abs(theta - log(r)) / 2.0);
  }

  void OptConjugateGradient(int steps = 100, double delta1 = 1E-3,
                            double delta2 = 1E-3, double tol = 1E-8) {
    mat gradientL = gradL(L, C);
    cube gradientC = gradC(L, C);

    directionL = gradientL;
    directionC = gradientC;

    double pre_loss = INFINITY;

    double cur_loss = expectedLoss(L, C);
    expectationW();

    double best_loss = INFINITY;

    // double delta1_decre = log(1E-8 / delta1) / (double)steps;
    // double delta2_decre = log(1E-8 / delta2) / (double)steps;

    for (int i = 0; i < steps; ++i) {
      {
        R_CheckUserInterrupt();
        mat gradient0 = gradientL;
        gradientL = gradL(L, C);
        // Polak–Ribière conjugate gradient
        double beta = accu(gradientL % (gradientL - gradient0)) /
                      accu(gradient0 % gradient0);

        if (beta < 0) beta = 0;

        if (!std::isnan(beta)) directionL = gradientL + beta * directionL;

        double cur_loss = expectedLoss(L, C);
        double loss_delta = lossAtDelta(delta1, 0);
        double opt_delta = lineSearch(0, delta1, cur_loss, loss_delta, 0);
        L -= opt_delta * directionL;
      }

      double cur_loss = expectedLoss(L, C);
      expectationW();

      {
        cube gradient0 = gradientC;
        gradientC = gradC(L, C);
        double beta = accu(gradientC % gradientC) / accu(gradient0 % gradient0);
        if (!std::isnan(beta)) directionC = gradientC + beta * directionC;

        double cur_loss = expectedLoss(L, C);
        double loss_delta = lossAtDelta(delta2, 1);
        double opt_delta = lineSearch(0, delta2, cur_loss, loss_delta, 1);
        C -= opt_delta * directionC;
      }

      cur_loss = expectedLoss(L, C);
      expectationW();
      // if ((abs(directionL)).max() + (abs(directionC)).max() < 1E-5)
      if (cur_loss < best_loss) {
        mapL = L;
        mapC = extractCdiag();
        mapCcube = C;
        best_loss = cur_loss;
      }

      // delta1 *= exp(delta1_decre);
      // delta2 *= exp(delta2_decre);

      if ((abs(directionL)).max() + (abs(directionC)).max() < tol) {
        // if (fabs((pre_loss - cur_loss) / cur_loss) < tol) {
        // if (delta1 < 1E-8 & delta2 < 1E-8)
        break;
        // else {
        // delta1 *= 0.1;
        // delta2 *= 0.1;
        // }
      }

      pre_loss = cur_loss;

      cout << "step " << i << ", current loglik " << cur_loss
           << ", best loglik: " << best_loss << endl;
      cout << delta1 << " " << delta2 << endl;
    }
  }

  mat extractCdiag() {
    mat Cdiag = zeros<mat>(p, k);
    for (int i = 0; i < p; ++i) {
      mat C_local = C.slice(i);
      Cdiag.row(i) = trans(C_local.diag());
    }

    return Cdiag;
  }

  // end of class definition
};
